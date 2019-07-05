#! /usr/bin/env python
# -*- coding: utf-8 -*-

import json
import sys
import os
import collections
import random
import argparse
import dendropy
from dendropy.model import protractedspeciation
import datetime
import socket
from spdw import spdwlib

def _log(msg):
    sys.stderr.write("- {}\n".format(msg))

### Example invocation:
# python ../../../bin/spdw-build-delineate-jobs.py -t test1 -n 1 --num-extant-lineages 20 joint-partition-probabilities --constrain-partitions topological --max-unconstrained-leaves 5

template = """\
#! /bin/bash

set -e -o pipefail

{jobs}

echo "Job $JOBID completed on $(date --rfc-3339='seconds')"
"""
species_partition_estimation_job_template="""\
delineate-estimate-species-partition.py {underflow_protection} -c {run_config_filepath} -t {tree_filepath} -I {speciation_completion_rate} -l num_lineages:{num_lineages} -l num_species:{num_species} -l src_filepath:{tree_filepath} -l true_speciation_completion_rate:{true_speciation_completion_rate} -l batch_id:{batch_id} > {delineate_results_filepath}
"""
species_partition_estimation_joint_probability_analysis_template="""\
{post_analysis_performance_assessment_command} {run_config_filepath} {delineate_results_filepath} > {joint_performance_assessment_results_filepath}
"""
species_partition_estimation_marginal_probability_analysis_template="""\
{post_analysis_performance_assessment_command} --marginal -t {tree_filepath} {run_config_filepath} {delineate_results_filepath} > {marginal_performance_assessment_results_filepath}
"""

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("population_tree_filepath",
            help="Path to population tree")
    parser.add_argument("-f", "--tree-format",
            default="nexus",
            dest="schema",
            choices=["nexus", "newick"],
            help="Input trees format [default: $(default)s].")
    parser.add_argument("-o", "--output-prefix",
            default="bpprun",
            help="Run title (default: '%(default)s')")
    parser.add_argument("-z", "--random-seed",
            type=int,
            default=None,
            help="Seed for random number generator engine.")
    parser.add_argument("-u", "--underflow-protection",
            action="store_true",
            default=False,
            help="Try to protect against underflow by using special number handling classes (slow).",)

    partition_estimation_test_group = parser.add_argument_group("Partition Estimation Test Settings")
    partition_estimation_test_group.add_argument("--constrain-partitions",
            dest="constrain_partitions",
            choices=["random", "topological", "user"],
            default="random",
            help="""
            Constrain partition sets by specifying the (true) species assignments of some lineages.
            Options are:
                'random': randomly select lineages from leaf set;
                'topological': select random internal node in tree to suppress (species assignments
                of all leaves not descending from this node will be known).
            """)
    partition_estimation_test_group.add_argument("--num-unconstrained-leaves",
            default=None,
            type=int,
            help="Exact number of leaves with unknown species assignments (overrides min/max below).")
    partition_estimation_test_group.add_argument("--min-unconstrained-leaves",
            default=None,
            type=int,
            help="Minimum number of leaves with unknown species assignments.")
    partition_estimation_test_group.add_argument("--max-unconstrained-leaves",
            default=None,
            type=int,
            help="Maximum number of leaves with unknown species assignments.")
    partition_estimation_test_group.add_argument("--specify-true-speciation-completion-rate",
            action="store_true",
            default=False,
            help="True speciation completion rate will be provided to partition probability calculator.")
    args = parser.parse_args()
    if args.random_seed is None:
        random_seed = random.randint(0, sys.maxsize-1)
    else:
        random_seed = args.random_seed

    rng = random.Random(random_seed)
    _log("Random seed: {}".format(random_seed))

    lineage_tree = dendropy.Tree.get(
            path=args.population_tree_filepath,
            schema=args.schema)

    if args.constrain_partitions is "user":
        raise NotImplementedError()
    else:
        species_lineage_label_map, lineage_species_label_map = spdwlib.build_tree_label_maps(lineage_tree=lineage_tree)

        # ensure if parent is collapsed, all children are collapsed
        for gnd in lineage_tree.preorder_node_iter():
            if gnd.parent_node is not None and gnd.parent_node.annotations["is_collapsed"].value:
                gnd.annotations["is_collapsed"] = True

        # identify the lowest nodes (closest to the tips) that are open, and
        # add its children if the children are (a) leaves; or (b) themselves
        # are closed --- indicating, essentially a tip
        terminal_clades = {}
        for nd in lineage_tree.postorder_node_iter():
            if nd.annotations["is_collapsed"].value:
                continue
            for child in nd.child_node_iter():
                if child.is_leaf():
                    terminal_clades[child] = set([child])
                elif child.annotations["is_collapsed"].value:
                    terminal_clades[child] = set([desc for desc in child.leaf_iter()])

        for nd in lineage_tree:
            if nd not in terminal_clades:
                nd.annotations["population"] = "0"
                continue
            if nd.is_leaf():
                nd.annotations["population"] = nd.taxon.label
            else:
                nd.annotations["population"] = "+".join(desc.taxon.label for desc in nd.leaf_iter())
            # print("{}: {}".format(nd.annotations["population"], len(terminal_clades[nd])))

        ## Come up with species labels that take into account:
        ## - a terminal clade may consist of a union of (true) species: S1, S2, etc.
        ## - the same species may be associated with different terminal clades
        ## - the same species may be mixed up with other species in different
        ## - al clades (e.g., terminal clade A may include lineages associated
        #         with species 1 and 2, while terminal clade B may include
        #         associated with species 1 and 3. Here we come up with a nominal
        #         species = S1+S2+S3 if S1, S2, and S3 are not found associated with
        ##        any other species in any other terminal clade.
        ## - a "species complex" are thus the most expansive set that does not exclude
        ##   any pair of nominal species that occur together in any terminal clade
        ## If each terminal clade consists of lineages of multiple nominal species;
        terminal_clade_found_species = {}
        species_complexes = {}
        for terminal_clade in terminal_clades:
            # print("---\n{}:".format(terminal_clade.annotations["population"]))
            found_species_labels = set()
            for lineage in terminal_clades[terminal_clade]:
                if not lineage.is_leaf():
                    continue
                # print(lineage.taxon.label)
                species_label, lineage_label = spdwlib.ProtractedSpeciationTreeGenerator.decompose_species_lineage_label(lineage.taxon.label)
                found_species_labels.add(species_label)
            # print(found_species_labels)
            for species_label in found_species_labels:
                if species_label in species_complexes:
                    for ex2 in list(species_complexes[species_label]):
                        species_complexes[ex2].update(found_species_labels)
                    # species_complexes[species_label].update(found_species_labels)
                else:
                    species_complexes[species_label] = set(found_species_labels)
            terminal_clade_found_species[terminal_clade] = found_species_labels
        # print(species_complexes)

        terminal_clade_species_names = {}
        for nd in terminal_clades:
            x1 = None
            for spp_label in terminal_clade_found_species[nd]:
                if x1 is None:
                    x1 = species_complexes[spp_label]
                else:
                    assert x1 == species_complexes[spp_label]
            terminal_clade_species_names[nd] = "+".join(x1)

        species_names = sorted(list(terminal_clade_species_names.values()))

        for nd in lineage_tree:
            if nd not in terminal_clades:
                nd.annotations["population"] = "0"
                continue
            if nd.is_leaf():
                nd.annotations["population"] = nd.taxon.label
            else:
                nd.annotations["population"] = "+".join(desc.taxon.label for desc in nd.leaf_iter())
            nd.annotations["species"] = terminal_clade_species_names[nd]

        while True:
            selected_lineages = set()
            unconstrained_tip_lineages = set()
            lineage_pool = list(terminal_clades.keys())
            if args.min_unconstrained_leaves is None and args.num_unconstrained_leaves is None:
                args.min_unconstrained_leaves = int(len(lineage_tree.taxon_namespace) / 2)
            if args.num_unconstrained_leaves is not None:
                condition_fn = lambda: len(unconstrained_tip_lineages) == args.num_unconstrained_leaves
            elif args.min_unconstrained_leaves is not None and args.max_unconstrained_leaves is not None:
                condition_fn = lambda: len(unconstrained_tip_lineages) >= args.min_unconstrained_leaves and len(unconstrained_tip_lineages) <= args.max_unconstrained_leaves
            elif args.min_unconstrained_leaves is not None:
                condition_fn = lambda: len(unconstrained_tip_lineages) >= args.min_unconstrained_leaves
            else:
                sys.exit("Need to specify '--min-unconstrained-leaves' or '--num-unconstrained-leaves'")
            while lineage_pool and not condition_fn():
                node = lineage_pool.pop(rng.randint(0, len(lineage_pool)-1))
                unconstrained_tip_lineages.update(terminal_clades[node])
                selected_lineages.add(node)
            if condition_fn:
                break
            else:
                print("Failed to meet condition")
        msg = ["Open terminal lineages:"]
        for nd in terminal_clades:
            msg.append("  - [{}] {} [SPP='{}']".format(
                "UNCONSTRAINED" if nd in selected_lineages else " CONSTRAINED ",
                nd.annotations["population"].value,
                terminal_clade_species_names[nd]))
        _log("\n".join(msg))
        sys.exit(0)
        # species_leafset_constraints, constrained_lineage_leaf_labels, unconstrained_lineage_leaf_labels, species_leafset_constraint_label_map = spdwlib.generate_constraints(
        #         lineage_tree=lineage_tree,
        #         orthospecies_tree=None,
        #         constraint_type=args.constrain_partitions,
        #         species_lineage_label_map=species_lineage_label_map,
        #         lineage_species_label_map=lineage_species_label_map,
        #         min_unconstrained_leaves=args.min_unconstrained_leaves,
        #         max_unconstrained_leaves=args.max_unconstrained_leaves,
        #         num_unconstrained_leaves=args.num_unconstrained_leaves,
        #         rng=rng,
        #         )
        # print(f"constrained: {constrained_lineage_leaf_labels}")
        # print(f"unconstrained: {unconstrained_lineage_leaf_labels}")
    assert args.constrain_partitions is not None
    config = {}
    config["species_leafset_constraints"] = species_leafset_constraints
    spdwlib.decorate_lineage_tree(
            lineage_tree=lineage_tree,
            lineage_species_label_map=lineage_species_label_map,
            unconstrained_lineage_leaf_labels=unconstrained_lineage_leaf_labels,
            )
    lineage_tree.write(path="{}.lineages.nex".format(args.output_prefix),
            supplemental_blocks=[lineage_tree.figtree_block],
            schema="nexus")

if __name__ == '__main__':
    main()

