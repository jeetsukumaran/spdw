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
        candidate_nodes = set()
        for nd in lineage_tree.postorder_node_iter():
            if nd.annotations["is_collapsed"].value:
                continue
            for child in nd.child_node_iter():
                if child.is_leaf() or child.annotations["is_collapsed"].value:
                    candidate_nodes.add(child)
        for nd in lineage_tree:
            if nd not in candidate_nodes:
                nd.annotations["tag"] = "0"
                continue
            if nd.is_leaf():
                nd.annotations["tag"] = nd.taxon.label
            else:
                nd.annotations["tag"] = "+".join(desc.taxon.label for desc in nd.leaf_iter())
            print("Identified: {}".format(nd.annotations["tag"].value))
        lineage_tree.write(path="x.tre", schema="nexus")
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

