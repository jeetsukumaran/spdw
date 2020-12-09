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
import itertools

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
            default="delineate-run",
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

    terminal_population_clades, lineage_population_clade_map = spdwlib.find_terminal_population_clades(lineage_tree)
    terminal_population_clade_species_identities = spdwlib.identify_terminal_population_clade_species(terminal_population_clades)
    if args.constrain_partitions == "user":
        raise NotImplementedError()
    else:
        constraints = spdwlib.generate_constraints_from_collapsed_tree(
                terminal_population_clades=terminal_population_clades,
                min_unconstrained_leaves=args.min_unconstrained_leaves,
                max_unconstrained_leaves=args.max_unconstrained_leaves,
                num_unconstrained_leaves=args.num_unconstrained_leaves,
                rng=rng,
                )
        msg = spdwlib.format_population_clade_constraint_report(
                constraints=constraints,
                terminal_population_clades=terminal_population_clades,
                terminal_population_clade_species_identities=terminal_population_clade_species_identities,
                )
        _log(msg)
    msg = spdwlib.format_lineage_constraint_report(
                constraints=constraints,
                lineage_population_clade_map=lineage_population_clade_map,
                terminal_population_clades=terminal_population_clades,
                terminal_population_clade_species_identities=terminal_population_clade_species_identities,
            )
    _log(msg)

    config_table = []
    for lineage_nd in itertools.chain(constraints["constrained_lineages"], constraints["unconstrained_lineages"]):
        population_nd = lineage_population_clade_map[lineage_nd]
        row = {
                "lineage": lineage_nd.taxon.label,
                "species": terminal_population_clade_species_identities[population_nd],
                "status": "1" if lineage_nd in constraints["constrained_lineages"] else "0"
        }
        config_table.append(row)
    out_fpath = "{}.tsv".format(args.output_prefix)
    delimiter = "\t"
    header_row = list(config_table[0].keys())
    with open(out_fpath, "w") as out:
        out.write("{}\n".format(delimiter.join(header_row)))
        for row in config_table:
            out.write("{}\n".format(delimiter.join(row[k] for k in header_row)))
    lineage_tree.write(
            path="{}.tree.nex".format(args.output_prefix),
            schema="nexus")

if __name__ == '__main__':
    main()

