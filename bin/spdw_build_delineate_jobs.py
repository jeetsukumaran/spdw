#! /usr/bin/env python
# -*- coding: utf-8 -*-

import json
import sys
import os
import collections
import random
import argparse
# from generate_trees import generate_source_trees_r
from dendropy.model import protractedspeciation
import datetime
import socket

SCRIPT_DIR = os.path.dirname(__file__)
LOGPATH = os.path.expandvars(os.path.expanduser(os.path.join("~", ".delineate-performance.log")))
START_DATETIME = datetime.datetime.now()

### Example invocation:
# python ../../../bin/spdw_build_delineate_jobs.py -t test1 -n 1 --num-extant-lineages 20 joint-partition-probabilities --constrain-partitions topological --max-unconstrained-leaves 5

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

class ColorMap(object):

    contrast_pairs = [
            ["#ffc20a", "#0c7bdc"],
            ["#1aff1a", "#4b0092"],
            ["#994f00", "#006cd1"],
            ["#fefe62", "#d35fb7"],
            ["#e1be6a", "#40b0a6"],
            ["#005ab5", "#dc3220"],
            ["#e66100", "#5d3a9b"],
            ["#1a85ff", "#d41159"],
    ]

    def __init__(self):
        self.label_color_map = {}
        self.color_label_map = {}
        self.colors = [
            "#000000",
            "#e69f00",
            "#56b4e9",
            "#009e73",
            "#f0e442",
            "#0072b2",
            "#d55e00",
            "#cc79a7"
        ]
        self.available_colors = list(self.colors)

    def __call__(self, label):
        try:
            return self.label_color_map[label]
        except KeyError:
            new_color = self.available_colors.pop()
            self.label_color_map[label] = new_color
            self.color_label_map[new_color] = label
        return new_color

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--title",
            action="store",
            default="delineate-test",
            help="Title of run [default=%(default)s].")
    parser.add_argument("--splitting-rate",
            action="store",
            type=float,
            default=0.10,
            help="Rate of origin of new lineages (population isolation rate) [default=%(default)s].")
    parser.add_argument("--speciation-completion-rate",
            action="store",
            type=float,
            default=0.01,
            help="Rate at which lineage develop into independent species [default=%(default)s].")
    parser.add_argument("-n", "--num-replicates",
            type=int,
            default=30,
            help="Number of replicates per combination of parameters (default=%(default)s).")
    parser.add_argument("-z", "--random-seed",
            default=None,
            help="Random seed.")
    parser.add_argument("-u", "--underflow-protection",
            action="store_true",
            default=False,
            help="Try to protect against underflow by using special number handling classes (slow).",)
    parser.add_argument("--clean",
            action="store_true",
            default=False,
            help="Clean up run *data* files after post-run analysis (ONLY summaries, job, and standard output will be kept).")
    parser.add_argument("--very-clean",
            action="store_true",
            default=False,
            help="Clean up run job and *data* files after post-run analysis (ONLY summaries and standard output will be kept).")
    parser.add_argument("--write-extra-for-demo",
            action="store_true",
            default=False,
            help="Write extra files for demonstration purposes.")

    regime_group = parser.add_argument_group("Regime")
    regime_group.add_argument("--max-time",
            default=None,
            type=float,
            help="Source trees generated with this crown age.")
    regime_group.add_argument("--num-extant-lineages",
            default=None,
            type=int,
            help="Source trees generated with exactly this number of tip lineages (incipient species + orthospecies).")
    regime_group.add_argument("--min-extant-lineages",
            default=None,
            type=int,
            help="Source trees generated with at least this number of tip lineages (incipient species + orthospecies).")
    regime_group.add_argument("--num-extant-orthospecies",
            default=None,
            type=int,
            help="Source trees generated with this number of orthospecies ('good' or true species).")
    regime_group.add_argument("--min-extant-orthospecies",
            default=2,
            type=int,
            help="Reject source trees with less than this number of orthospecies ('good' or true species).")
    partition_estimation_test_group = parser.add_argument_group("Partition Estimation Test Settings")
    partition_estimation_test_group.add_argument("--constrain-partitions",
            dest="constrain_partitions",
            choices=["random", "topological"],
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
    command_kwargs = {}

    selected_condition = None
    for kw in ("max_time", "num_extant_lineages", "num_extant_orthospecies"):
        if getattr(args, kw) is not None:
            if selected_condition:
                sys.exit("Need to specify only one of: '--max-time', '--num-extant-lineages', '--num-extant-orthospecies'")
            selected_condition = kw
            command_kwargs[kw] = getattr(args, kw)
    if selected_condition is None:
        sys.exit("Need to specify at least one of: '--max-time', '--num-extant-lineages', '--num-extant-orthospecies'")
    if args.random_seed is None:
        args.random_seed = random.randrange(sys.maxsize)
    rng = random.Random(args.random_seed)
    batch_id = "::".join([
            socket.gethostname(),
            START_DATETIME.strftime("%Y%m%d%H%M%S"),
            str(args.random_seed),
            ])
    log_message = [
        batch_id,
        START_DATETIME.strftime("%Y-%m-%d"),
        START_DATETIME.strftime("%H:%M:%S"),
        "'" + os.path.abspath(os.getcwd()) + "'",
        "'" + " ".join(sys.argv) + "'",
        str(args.random_seed),
    ]
    with open(LOGPATH, "a") as dest:
        dest.write("\t".join(log_message) + "\n")
    color_map = ColorMap()

    true_speciation_completion_rate = args.speciation_completion_rate
    splitting_rate = args.splitting_rate
    extinction_rate = 0.0
    data = collections.OrderedDict()
    data["params"] = collections.OrderedDict()
    data["params"]["good_species_speciation_initiation_rate"] = splitting_rate
    data["params"]["true_speciation_completion_rate"] = true_speciation_completion_rate
    data["params"]["incipient_species_speciation_initiation_rate"] = splitting_rate
    data["params"]["good_species_extinction_rate"] = extinction_rate
    data["params"]["incipient_species_extinction_rate"] = extinction_rate
    data["condition"] = selected_condition
    data["condition_value"] = command_kwargs[selected_condition]
    data["trees"] = []
    psm = protractedspeciation.ProtractedSpeciationProcess(
        speciation_initiation_from_orthospecies_rate=splitting_rate,
        speciation_initiation_from_incipient_species_rate=splitting_rate,
        speciation_completion_rate=true_speciation_completion_rate,
        orthospecies_extinction_rate=extinction_rate,
        incipient_species_extinction_rate=extinction_rate,
        rng=rng)
    output_prefix = "{}_spr{:0.3f}_".format(args.title, true_speciation_completion_rate)
    tree_idx = 0
    for tree_idx in range(args.num_replicates):
        while True:
            # make sure that the tree we generate has enough species
            lineage_tree, orthospecies_tree = psm.generate_sample(**command_kwargs)
            if len(orthospecies_tree.taxon_namespace) >= args.min_extant_orthospecies:
                ok = []
                if args.min_unconstrained_leaves:
                    if len(lineage_tree.taxon_namespace) >= args.min_unconstrained_leaves:
                        ok.append(True)
                    else:
                        ok.append(False)
                if args.min_extant_lineages:
                    if len(lineage_tree.taxon_namespace) >= args.min_extant_lineages:
                        ok.append(True)
                    else:
                        ok.append(False)
                if all(ok):
                    break
        sorted_species_lineages_map = collections.OrderedDict()
        lineage_label_species_label_map = {}
        for k in sorted([t.label for t in orthospecies_tree.taxon_namespace]):
            sorted_species_lineages_map[k] = []
        for ond in orthospecies_tree.leaf_node_iter():
            sorted_species_lineages_map[ond.taxon.label] = sorted([lnd.taxon.label for lnd in ond.lineage_tree_nodes])
            for lnd in ond.lineage_tree_nodes:
                lineage_label_species_label_map[lnd.taxon.label] = ond.taxon.label
        true_species_leafsets = sorted(sorted_species_lineages_map.values())
        entry = collections.OrderedDict()
        entry["tree_filepath"] = "{}.{:04d}.nex".format(output_prefix, tree_idx+1)
        entry["run_config_filepath"] = "{}.{:04d}.json".format(output_prefix, tree_idx+1)
        entry["lineage_taxon_namespace"] = [t.label for t in lineage_tree.taxon_namespace]
        entry["lineage_tree"] = lineage_tree.as_string("newick").replace("\n", "")
        entry["species_taxon_namespace"] = sorted(sorted_species_lineages_map.keys())
        entry["species_tree"] = orthospecies_tree.as_string("newick").replace("\n", "")
        entry["species_lineages_map"] = sorted_species_lineages_map
        data["trees"].append(entry)
        lineage_tree.write(path=entry["tree_filepath"], schema="nexus")

        config = collections.OrderedDict()
        # for speciation rate estimation; ignored by species partition estimation
        # for species partition estimation
        species_leafset_constraints = None
        if args.constrain_partitions is not None:
            if args.constrain_partitions == "topological":
                lineage_tree_internal_nodes = [lnd for lnd in lineage_tree.postorder_internal_node_iter() if lnd is not lineage_tree.seed_node]
                rng.shuffle(lineage_tree_internal_nodes)
                true_unconstrained_lineage_leaf_labels = None
                for lineage_tree_internal_node in lineage_tree_internal_nodes:
                    true_unconstrained_lineage_leaf_labels = set([lineage_tree_leaf_node.taxon.label for lineage_tree_leaf_node in lineage_tree_internal_node.leaf_iter()])
                    if not args.num_unconstrained_leaves and not args.min_unconstrained_leaves and not args.max_unconstrained_leaves:
                        break
                    elif args.num_unconstrained_leaves:
                        if len(true_unconstrained_lineage_leaf_labels) == args.num_unconstrained_leaves:
                            break
                    elif args.min_unconstrained_leaves and args.max_unconstrained_leaves:
                        if len(true_unconstrained_lineage_leaf_labels) >= args.min_unconstrained_leaves and len(true_unconstrained_lineage_leaf_labels) <= args.max_unconstrained_leaves:
                            break
                    elif args.min_unconstrained_leaves and len(true_unconstrained_lineage_leaf_labels) >= args.min_unconstrained_leaves:
                        break
                    elif args.max_unconstrained_leaves and len(true_unconstrained_lineage_leaf_labels) <= args.max_unconstrained_leaves:
                        break
                    true_unconstrained_lineage_leaf_labels = None
                else:
                    raise ValueError("Unable to meet min/max unconstrained leaves criteria.")
            elif args.constrain_partitions == "random":
                lineage_leaf_labels = [taxon.label for taxon in lineage_tree.taxon_namespace]
                if args.num_unconstrained_leaves:
                    num_to_sample = args.num_unconstrained_leaves
                else:
                    if args.min_unconstrained_leaves:
                        min_count = args.min_unconstrained_leaves
                    else:
                        min_count = 1
                    if args.max_unconstrained_leaves:
                        max_count = args.max_unconstrained_leaves
                    else:
                        max_count = len(lineage_leaf_labels)
                    num_to_sample = rng.randint(min_count, max_count)
                true_unconstrained_lineage_leaf_labels = rng.sample(lineage_leaf_labels, num_to_sample)
            true_constrained_lineage_leaf_labels = sorted([lineage_tree_leaf_node.taxon.label for lineage_tree_leaf_node in lineage_tree.leaf_node_iter() if lineage_tree_leaf_node.taxon.label not in true_unconstrained_lineage_leaf_labels])
            species_leafset_constraint_label_map = collections.OrderedDict()
            for lineage_leaf_label in true_constrained_lineage_leaf_labels:
                true_species_label = lineage_label_species_label_map[lineage_leaf_label]
                try:
                    species_leafset_constraint_label_map[true_species_label].append(lineage_leaf_label)
                except KeyError:
                    species_leafset_constraint_label_map[true_species_label] = [lineage_leaf_label]
            for species_label in species_leafset_constraint_label_map:
                species_leafset_constraint_label_map[species_label].sort()
            species_leafset_constraints = []
            for sp in sorted_species_lineages_map:
                if sp in species_leafset_constraint_label_map:
                    species_leafset_constraints.append(sorted(species_leafset_constraint_label_map[sp]))
        else:
            assert args.constrain_partitions is None
            assert args.num_unconstrained_leaves is None
            assert args.min_unconstrained_leaves is None
            assert args.max_unconstrained_leaves is None
            true_constrained_lineage_leaf_labels = []
            true_unconstrained_lineage_leaf_labels = [lineage_tree_leaf_node.taxon.label for lineage_tree_leaf_node in lineage_tree.leaf_node_iter()]
            species_leafset_constraint_label_map = {}
            species_leafset_constraints = None
        if species_leafset_constraints is not None:
            # this is actually used by the DELINEATE program
            assert args.constrain_partitions is not None
            config["species_leafset_constraints"] = species_leafset_constraints
        else:
            assert args.constrain_partitions is None
            try:
                del config["species_leafset_constraints"]
            except KeyError:
                pass

        # extra files for demo
        if args.write_extra_for_demo:
            constrained_color, unconstrained_color = ColorMap.contrast_pairs[0]
            demo_output_prefix = "{}.{:04d}.demo".format(output_prefix, tree_idx+1)
            true_unconstrained_lineage_leaf_label_set = set(true_unconstrained_lineage_leaf_labels)
            for nd in lineage_tree:
                if nd.is_leaf():
                    is_constrained = nd.taxon.label in true_unconstrained_lineage_leaf_label_set
                    nd.annotations["constrained"] = is_constrained
                    species_label = lineage_label_species_label_map[nd.taxon.label]
                    nd.annotations["species"] = species_label
                    if False: #args.color_by_species:
                        nd.annotations["!color"] = color_map(species_label)
                    else:
                        if is_constrained:
                            nd.annotations["!color"] = constrained_color
                        else:
                            nd.annotations["!color"] = unconstrained_color
                else:
                    nd.annotations["!color"] = "#aaaaaa"
            lineage_tree_figtree_block = [
                    'set appearance.branchLineWidth=5.0;',
                    'set scaleBar.isShown=false;',
            ]
            # lineage_tree_figtree_block.extend([
            #         'set nodeShapeExternal.colourAttribute="constrained";',
            #         'set nodeShapeExternal.isShown=true;',
            #         'set nodeShapeExternal.minSize=10.0;',
            #         'set nodeShapeExternal.scaleType=Width;',
            #         'set nodeShapeExternal.shapeType=Circle;',
            #         'set nodeShapeExternal.size=10.0;',
            #         'set nodeShapeExternal.sizeAttribute="Fixed";',
            #     ])
            lineage_tree_figtree_block.extend([
                'set tipLabels.colorAttribute="species";',
                'set tipLabels.displayAttribute="species";',
                'set tipLabels.fontName="sansserif";',
                'set tipLabels.fontSize=30;',
                'set tipLabels.fontStyle=0;',
                'set tipLabels.isShown=true;',
                'set tipLabels.significantDigits=4;',
                ])
            # lineage_tree_figtree_block.extend([
                # 'set legend.attribute="constrained";',
                # 'set legend.fontSize=10.0;',
                # 'set legend.isShown=true;',
                # 'set legend.significantDigits=4;',
            #     ])
            lineage_tree_figtree_block = "begin figtree;\n{}\nend;\n".format("\n".join(lineage_tree_figtree_block))
            lineage_tree.write(path="{}.lineages.nex".format(demo_output_prefix),
                    supplemental_blocks=[lineage_tree_figtree_block],
                    schema="nexus")

            orthospecies_tree_figtree_block = [
                    'set appearance.branchLineWidth=5.0;',
                    'set scaleBar.isShown=false;',
            ]
            # orthospecies_tree_figtree_block.extend([
            #         'set nodeShapeExternal.colourAttribute="constrained";',
            #         'set nodeShapeExternal.isShown=true;',
            #         'set nodeShapeExternal.minSize=10.0;',
            #         'set nodeShapeExternal.scaleType=Width;',
            #         'set nodeShapeExternal.shapeType=Circle;',
            #         'set nodeShapeExternal.size=10.0;',
            #         'set nodeShapeExternal.sizeAttribute="Fixed";',
            #     ])
            orthospecies_tree_figtree_block.extend([
                # 'set tipLabels.colorAttribute="Name";',
                'set tipLabels.displayAttribute="Name";',
                'set tipLabels.fontName="sansserif";',
                'set tipLabels.fontSize=30;',
                'set tipLabels.fontStyle=0;',
                'set tipLabels.isShown=true;',
                'set tipLabels.significantDigits=4;',
                ])
            # orthospecies_tree_figtree_block.extend([
                # 'set legend.attribute="constrained";',
                # 'set legend.fontSize=10.0;',
                # 'set legend.isShown=true;',
                # 'set legend.significantDigits=4;',
            #     ])
            orthospecies_tree_figtree_block = "begin figtree;\n{}\nend;\n".format("\n".join(orthospecies_tree_figtree_block))
            orthospecies_tree.write(path="{}.species.nex".format(demo_output_prefix),
                    supplemental_blocks=[orthospecies_tree_figtree_block],
                    schema="nexus")

        # for post analysis assessment (not used by the inference program)
        config["test_info"] = collections.OrderedDict()
        config["test_info"]["species_leafsets"] = true_species_leafsets
        config["test_info"]["constrained_lineages"] = sorted(true_constrained_lineage_leaf_labels)
        config["test_info"]["unconstrained_lineages"] = sorted(true_unconstrained_lineage_leaf_labels)
        config["test_info"]["species_leafset_constraint_label_map"] = species_leafset_constraint_label_map
        config["test_info"]["species_partition_estimation_num_constrained_species"] = len(species_leafset_constraint_label_map)
        config["test_info"]["species_partition_estimation_num_constrained_lineages"] = len(true_constrained_lineage_leaf_labels)
        config["test_info"]["species_partition_estimation_num_unconstrained_lineages"] = len(true_unconstrained_lineage_leaf_labels)
        config["test_info"]["species_partition_estimation_num_unconstrained_lineages"] = len(true_unconstrained_lineage_leaf_labels)
        config["test_info"]["true_speciation_completion_rate"] = true_speciation_completion_rate # not actually used?
        config["test_info"]["true_species_leafsets"] = true_species_leafsets
        config["test_info"]["true_num_species"] = len(true_species_leafsets)
        with open(entry["run_config_filepath"], "w") as dest:
            json.dump(config, dest, indent=2)
        job_prefix = entry["tree_filepath"].replace(".nex", "")
        print(job_prefix)
        job_commands = []
        to_clean = []
        common_settings = {
                "run_config_filepath": entry["run_config_filepath"],
                "tree_filepath": entry["tree_filepath"],
                "num_lineages": len(entry["lineage_taxon_namespace"]),
                "num_species": len(entry["species_taxon_namespace"]),
                "true_speciation_completion_rate": true_speciation_completion_rate,
                "batch_id": batch_id,
                }
        to_clean.append(common_settings["run_config_filepath"])
        to_clean.append(common_settings["tree_filepath"])
        if args.underflow_protection:
            underflow_protection = "--underflow-protection"
        else:
            underflow_protection = ""
        job_kwargs = dict(common_settings)
        job_kwargs["underflow_protection"] = underflow_protection
        job_kwargs["delineate_results_filepath"] = job_prefix + ".partition-probs.json"
        to_clean.append(job_kwargs["delineate_results_filepath"])
        if args.specify_true_speciation_completion_rate:
            job_kwargs["speciation_completion_rate"] = "--speciation-completion-rate {}".format(true_speciation_completion_rate)
        else:
            job_kwargs["speciation_completion_rate"] = ""
        job_kwargs["post_analysis_performance_assessment_command"] = "python3 {}/spwd_evaluate_delineate_jobs.py".format(SCRIPT_DIR)
        job_commands.append(species_partition_estimation_job_template.format(**job_kwargs))
        job_kwargs["joint_performance_assessment_results_filepath"] = job_prefix + ".joint-partition-est-perf.tsv"
        job_commands.append(species_partition_estimation_joint_probability_analysis_template.format(**job_kwargs))
        job_filepath = job_prefix + ".job"
        if args.clean or args.very_clean:
            clean_command = ["rm", "-f"]
            clean_command.extend(to_clean)
            if args.very_clean:
                clean_command.append(job_filepath)
            job_commands.append(" ".join(clean_command))
        with open(job_filepath, "w") as dest:
            dest.write(template.format(jobs="\n".join(job_commands)))
    if not args.clean and not args.very_clean:
        with open(output_prefix + ".json", "w") as dest:
            json.dump(data, dest, indent=2)

if __name__ == '__main__':
    main()
