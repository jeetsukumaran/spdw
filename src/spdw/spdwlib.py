#! /usr/bin/env python
# -*- coding: utf-8 -*-

import random
from dendropy.model import protractedspeciation
import collections

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

class ProtractedSpeciationTreeGenerator(object):

    def __init__(self, **kwargs):
        self.speciation_initiation_from_orthospecies_rate = kwargs.pop("splitting_rate", 0.01)
        self.speciation_initiation_from_incipient_species_rate = self.speciation_initiation_from_orthospecies_rate
        self.speciation_completion_rate = kwargs.pop("speciation_completion_rate", 0.01)
        self.orthospecies_extinction_rate = kwargs.pop("extinction_rate", 0.0)
        self.incipient_species_extinction_rate = self.orthospecies_extinction_rate
        self.max_time = kwargs.pop("max_time", None)
        self.num_extant_lineages = kwargs.pop("num_extant_lineages", None)
        self.min_extant_lineages = kwargs.pop("min_extant_lineages", None)
        self.num_extant_orthospecies = kwargs.pop("num_extant_orthospecies", None)
        self.min_extant_orthospecies = kwargs.pop("min_extant_orthospecies", 2)
        self.min_unconstrained_leaves = kwargs.pop("min_unconstrained_leaves", None)
        self.max_unconstrained_leaves = kwargs.pop("max_unconstrained_leaves", None)
        self.num_unconstrained_leaves = kwargs.pop("num_unconstrained_leaves", None)
        self.rng = kwargs.pop("rng", random.Random())
        self.psm = protractedspeciation.ProtractedSpeciationProcess(
            speciation_initiation_from_orthospecies_rate=self.speciation_initiation_from_orthospecies_rate,
            speciation_initiation_from_incipient_species_rate=self.speciation_initiation_from_incipient_species_rate,
            speciation_completion_rate=self.speciation_completion_rate,
            orthospecies_extinction_rate=self.orthospecies_extinction_rate,
            incipient_species_extinction_rate=self.incipient_species_extinction_rate,
            rng=self.rng,
            )

    def generate_sample(self,):
        selected_condition = None
        for kw in (self.max_time, self.num_extant_lineages, self.num_extant_orthospecies):
            if kw is not None:
                if selected_condition:
                    sys.exit("Need to specify only one of: 'max_time', 'num_extant_lineages', 'num_extant_orthospecies'")
                selected_condition = kw
        if selected_condition is None:
            sys.exit("Need to specify at least one of: 'max_time', 'num_extant_lineages', 'num_extant_orthospecies'")

        while True:
            # make sure that the tree we generate has enough species
            lineage_tree, orthospecies_tree = self.psm.generate_sample(
                    max_time=self.max_time,
                    num_extant_lineages=self.num_extant_lineages,
                    num_extant_orthospecies=self.num_extant_orthospecies,
                    )
            if len(orthospecies_tree.taxon_namespace) >= self.min_extant_orthospecies:
                ok = []
                if self.min_unconstrained_leaves:
                    if len(lineage_tree.taxon_namespace) >= self.min_unconstrained_leaves:
                        ok.append(True)
                    else:
                        ok.append(False)
                if self.min_extant_lineages:
                    if len(lineage_tree.taxon_namespace) >= self.min_extant_lineages:
                        ok.append(True)
                    else:
                        ok.append(False)
                if all(ok):
                    break
        return lineage_tree, orthospecies_tree

    def build_label_maps(self, orthospecies_tree):
        species_lineage_label_map = collections.OrderedDict()
        lineage_species_label_map = {}
        for k in sorted([t.label for t in orthospecies_tree.taxon_namespace]):
            species_lineage_label_map[k] = []
        for ond in orthospecies_tree.leaf_node_iter():
            species_lineage_label_map[ond.taxon.label] = sorted([lnd.taxon.label for lnd in ond.lineage_tree_nodes])
            for lnd in ond.lineage_tree_nodes:
                lineage_species_label_map[lnd.taxon.label] = ond.taxon.label
        return species_lineage_label_map, lineage_species_label_map

    def generate_constraints(self,
            lineage_tree,
            constraint_type,
            species_lineage_label_map,
            lineage_species_label_map,
            ):
        true_species_leafsets = sorted(species_lineage_label_map.values())
        species_leafset_constraints = None
        if constraint_type == "topological":
            lineage_tree_internal_nodes = [lnd for lnd in lineage_tree.postorder_internal_node_iter() if lnd is not lineage_tree.seed_node]
            rng.shuffle(lineage_tree_internal_nodes)
            unconstrained_lineage_leaf_labels = None
            for lineage_tree_internal_node in lineage_tree_internal_nodes:
                unconstrained_lineage_leaf_labels = set([lineage_tree_leaf_node.taxon.label for lineage_tree_leaf_node in lineage_tree_internal_node.leaf_iter()])
                if not self.num_unconstrained_leaves and not self.min_unconstrained_leaves and not self.max_unconstrained_leaves:
                    break
                elif self.num_unconstrained_leaves:
                    if len(unconstrained_lineage_leaf_labels) == self.num_unconstrained_leaves:
                        break
                elif self.min_unconstrained_leaves and self.max_unconstrained_leaves:
                    if len(unconstrained_lineage_leaf_labels) >= self.min_unconstrained_leaves and len(unconstrained_lineage_leaf_labels) <= self.max_unconstrained_leaves:
                        break
                elif self.min_unconstrained_leaves and len(unconstrained_lineage_leaf_labels) >= self.min_unconstrained_leaves:
                    break
                elif self.max_unconstrained_leaves and len(unconstrained_lineage_leaf_labels) <= self.max_unconstrained_leaves:
                    break
                unconstrained_lineage_leaf_labels = None
            else:
                raise ValueError("Unable to meet min/max unconstrained leaves criteria.")
        elif constraint_type == "random":
            lineage_leaf_labels = [taxon.label for taxon in lineage_tree.taxon_namespace]
            if self.num_unconstrained_leaves:
                num_to_sample = self.num_unconstrained_leaves
            else:
                if self.min_unconstrained_leaves:
                    min_count = self.min_unconstrained_leaves
                else:
                    min_count = 1
                if self.max_unconstrained_leaves:
                    max_count = self.max_unconstrained_leaves
                else:
                    max_count = len(lineage_leaf_labels)
                num_to_sample = self.rng.randint(min_count, max_count)
            unconstrained_lineage_leaf_labels = self.rng.sample(lineage_leaf_labels, num_to_sample)
        constrained_lineage_leaf_labels = sorted([lineage_tree_leaf_node.taxon.label for lineage_tree_leaf_node in lineage_tree.leaf_node_iter() if lineage_tree_leaf_node.taxon.label not in unconstrained_lineage_leaf_labels])
        species_leafset_constraint_label_map = collections.OrderedDict()
        for lineage_leaf_label in constrained_lineage_leaf_labels:
            true_species_label = lineage_species_label_map[lineage_leaf_label]
            try:
                species_leafset_constraint_label_map[true_species_label].append(lineage_leaf_label)
            except KeyError:
                species_leafset_constraint_label_map[true_species_label] = [lineage_leaf_label]
        for species_label in species_leafset_constraint_label_map:
            species_leafset_constraint_label_map[species_label].sort()
        species_leafset_constraints = []
        for sp in species_lineage_label_map:
            if sp in species_leafset_constraint_label_map:
                species_leafset_constraints.append(sorted(species_leafset_constraint_label_map[sp]))
        return species_leafset_constraints, constrained_lineage_leaf_labels, unconstrained_lineage_leaf_labels, species_leafset_constraint_label_map

        # extra files for demo
        if args.write_extra_for_demo:
            constrained_color, unconstrained_color = spdwlib.ColorMap.contrast_pairs[0]
            demo_output_prefix = "{}.{:04d}.demo".format(output_prefix, tree_idx+1)
            true_unconstrained_lineage_leaf_label_set = set(unconstrained_lineage_leaf_labels)
            for nd in lineage_tree:
                if nd.is_leaf():
                    is_constrained = nd.taxon.label in true_unconstrained_lineage_leaf_label_set
                    nd.annotations["constrained"] = is_constrained
                    species_label = lineage_species_label_map[nd.taxon.label]
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
        config["test_info"]["constrained_lineages"] = sorted(constrained_lineage_leaf_labels)
        config["test_info"]["unconstrained_lineages"] = sorted(unconstrained_lineage_leaf_labels)
        config["test_info"]["species_leafset_constraint_label_map"] = species_leafset_constraint_label_map
        config["test_info"]["species_partition_estimation_num_constrained_species"] = len(species_leafset_constraint_label_map)
        config["test_info"]["species_partition_estimation_num_constrained_lineages"] = len(constrained_lineage_leaf_labels)
        config["test_info"]["species_partition_estimation_num_unconstrained_lineages"] = len(unconstrained_lineage_leaf_labels)
        config["test_info"]["species_partition_estimation_num_unconstrained_lineages"] = len(unconstrained_lineage_leaf_labels)
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
        job_kwargs["post_analysis_performance_assessment_command"] = "spdw-evaluate-delineate-jobs.py".format(SCRIPT_DIR)
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
