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

    SPECIES_LABEL_FORMAT_TEMPLATE = "S{species_id}"
    LINEAGE_LABEL_FORMAT_TEMPLATE = "S{species_id}.L{lineage_id}"

    @staticmethod
    def decompose_species_lineage_label(label):
        parts = label.split(".")
        return parts[0], parts[1]

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
        self.rng = kwargs.pop("rng", random.Random())
        self.psm = protractedspeciation.ProtractedSpeciationProcess(
            speciation_initiation_from_orthospecies_rate=self.speciation_initiation_from_orthospecies_rate,
            speciation_initiation_from_incipient_species_rate=self.speciation_initiation_from_incipient_species_rate,
            speciation_completion_rate=self.speciation_completion_rate,
            orthospecies_extinction_rate=self.orthospecies_extinction_rate,
            incipient_species_extinction_rate=self.incipient_species_extinction_rate,
            species_label_format_template=ProtractedSpeciationTreeGenerator.SPECIES_LABEL_FORMAT_TEMPLATE,
            lineage_label_format_template=ProtractedSpeciationTreeGenerator.LINEAGE_LABEL_FORMAT_TEMPLATE,
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
                # if self.min_unconstrained_leaves:
                #     if len(lineage_tree.taxon_namespace) >= self.min_unconstrained_leaves:
                #         ok.append(True)
                #     else:
                #         ok.append(False)
                if self.min_extant_lineages:
                    if len(lineage_tree.taxon_namespace) >= self.min_extant_lineages:
                        ok.append(True)
                    else:
                        ok.append(False)
                if all(ok):
                    break
        return lineage_tree, orthospecies_tree

def build_tree_label_maps(orthospecies_tree=None, lineage_tree=None):
    species_lineage_label_map = collections.OrderedDict()
    lineage_species_label_map = {}
    if orthospecies_tree is not None:
        for k in sorted([t.label for t in orthospecies_tree.taxon_namespace]):
            species_lineage_label_map[k] = []
        for ond in orthospecies_tree.leaf_node_iter():
            species_lineage_label_map[ond.taxon.label] = sorted([lnd.taxon.label for lnd in ond.lineage_tree_nodes])
            for lnd in ond.lineage_tree_nodes:
                lineage_species_label_map[lnd.taxon.label] = ond.taxon.label
    elif lineage_tree is not None:
        for taxon in lineage_tree.taxon_namespace:
            label = taxon.label
            species_id, lineage_id = ProtractedSpeciationTreeGenerator.decompose_species_lineage_label(label)
            try:
                species_lineage_label_map[species_id].append(label)
            except KeyError:
                species_lineage_label_map[species_id] = [label]
            lineage_species_label_map[label] = species_id
    else:
        raise ValueError()
    return species_lineage_label_map, lineage_species_label_map

def generate_constraints(
        lineage_tree,
        orthospecies_tree,
        constraint_type,
        min_unconstrained_leaves,
        max_unconstrained_leaves,
        num_unconstrained_leaves,
        rng,
        species_lineage_label_map=None,
        lineage_species_label_map=None,
        ):
    if species_lineage_label_map is None or lineage_species_label_map is None:
        species_lineage_label_map, lineage_species_label_map = build_tree_label_maps(
                orthospecies_tree=orthospecies_tree)
    true_species_leafsets = sorted(species_lineage_label_map.values())
    species_leafset_constraints = None
    if constraint_type == "topological":
        lineage_tree_internal_nodes = [lnd for lnd in lineage_tree.postorder_internal_node_iter() if lnd is not lineage_tree.seed_node]
        rng.shuffle(lineage_tree_internal_nodes)
        unconstrained_lineage_leaf_labels = None
        for lineage_tree_internal_node in lineage_tree_internal_nodes:
            unconstrained_lineage_leaf_labels = set([lineage_tree_leaf_node.taxon.label for lineage_tree_leaf_node in lineage_tree_internal_node.leaf_iter()])
            if not num_unconstrained_leaves and not min_unconstrained_leaves and not max_unconstrained_leaves:
                break
            elif num_unconstrained_leaves:
                if len(unconstrained_lineage_leaf_labels) == num_unconstrained_leaves:
                    break
            elif min_unconstrained_leaves and max_unconstrained_leaves:
                if len(unconstrained_lineage_leaf_labels) >= min_unconstrained_leaves and len(unconstrained_lineage_leaf_labels) <= max_unconstrained_leaves:
                    break
            elif min_unconstrained_leaves and len(unconstrained_lineage_leaf_labels) >= min_unconstrained_leaves:
                break
            elif max_unconstrained_leaves and len(unconstrained_lineage_leaf_labels) <= max_unconstrained_leaves:
                break
            unconstrained_lineage_leaf_labels = None
        else:
            raise ValueError("Unable to meet min/max unconstrained leaves criteria.")
    elif constraint_type == "random":
        lineage_leaf_labels = [taxon.label for taxon in lineage_tree.taxon_namespace]
        if num_unconstrained_leaves:
            num_to_sample = num_unconstrained_leaves
        else:
            if min_unconstrained_leaves:
                min_count = min_unconstrained_leaves
            else:
                min_count = 1
            if max_unconstrained_leaves:
                max_count = max_unconstrained_leaves
            else:
                max_count = len(lineage_leaf_labels)
            num_to_sample = rng.randint(min_count, max_count)
        unconstrained_lineage_leaf_labels = rng.sample(lineage_leaf_labels, num_to_sample)
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

def decorate_lineage_tree(
        lineage_tree,
        lineage_species_label_map,
        unconstrained_lineage_leaf_labels,
        ):
    constrained_color, unconstrained_color = ColorMap.contrast_pairs[0]
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
    lineage_tree.figtree_block = "begin figtree;\n{}\nend;\n".format("\n".join(lineage_tree_figtree_block))

def decorate_orthospecies_tree(orthospecies_tree):
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
    orthospecies_tree.figtree_block = "begin figtree;\n{}\nend;\n".format("\n".join(orthospecies_tree_figtree_block))
