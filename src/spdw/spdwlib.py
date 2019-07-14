#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import random
from dendropy.utility import textprocessing
from dendropy.model import protractedspeciation
import collections
import itertools

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

class TableRowFormatter(object):

    def __init__(self,
            field_names,
            field_lengths=None,
            field_value_templates=None,
            field_templates=None,
            field_separator="  ",
            ):
        self.field_names = field_names
        self.field_lengths = field_lengths
        self.field_templates = field_templates
        self.field_value_templates = field_value_templates
        self.field_separator = field_separator

    @property
    def field_lengths(self):
        if self._field_lengths is None:
            self._field_lengths = [len(f) for f in self.field_names]
        return self._field_lengths
    @field_lengths.setter
    def field_lengths(self, value):
        self._field_lengths = value
    @field_lengths.deleter
    def field_lengths(self):
        del self._field_lengths

    @property
    def field_templates(self):
        if self._field_templates is None:
            self._field_templates = []
            for field_name_idx, field_name in enumerate(self.field_names):
                field_template = "{{:{}}}".format(self.field_lengths[field_name_idx])
                self._field_templates.append(field_template)
        return self._field_templates
    @field_templates.setter
    def field_templates(self, value):
        self._field_templates = value
    @field_templates.deleter
    def field_templates(self):
        del self._field_templates

    @property
    def field_value_templates(self):
        if self._field_value_templates is None:
            self._field_value_templates = ["{}" for i in range(len(self.field_names))]
        return self._field_value_templates
    @field_value_templates.setter
    def field_value_templates(self, value):
        self._field_value_templates = value
    @field_value_templates.deleter
    def field_value_templates(self):
        del self._field_value_templates

    def format_header_row(self):
        fields = []
        for field_idx, field_name in enumerate(self.field_names):
            fields.append(self.field_templates[field_idx].format(field_name))
        return self.field_separator.join(fields)

    def format_fields(self, field_values):
        fields = []
        for field_idx, field_value in enumerate(field_values):
            formatted_field_value = self.field_value_templates[field_idx].format(field_value)
            fields.append(self.field_templates[field_idx].format(formatted_field_value))
        return fields

    def format_row(self, field_values):
        return self.field_separator.join(self.format_fields(field_values))

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
        self.max_extant_lineages = kwargs.pop("max_extant_lineages", None)
        self.num_extant_orthospecies = kwargs.pop("num_extant_orthospecies", None)
        self.min_extant_orthospecies = kwargs.pop("min_extant_orthospecies", None)
        self.max_extant_orthospecies = kwargs.pop("max_extant_orthospecies", None)
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
        # make sure that the tree we generate has enough species
        lineage_tree, orthospecies_tree = self.psm.generate_sample(
                max_time=self.max_time,
                num_extant_lineages=self.num_extant_lineages,
                min_extant_lineages=self.min_extant_lineages,
                max_extant_lineages=self.max_extant_lineages,
                num_extant_orthospecies=self.num_extant_orthospecies,
                min_extant_orthospecies=self.min_extant_orthospecies,
                max_extant_orthospecies=self.max_extant_orthospecies,
                )
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

def generate_constraints_from_psm_trees(
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

def find_terminal_population_clades(tree):
    """
    Takes a tree with nodes *annotated" to indicate whether or not they should
    be collapsed and identifies the set of 'terminal population clades'.

    'Terminal population clades' are substrees descending from nodes with no
    ancestors in a collapsed state AND either: (a) no child nodes [i.e., a
    leaf] or (b) no descendent nodes in a non-collapsed or open states.
    Nodes corresponding to terminal populations themselves are, by definition,
    need to in a collapsed state unless they are leaves.

    Returns a dictionary mapping terminal population nodes to leaf nodes descending
    from the terminal population nodes.
    """
    # ensure if parent is collapsed, all children are collapsed
    for gnd in tree.preorder_node_iter():
        if gnd.parent_node is not None and gnd.parent_node.annotations["is_collapsed"].value:
            gnd.annotations["is_collapsed"] = True
    # identify the lowest nodes (closest to the tips) that are open, and
    # add its children if the children are (a) leaves; or (b) themselves
    # are closed
    terminal_population_clades = {}
    lineage_population_clade_map = {}
    for nd in tree.postorder_node_iter():
        if nd.annotations["is_collapsed"].value:
            continue
        for child in nd.child_node_iter():
            if child.is_leaf():
                terminal_population_clades[child] = set([child])
                lineage_population_clade_map[child] = child
            elif child.annotations["is_collapsed"].value:
                terminal_population_clades[child] = set([desc for desc in child.leaf_iter()])
                for lnd in terminal_population_clades[child]:
                    lineage_population_clade_map[lnd] = child
    for nd in tree:
        if nd not in terminal_population_clades:
            nd.annotations["population_id"] = "0"
            continue
        if nd.is_leaf():
            nd.annotations["population_id"] = nd.taxon.label
        else:
            nd.annotations["population_id"] = "+".join(desc.taxon.label for desc in nd.leaf_iter())
        # print("{}: {}".format(nd.annotations["population"], len(terminal_population_clades[nd])))
    return terminal_population_clades, lineage_population_clade_map

def identify_terminal_population_clade_species(terminal_population_clades):
    """
    Following population lineages being collapsed in an upstream analysis due
    to the structuring not being supported by the data, we may find lineages
    from two or more (true) nominal species have been collapsed together due to
    estimation error. This means that the original labeling / species
    identities cannot be used, and we need to establish new species
    identities.

    The criteria we shall use is that the new species identities comprise of
    the union of all species identities that are found together in the same
    collapsed population across all population.

    If lineages assigned to species "S1" and "S2" are collapsed into one
    population, while in another population we find lineages associated with
    species "S1" and "S3", and in no other population does "S1", "S2", and "S3"
    found with any other species except each other, then we establish a new
    species identity, "S1+S2+S3". All lineages associated with "S1", "S2", or
    "S3" are now referred to "S1+S2+S3". If it turns out that there is another
    population where are lineage assigned to "S3" and "S5" co-occurs, then the
    new species identity is "S1+S2+S3+S5", and all lineages associated with
    "S1", "S2", "S3", or "S5" are assigned to this new species.
    """
    terminal_population_clade_found_species = {}
    species_complexes = {}
    for terminal_population_clade in terminal_population_clades:
        found_species_labels = set()
        for lineage in terminal_population_clades[terminal_population_clade]:
            if not lineage.is_leaf():
                continue
            species_label, lineage_label = ProtractedSpeciationTreeGenerator.decompose_species_lineage_label(lineage.taxon.label)
            found_species_labels.add(species_label)
        for species_label in found_species_labels:
            if species_label in species_complexes:
                for ex2 in list(species_complexes[species_label]):
                    species_complexes[ex2].update(found_species_labels)
            else:
                species_complexes[species_label] = set(found_species_labels)
        terminal_population_clade_found_species[terminal_population_clade] = found_species_labels
    terminal_population_clade_species_identities = {}
    for nd in terminal_population_clades:
        x1 = None
        for spp_label in terminal_population_clade_found_species[nd]:
            if x1 is None:
                x1 = species_complexes[spp_label]
            else:
                assert x1 == species_complexes[spp_label]
        terminal_population_clade_species_identities[nd] = "+".join(x1)
    species_names = sorted(list(terminal_population_clade_species_identities.values()))
    return terminal_population_clade_species_identities

def generate_constraints_from_collapsed_tree(
        terminal_population_clades,
        min_unconstrained_leaves=None,
        max_unconstrained_leaves=None,
        num_unconstrained_leaves=None,
        max_tries=100,
        rng=None,
        ):
    if rng is None:
        rng = random.Random()
    num_tries = 0
    while True:
        num_tries += 1
        unconstrained_population_clades = set()
        unconstrained_lineages = set()
        lineage_pool = list(terminal_population_clades.keys())
        if min_unconstrained_leaves is None and num_unconstrained_leaves is None:
            lineages = list(itertools.chain.from_iterable(terminal_population_clades.values()))
            min_unconstrained_leaves = int(len(lineages) / 2)
        if num_unconstrained_leaves is not None:
            condition_fn = lambda: len(unconstrained_lineages) == num_unconstrained_leaves
        elif min_unconstrained_leaves is not None and max_unconstrained_leaves is not None:
            condition_fn = lambda: len(unconstrained_lineages) >= min_unconstrained_leaves and len(unconstrained_lineages) <= max_unconstrained_leaves
        elif min_unconstrained_leaves is not None:
            condition_fn = lambda: len(unconstrained_lineages) >= min_unconstrained_leaves
        else:
            sys.exit("Need to specify '--min-unconstrained-leaves' or '--num-unconstrained-leaves'")
        while lineage_pool and not condition_fn():
            node = lineage_pool.pop(rng.randint(0, len(lineage_pool)-1))
            unconstrained_lineages.update(terminal_population_clades[node])
            unconstrained_population_clades.add(node)
        if condition_fn:
            break
        else:
            print("Failed to meet condition (try {} of {})".format(num_tries, max_tries))
            if num_tries == max_tries:
                raise RunTimeError("Failed to meet constraint conditions within maximum try limit")
    constrained_population_clades = set([nd for nd in terminal_population_clades if nd not in unconstrained_population_clades])
    constrained_lineages = set(itertools.chain.from_iterable(terminal_population_clades[nd] for nd in constrained_population_clades))
    return {
            "constrained_population_clades": constrained_population_clades,
            "constrained_lineages": constrained_lineages,
            "unconstrained_population_clades": unconstrained_population_clades,
            "unconstrained_lineages": unconstrained_lineages,
            }

def format_population_clade_constraint_report(
        constraints,
        terminal_population_clades,
        terminal_population_clade_species_identities,
        ):
    population_nodes = sorted([nd for nd in terminal_population_clades], key=lambda nd: nd.annotations["population_id"].value)
    table = []
    for nd in population_nodes:
        row = {
                "Population": nd.annotations["population_id"].value,
                "Species": terminal_population_clade_species_identities[nd],
                "Status": "constrained" if nd in constraints["constrained_population_clades"] else "unconstrained",
        }
        table.append(row)
    msg = ["{} terminal population clades, {} organized into {} species and {} of unknown identity:".format(
        len(population_nodes),
        len(constraints["constrained_population_clades"]),
        len(set(terminal_population_clade_species_identities.values())),
        len(constraints["unconstrained_population_clades"]),
        )]
    formatted_table_rows = textprocessing.format_dict_table_rows(rows=table)
    for row in formatted_table_rows:
        msg.append("  {}".format(row))
    return "\n".join(msg)

def format_lineage_constraint_report(
        constraints,
        lineage_population_clade_map,
        terminal_population_clades,
        terminal_population_clade_species_identities,
        ):
    table = []
    constrained_lineages = constraints["constrained_lineages"]
    unconstrained_lineages = constraints["unconstrained_lineages"]
    for lineage_nd in lineage_population_clade_map:
        population_nd = lineage_population_clade_map[lineage_nd]
        if lineage_nd in constraints["constrained_lineages"]:
            status = "constrained"
        else:
            status = "unconstrained"
        row = {
                "Deme": lineage_nd.taxon.label,
                "Population":  population_nd.annotations["population_id"].value,
                "Species": terminal_population_clade_species_identities[population_nd],
                "Status": status,
                }
        table.append(row)
    formatted_table_rows= textprocessing.format_dict_table_rows(rows=table)
    msg = ["{} sub-population lineages organized into {} populations and {} species, with {} sub-population lineages of unknown species affinities:".format(
        len(lineage_population_clade_map),
        len(terminal_population_clades),
        len(set(terminal_population_clade_species_identities.values())),
        len(constrained_lineages),
        len(unconstrained_lineages),
        )]
    for row in formatted_table_rows:
        msg.append("  {}".format(row))
    return "\n".join(msg)

def decorate_constrained_collapsed_lineage_tree(
        lineage_tree,
        terminal_clades):
    for nd in lineage_tree:
        if nd not in terminal_clades:
            nd.annotations["population"] = "0"
            continue
        if nd.is_leaf():
            nd.annotations["population"] = nd.taxon.label
        else:
            nd.annotations["population"] = "+".join(desc.taxon.label for desc in nd.leaf_iter())
        nd.annotations["species"] = terminal_clade_species_identities[nd]

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
