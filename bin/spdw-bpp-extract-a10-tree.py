#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import csv
import argparse
import dendropy
from dendropy.calculate import treecompare

"""
Species Delimitation Workshop: Extract BPP A10 Analysis Tree
"""

def extract_bpp_output_posterior_guide_tree_string(source_text):
    # pretty dumb logic for now: just locates the last non-blank line
    # works with: bp&p Version 3.1, April 2015, in "11" mode, i.e. infer tree and species delimitaiton
    lines = source_text.split("\n")
    for idx, line in enumerate(lines[-1::-1]):
        if line:
            break
    if not line:
        raise dendropy.DataParseError
    return line

def calculate_bpp_full_species_tree(
        src_tree_string,
        guide_tree,
        population_probability_threshold=0.95):
    # Logic:
    # - Any internal node label is assumed to be a bpp annotation in the
    #   form of "#<float>" indicating the posterior probability of the node.
    # - If not found, assigned posterior probability of 0 if internal node or 1 if leaf.
    # - In post-order, collapse all nodes with less than threshold pp
    # - Revisit tree, assign taxon labels
    tree0 = dendropy.Tree.get(
            data=src_tree_string,
            schema="newick",
            rooting="force-rooted",
            suppress_external_node_taxa=False,
            suppress_internal_node_taxa=True,
            taxon_namespace=guide_tree.taxon_namespace,
            )
    tree0.encode_bipartitions()
    guide_tree.encode_bipartitions()
    try:
        assert treecompare.symmetric_difference(tree0, guide_tree) == 0
    except AssertionError:
        print("[BPP tree]{}".format(tree0.as_string("newick")))
        print("[Guide tree]{}".format(guide_tree.as_string("newick")))
        raise
    for nd in tree0:
        edge_len = guide_tree.bipartition_edge_map[nd.edge.bipartition].length
        nd.edge.length = edge_len
        if nd.is_leaf():
            nd.pp = 1.0
            nd.label = nd.taxon.label
        elif nd.label:
            nd.pp = float(nd.label[1:])
            nd.label = None
        else:
            nd.pp = 0.0

    tree1 = dendropy.Tree(tree0)
    tree1.taxon_namespace = dendropy.TaxonNamespace()

    nodes_to_process = [tree1.seed_node]
    while nodes_to_process:
        nd = nodes_to_process.pop(0)
        if nd.is_leaf():
            pass
        elif nd.pp < population_probability_threshold:
            desc_tips = []
            for sub_nd in nd.leaf_iter():
                desc_tips.append(sub_nd)
            # nd.set_child_nodes(new_children)
            len_to_add = 0.0
            subnode = desc_tips[0]
            while subnode is not nd:
                len_to_add += subnode.edge.length
                subnode = subnode.parent_node
            child_labels = [c.label for c in desc_tips]
            nd.label = "+".join(child_labels)
            nd.edge.length += len_to_add
            nd.child_labels = child_labels
            nd.clear_child_nodes()
            guide_tree_edge = guide_tree.bipartition_edge_map[nd.edge.bipartition]
            guide_tree_edge.head_node.is_collapsed = True
        else:
            for sub_nd in nd.child_node_iter():
                nodes_to_process.append(sub_nd)
    collapsed_nodes = set()
    for gnd in guide_tree.postorder_node_iter():
        if getattr(gnd, "is_collapsed", False):
            len_to_add = 0.0
            desc_nd = gnd
            while True:
                try:
                    desc_nd = desc_nd._child_nodes[0]
                    len_to_add += desc_nd.edge.length
                except IndexError:
                    break
            gnd.edge.length += len_to_add
            for xnd in gnd.postorder_iter():
                if xnd is gnd:
                    continue
                xnd.edge.length = 1e-10
            gnd.annotations["is_collapsed"] = True
        else:
            gnd.annotations["is_collapsed"] = False
    for gnd in guide_tree.preorder_node_iter():
        if gnd.parent_node is not None and gnd.parent_node.annotations["is_collapsed"].value:
            gnd.annotations["is_collapsed"] = True
    for nd in tree1.leaf_node_iter():
        nd.taxon = tree1.taxon_namespace.require_taxon(label=nd.label)
        nd.label = None

    return guide_tree, tree1

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bpp_results_path",
            action="store",
            help="Path to input file (BPP '.out.txt' file).")
    parser.add_argument("guide_tree_path",
            action="store",
            help="Path to original guide_tree.")
    parser.add_argument("-o", "--output-prefix",
            action="store",
            default="bpp_a01",
            help="Prefix for output files (default=%(default)s).")
    parser.add_argument("-p", "--population-probability-threshold",
            action="store",
            type=float,
            default=0.95,
            help="Mininum probability of splits to include (default=%(default)s)")

    args = parser.parse_args()
    with open(args.bpp_results_path) as src:
        bpp_out = src.read()
    posterior_guide_tree_string = extract_bpp_output_posterior_guide_tree_string(bpp_out)
    guide_tree = dendropy.Tree.get(path=args.guide_tree_path, schema="nexus")
    guide_tree, relabeled_tree = calculate_bpp_full_species_tree(
            src_tree_string=posterior_guide_tree_string,
            guide_tree=guide_tree,
            population_probability_threshold=args.population_probability_threshold)
    guide_tree.write(path="{}.collapsed.nex".format(args.output_prefix),
            schema="nexus")
    relabeled_tree.write(path="{}.relabeled.nex".format(args.output_prefix),
            schema="nexus")

if __name__ == '__main__':
    main()

