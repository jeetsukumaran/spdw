#! /usr/bin/env python
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
    for line in lines[-1:0:-1]:
        if line:
            break
    if not line:
        raise dendropy.DataError
    return line

def calculate_bpp_full_species_tree(
        src_tree_string,
        guide_tree,
        minimum_species_probability_threshold=0.95):
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
    assert treecompare.symmetric_difference(tree0, guide_tree) == 0
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
        elif nd.pp < minimum_species_probability_threshold:
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
            nd.label = "_".join(child_labels)
            nd.edge.length += len_to_add
            nd.child_labels = child_labels
            nd.clear_child_nodes()
            guide_tree_edge = guide_tree.bipartition_edge_map[nd.edge.bipartition]
            guide_tree_edge.collapse(adjust_collapsed_head_children_edge_lengths=True)
        else:
            for sub_nd in nd.child_node_iter():
                nodes_to_process.append(sub_nd)
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
    parser.add_argument("-t", "--minimum_species_probability_threshold",
            action="store",
            type=float,
            default=0.96,
            help="Mininum probability of splits to include (default=%(default)s)")

    args = parser.parse_args()
    with open(args.bpp_results_path) as src:
        bpp_out = src.read()
    guide_tree = dendropy.Tree.get(path=args.guide_tree_path, schema="nexus")

    posterior_guide_tree_string = extract_bpp_output_posterior_guide_tree_string(bpp_out)
    guide_tree, relabeled_tree = calculate_bpp_full_species_tree(
            src_tree_string=posterior_guide_tree_string,
            guide_tree=guide_tree,
            minimum_species_probability_threshold=args.minimum_species_probability_threshold)
    guide_tree.write(path="{}.collapsed.nex".format(args.output_prefix),
            schema="nexus")
    relabeled_tree.write(path="{}.relabeled.nex".format(args.output_prefix),
            schema="nexus")

if __name__ == '__main__':
    main()

