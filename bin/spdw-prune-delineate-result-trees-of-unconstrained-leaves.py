#! /usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
##
##  Copyright 2019 Jeet Sukumaran.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
Species Delimitation Workshop: Show what the induced tree looks like.
"""

import sys
import os
import random
import argparse

import dendropy
from dendropy.simulate import treesim

__prog__ = os.path.basename(__file__)
__version__ = "1.0.0"
__description__ = __doc__
__author__ = 'Jeet Sukumaran'
__copyright__ = 'Copyright (C) 2019 Jeet Sukumaran.'

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)
    parser.add_argument("delineate_results_treefile_path")
    parser.add_argument(
            "--max-trees",
            action="store",
            default=None,
            help="Restrict processing to this number of trees.")
    parser.add_argument(
            "-o", "--output-prefix",
            action="store",
            default="constraint-info",
            help="Prefix for output files [default=%(default)s].")
    args = parser.parse_args()
    src = args.delineate_results_treefile_path
    dataset = dendropy.DataSet.get(path=src, schema="nexus", store_ignored_blocks=True)
    trees = dataset.tree_lists[0]
    filter_fn = lambda x: x.annotations["status"].value == "constrained" or x.annotations["constrained"].value == True
    supplemental_blocks = dataset.annotations["ignored_nexus_blocks"].value
    for nd in trees[0].leaf_node_iter():
        if filter_fn(nd):
            nd.annotations["known_species_label"] = "?"
        else:
            nd.annotations["known_species_label"] = nd.annotations["species"]
    if args.max_trees:
        trees = trees[:args.max_trees]
    trees.write(
            path=args.output_prefix + ".suppressed-labels.nex",
            schema="nexus",
            supplemental_blocks=supplemental_blocks)
    # trees = dendropy.TreeList.get(path=src, schema="nexus", store_ignored_blocks=True)
    new_tree_list = dendropy.TreeList(taxon_namespace=trees.taxon_namespace)
    for tidx, tree in enumerate(trees):
        # for leaf in tree.leaf_node_iter():
        #     print(leaf.annotations["status"].value)
        #     print(leaf.annotations["constrained"].value)
        filtered = tree.filter_leaf_nodes(filter_fn=filter_fn)
        # sys.stderr.write("Filtered: {}\n".format(filtered))
        new_tree_list.append(tree)
        if args.max_trees and idx >= args.max_trees:
            break
    new_tree_list.write(
            path=args.output_prefix + ".induced-trees.nex",
            schema="nexus",
            translate_tree_taxa=True,
            supplemental_blocks=supplemental_blocks,
            )

if __name__ == '__main__':
    main()


