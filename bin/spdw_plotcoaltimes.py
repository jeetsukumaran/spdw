#! /usr/bin/env python
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
Species Delimitation Workshop: Plot coalescent times.
"""

import sys
import os
import random
import argparse
import dendropy

import pandas as pd
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spdw

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
    parser.add_argument("tree_files",
            action="store",
            type=str,
            nargs="+",
            metavar="TREEFILE",
            help="Path to tree files (default: read from standard input).")
    parser.add_argument("-f", "--input-format",
            type=str,
            default="newick",
            choices=["nexus", "newick"],
            help="Output data format (default='%(default)s')")

    args = parser.parse_args()
    args.output_prefix = None
    args.show_plot_on_screen = True
    data = []
    for src_idx, src_path in enumerate(args.tree_files):
        if src_path == "-":
            src = sys.stdin
        else:
            src = open(src_path)
        try:
            src_id = src.name
        except AttributeError:
            src_id = "<stdin>"
        with src:
            for tree in dendropy.Tree.yield_from_files(
                    files=[src],
                    schema=args.input_format):
                ages = tree.calc_node_ages(is_return_internal_node_ages_only=True)
                coalescence_events = sorted([nd for nd in tree if not nd.is_leaf()], key=lambda nd:nd.age)
                num_genes = len(coalescence_events) + 1
                # assert num_genes == len(tree.taxon_namespace)
                previous_age = 0.0
                coalescent_event_idx = 0
                while coalescence_events:
                    num_genes -= 1
                    coalescent_event_idx += 1
                    nd = coalescence_events.pop()
                    age = nd.age
                    waiting_time = nd.age - previous_age
                    data.append({
                        # "src_id": "I{:03d}".format(src_idx+1),
                        "src_id": src_id,
                        "num_genes": num_genes,
                        "coalescent_event_idx": coalescent_event_idx,
                        "age": age,
                        "waiting_time": waiting_time,
                        })
    df = pd.DataFrame(data)
    kwargs = {}
    if len(args.tree_files) > 1:
        kwargs["hue"] = "src_id"
    ax = sns.scatterplot(
            x="coalescent_event_idx",
            y="waiting_time",
            data=df,
            **kwargs
            )
    spdw.render_output(args, "Age")

if __name__ == '__main__':
    main()


