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
Species Delimitation Workshop: Build FastSimCoal2 parameter file.
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

fsc_parfile_template = """\
//Number of population samples (demes)
{num_samples} samples to simulate :
//Population effective sizes (number of genes)
{population_effective_sizes}
//Samples sizes
{sample_sizes}
//Growth rates: negative growth implies population expansion
{growth_rates}
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, growth rate, migr. matrix
{num_historical_events} historical event
{historical_events}
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
DNA {num_loci} 0.0000 0.0005 0
"""

historical_events_template="""\
{div_time} {source} {sink} 1 1 0 0"""

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("tree_file")
    parser.add_argument("-o", "--output-prefix",
            action="store",
            default="fscrun",
            help="Prefix for otuput files (default=%(default)s).")
    parser.add_argument("-k", "--num-samples-per-pop",
            action="store",
            type=int,
            default=10,
            metavar="NUM-SAMPLES",
            help="Number of samples (default=%(default)s)")
    parser.add_argument("-N", "--pop-size", "--population-size",
            action="store",
            type=float,
            default=1000,
            metavar="POP-SIZE",
            help="Population size (default=%(default)s)")
    parser.add_argument("--num-loci",
            action="store",
            type=int,
            default=1,
            help="Population size (default=%(default)s)")
    parser.add_argument("-f", "--input-format",
            type=str,
            default="newick",
            choices=["nexus", "newick"],
            help="Input data format (default='%(default)s')")
    parser.add_argument("-s", "--scale-branch-lengths",
            action="store",
            type=float,
            default=1.0,
            help="Scale branch lengths by this factor [default=%(default)s].")
    args = parser.parse_args()
    tree = dendropy.Tree.get(path=args.tree_file, schema=args.input_format, rooting="force-rooted")
    for nd in tree:
        nd.edge.length = nd.edge.length * args.scale_branch_lengths
    tree.calc_node_ages()
    # convention: right child merges into left child
    historical_events = []
    leaf_idx = 0
    for nd_idx, nd in enumerate(tree.postorder_node_iter()):
        if nd.is_leaf():
            nd.pop_idx = leaf_idx
            leaf_idx += 1
        else:
            source_pop, sink_pop = nd.child_nodes() # assumes 2 children only!
            nd.pop_idx = sink_pop.pop_idx
            print("{}: {}=>{}".format(nd.pop_idx, source_pop.pop_idx, sink_pop.pop_idx))
            historical_event = historical_events_template.format(
                    div_time=nd.age,
                    source=source_pop.pop_idx,
                    sink=sink_pop.pop_idx)
            historical_events.append(historical_event)
    num_samples = len(tree.taxon_namespace)
    fsc_parfile_str = fsc_parfile_template.format(
            num_samples=num_samples,
            population_effective_sizes="\n".join(["{}".format(args.pop_size) for i in range(num_samples)]),
            sample_sizes="\n".join(["{}".format(args.num_samples_per_pop) for i in range(num_samples)]),
            growth_rates="\n".join(["0" for i in range(num_samples)]),
            num_historical_events=len(historical_events),
            historical_events="\n".join(historical_events),
            num_loci=args.num_loci)
    out = open("{}.par".format(args.output_prefix), "w")
    with out:
        out.write(fsc_parfile_str)

if __name__ == '__main__':
    main()



