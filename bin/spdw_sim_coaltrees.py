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
Species Delimitation Workshop: Simulate trees under Kingman's neutral coalescent.
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
    parser.add_argument("output_prefix")
    parser.add_argument("-k", "--num-tips",
            action="store",
            type=int,
            default=10,
            metavar="NUM-TIPS",
            help="Number of samples (default=%(default)s)")
    parser.add_argument("-N", "--pop-size", "--population-size",
            action="store",
            type=float,
            default=1.0,
            metavar="POP-SIZE",
            help="Population size (default=%(default)s)")
    parser.add_argument("--num-reps",
            action="store",
            type=int,
            default=10,
            metavar="NUM-REPS",
            help="Number of replicates (default=%(default)s)")
    parser.add_argument("-z", "--random-seed",
            type=int,
            default=None,
            help="Random seed.")

    args = parser.parse_args()

    if args.random_seed is None:
        args.random_seed = random.randint(0, sys.maxsize)
    rng = random.Random(args.random_seed)
    tns = dendropy.TaxonNamespace(["G{:03d}".format(i+1) for i in range(args.num_tips)])
    trees = dendropy.TreeList(taxon_namespace=tns)
    if args.output_prefix == "-":
        coal_out = sys.stdout
    else:
        coal_out = open("{}.coal.trees".format(args.output_prefix), "w")
    for rep_id in range(args.num_reps):
        tree = treesim.pure_kingman_tree(
                taxon_namespace=tns,
                pop_size=args.pop_size,
                rng=rng)
        tree.write(file=coal_out, schema="newick")

if __name__ == '__main__':
    main()


