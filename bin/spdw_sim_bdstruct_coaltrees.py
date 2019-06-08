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
Species Delimitation Workshop: Simulate trees under the censored coalescent structured by a birth-death process.
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
    parser.add_argument("-k", "--num-genes-per-pop",
            action="store",
            type=int,
            default=10,
            metavar="NUM-GENES-PER-POP",
            help="Number of samples per population (default=%(default)s)")
    parser.add_argument("-p", "--num-pops",
            action="store",
            type=int,
            default=10,
            metavar="NUM-POPS",
            help="Number of populations (default=%(default)s)")
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
    tns = dendropy.TaxonNamespace()
    if args.output_prefix == "-":
        coal_out = sys.stdout
        bd_out = sys.stderr
    else:
        coal_out = open("{}.coal.trees".format(args.output_prefix), "w")
        bd_out = open("{}.bd.trees".format(args.output_prefix), "w")
    for rep_id in range(args.num_reps):
        pop_tree = treesim.birth_death_tree(
                birth_rate=0.10,
                death_rate=0.00,
                num_extant_tips=args.num_pops,
                gsa_ntax=100,
                rng=rng,
                )
        # sys.stderr.write("{}\n".format(pop_tree.seed_node.age))
        for nd in pop_tree.postorder_node_iter():
            if nd.is_leaf():
                nd.num_genes = args.num_genes_per_pop
                nd.edge.pop_size = args.pop_size / args.num_pops
            else:
                nd.edge.pop_size = sum([ch.edge.pop_size for ch in nd.child_nodes()])
        pop_tree.calc_node_ages()
        pop_tree.write(file=bd_out, schema="newick")
        ctree, ptree = treesim.constrained_kingman_tree(
                pop_tree=pop_tree,
                rng=rng)
        ctree.write(file=coal_out, schema="newick")

if __name__ == '__main__':
    main()


