#! /usr/bin/
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
Species Delimitation Workshop: Generate sequences.
"""

import sys
import os
import random
import argparse
import dendropy
from dendropy.interop import seqgen

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
    parser.add_argument("-t", "--tree-files",
            action="append",
            type=str,
            metavar="TREEFILE",
            help="Path to tree files (default: read from standard input).")
    parser.add_argument("-f", "--input-format",
            type=str,
            default="newick",
            choices=["nexus", "newick"],
            help="Input data format (default='%(default)s')")
    parser.add_argument("-n", "--num-characters-per-locus",
            type=int,
            default=1000,
            help="Number of characters sampled per locus (default: %(default)s).")
    parser.add_argument("--mutation-rate-per-site",
            type=float,
            default=0.00001,
            help="Per-site mutation rate (default: %(default)s).")
    parser.add_argument("-s", "--scale-branch-lengths",
            action="store",
            type=float,
            default=1.0,
            help="Scale branch lengths by this factor [default=%(default)s].")
    parser.add_argument("--num-replicates",
            type=int,
            default=1,
            help="Number of replicates (default: %(default)s).")
    parser.add_argument("-F", "--output-format",
            type=str,
            default="bpp",
            choices=["bpp", "nexus", "phylip", "fasta"],
            help="Input data format (default='%(default)s')")
    parser.add_argument("--concatenate",
            action="store_true",
            default=False,
            help="Concatenate the alignments across all genealogies")
    parser.add_argument("-z", "--random-seed",
            type=int,
            default=None,
            help="Seed for random number generator engine.")

    args = parser.parse_args()
    if not args.tree_files:
        sys.exit("Please specify path(s) to genealogy tree file(s)")
    sg = seqgen.SeqGen()
    sg.seq_len = args.num_characters_per_locus
    sg.scale_branch_lens = args.mutation_rate_per_site
    gene_trees = dendropy.TreeList()
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
            data = []
            gene_trees.read(
                    file=src,
                    schema=args.input_format,
                    rooting="force-rooted")
    if args.output_format == "bpp":
        for t in gene_trees.taxon_namespace:
            t.label = "^{}".format(t.label)
    for rep_idx in range(args.num_replicates):
        d0 = sg.generate(gene_trees)
        chars_filepath = "{}.{:03d}.chars".format(args.output_prefix, rep_idx+1)
        if args.output_format == "nexus":
            chars_filepath += ".nex"
            d0.write(path=chars_filepath, schema="nexus")
        elif args.output_format == "phylip":
            chars_filepath += ".phylip"
            d0.write(path=chars_filepath, schema="phylip")
        elif args.output_format == "fasta":
            chars_filepath += ".fasta"
            d0.write(path=chars_filepath, schema="fasta")
        elif args.output_format == "bpp":
            chars_filepath += ".txt"
            f = open(chars_filepath, "w")
            for cm in d0.char_matrices:
                cm.write(file=f, schema="phylip")
                f.write("\n")
        else:
            raise NotImplementedError

if __name__ == '__main__':
    main()



