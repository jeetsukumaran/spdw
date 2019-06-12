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
Species Delimitation Workshop: Fix BPP traces
"""

import sys
import os
import argparse

__prog__ = os.path.basename(__file__)
__version__ = "1.0.0"
__description__ = __doc__
__author__ = 'Jeet Sukumaran'
__copyright__ = 'Copyright (C) 2019 Jeet Sukumaran.'

def get_max_cols(src_path):
    max_cols = 0
    with open(src_path) as src:
        for row in src:
            row = row.strip()
            cols = row.split("\t")
            if len(cols) > max_cols:
                max_cols = len(cols)
    return max_cols

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("input_files",
            action="store",
            type=str,
            nargs="+",
            help="Path to BPP 'mcmc.txt' files")
    args = parser.parse_args()
    for src_path in args.input_files:
        output_path = src_path + ".traces"
        out = open(output_path, "w")
        max_cols = get_max_cols(src_path)
        with open(src_path) as src:
            first_row = src.readline()
            col_names = first_row.strip().split("\t")
            for new_col_idx in range(max_cols - len(col_names)):
                col_names.insert(-2, "col{}".format(new_col_idx+1))
            out.write("\t".join(col_names))
            out.write("\n")
            for row in src:
                row = row.strip()
                cols = row.split("\t")
                for new_col_idx in range(max_cols - len(cols)):
                    cols.insert(-2, "0.0000")
                new_row = "\t".join(cols)
                out.write("{}".format(new_row))
                out.write("\n")

if __name__ == '__main__':
    main()


