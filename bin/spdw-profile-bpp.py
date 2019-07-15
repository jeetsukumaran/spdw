#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import collections
import csv
import re
import platform
import random
import collections
import argparse
from dendropy.calculate import popgenstat
from dendropy.utility import textprocessing
from spdw import spdwlib
import dendropy

def _log(msg):
    sys.stderr.write("- {}\n".format(msg))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sequences_file",
            metavar="FILEPATH",
            help="Path to sequence data file.")
    # parser.add_argument("control_file",
    #         metavar="FILEPATH",
    #         help="Path to control file.")
    args = parser.parse_args()
    # control_filepath = os.path.expanduser(os.path.expandvars(args.control_file))
    # root_dir = os.path.parentdir(control_filepath)
    seqs_filepath = os.path.expanduser(os.path.expandvars(args.sequences_file))
    data = dendropy.DataSet.get(path=seqs_filepath,
            schema="multiphylip",
            data_type="dna")
    field_names = (
            "Locus",
            "pi",
            "theta",
            "TajD",
            )
    table_row_formatter = spdwlib.TableRowFormatter(
            field_names=field_names,
            field_lengths=[5] + [11] * (len(field_names)-1),
            field_value_templates = ["{}"] + ["{:<11.8f}"] * (len(field_names)-1),
            field_separator = "    ",
            )
    print("{}".format(table_row_formatter.format_header_row()))
    locus_thetas = []
    locus_profile_info = []
    for cm_idx, cm in enumerate(data.char_matrices):
        try:
            td = popgenstat.tajimas_d(cm)
        except ZeroDivisionError as e:
            td = "N/A"
        nucleotide_diversity = popgenstat.nucleotide_diversity(cm)
        # locus_thetas.append(nucleotide_diversity)
        wtheta = popgenstat.wattersons_theta(cm) / cm.sequence_size
        locus_thetas.append(wtheta)
        field_values = [
            cm_idx+1,
            nucleotide_diversity,
            locus_thetas[-1],
            td,
        ]
        assert len(field_values) == len(field_names)
        profile = table_row_formatter.format_row(field_values=field_values)
        print("{}".format(profile))
        locus_profile_info.append(dict(zip(field_names, field_values)))

if __name__ == "__main__":
    main()


