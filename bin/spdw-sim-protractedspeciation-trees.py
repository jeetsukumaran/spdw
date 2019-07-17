#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import random
import argparse
from dendropy.model import protractedspeciation
from spdw import spdwlib

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output-prefix",
            action="store",
            default="pbd",
            help="Prefix for output files [default=%(default)s].")
    parser.add_argument("-f", "--output-format",
            action="store",
            default="nexus",
            choices=["nexus", "newick"],
            help="Format for output files [default=%(default)s].")
    parser.add_argument("--splitting-rate",
            action="store",
            type=float,
            default=0.10,
            help="Rate of origin of new lineages (population isolation rate) [default=%(default)s].")
    parser.add_argument("--extinction-rate",
            action="store",
            type=float,
            default=0.00,
            help="Rate of death of lineages [default=%(default)s].")
    parser.add_argument("--speciation-completion-rate",
            action="store",
            type=float,
            default=0.01,
            help="Rate at which lineage develop into independent species [default=%(default)s].")
    parser.add_argument("-n", "--num-replicates",
            type=int,
            default=30,
            help="Number of replicates per combination of parameters (default=%(default)s).")
    parser.add_argument("-z", "--random-seed",
            default=None,
            help="Random seed.")

    regime_group = parser.add_argument_group("Regime")
    regime_group.add_argument("--max-time",
            default=None,
            type=float,
            help="Source trees generated with this crown age.")
    regime_group.add_argument("--num-extant-lineages",
            default=None,
            type=int,
            help="Trees generated with exactly this number of tip lineages (incipient species + orthospecies).")
    regime_group.add_argument("--min-extant-lineages",
            default=None,
            type=int,
            help="Trees generated with at least this number of tip lineages (incipient species + orthospecies).")
    regime_group.add_argument("--max-extant-lineages",
            default=None,
            type=int,
            help="Trees generated with no more than this number of tip lineages (incipient species + orthospecies).")
    regime_group.add_argument("--num-extant-orthospecies",
            default=None,
            type=int,
            help="Trees generated with this number of orthospecies ('good' or true species).")
    regime_group.add_argument("--min-extant-orthospecies",
            default=2,
            type=int,
            help="Reject source trees with less than this number of orthospecies ('good' or true species).")
    regime_group.add_argument("--max-extant-orthospecies",
            default=None,
            type=int,
            help="Reject source trees with more than this number of orthospecies ('good' or true species).")
    args = parser.parse_args()
    selected_conditions = {}
    for kw in ("max_time", "num_extant_lineages", "num_extant_orthospecies"):
        selected_conditions[kw] = getattr(args, kw)
    if selected_conditions is None:
        sys.exit("Need to specify at least one of: '--max-time', '--num-extant-lineages', '--num-extant-orthospecies'")
    if args.random_seed is None:
        args.random_seed = random.randrange(sys.maxsize)
    rng = random.Random(args.random_seed)
    psm = spdwlib.ProtractedSpeciationTreeGenerator(
            splitting_rate=args.splitting_rate,
            extinction_rate=args.extinction_rate,
            speciation_completion_rate=args.speciation_completion_rate,
            max_time=args.max_time,
            num_extant_lineages=args.num_extant_lineages,
            min_extant_lineages=args.min_extant_lineages,
            max_extant_lineages=args.max_extant_lineages,
            num_extant_orthospecies=args.num_extant_orthospecies,
            min_extant_orthospecies=args.min_extant_orthospecies,
            max_extant_orthospecies=args.max_extant_orthospecies,
            rng=rng,
            )
    root_ages = []
    num_lineages = []
    num_species = []
    for tree_idx in range(args.num_replicates):
        lineage_tree, orthospecies_tree = psm.generate_sample()
        lineage_tree.calc_node_ages()
        print("Root age: {} ({} lineages, {} species))".format(lineage_tree.seed_node.age, len(lineage_tree.taxon_namespace), len(orthospecies_tree.taxon_namespace)))
        root_ages.append(lineage_tree.seed_node.age)
        num_lineages.append(len(lineage_tree.taxon_namespace))
        num_species.append(len(orthospecies_tree.taxon_namespace))
        output_prefix = "{}{:03d}".format(args.output_prefix, tree_idx+1)
        lineage_tree.write(
                path="{}.lineage.tre".format(output_prefix),
                schema=args.output_format)
        orthospecies_tree.write(
                path="{}.orthospecies.tre".format(output_prefix),
                schema=args.output_format)
    print("---")
    print("          Mean root age: {}".format(sum(root_ages)/len(root_ages)))
    print("Mean number of lineages: {}".format(sum(num_lineages)/len(num_lineages)))
    print(" Mean number of species: {}".format(sum(num_species)/len(num_species)))

if __name__ == '__main__':
    main()