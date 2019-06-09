#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import csv
import argparse
import dendropy

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

def calculate_bpp_full_species_tree(src_tree_string, minimum_species_probability_threshold=0.95):
    # Logic:
    # - Any internal node label is assumed to be a bpp annotation in the
    #   form of "#<float>" indicating the posterior probability of the node.
    # - If not found, assigned posterior probability of 0 if internal node or 1 if leaf.
    # - In post-order, collapse all nodes with less than threshold pp
    # - Revisit tree, assign taxon labels
    tree0 = dendropy.Tree.get(
            data=src_tree_string,
            schema="newick",
            suppress_external_node_taxa=True,
            suppress_internal_node_taxa=True,
            )
    for nd in tree0:
        if nd.is_leaf():
            nd.pp = 1.0
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
            nd.label = "_".join(c.label for c in desc_tips)
            nd.clear_child_nodes()
        else:
            for sub_nd in nd.child_node_iter():
                nodes_to_process.append(sub_nd)
    for nd in tree1.leaf_node_iter():
        nd.taxon = tree1.taxon_namespace.require_taxon(label=nd.label)
        nd.label = None

    return tree1

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bpp_results_path",
            action="store",
            help="Path to input file (BPP '.out.txt' file).")
    parser.add_argument("tree_output_path",
            action="store",
            help="Path to save tree in.")
    parser.add_argument("-t", "--minimum_species_probability_threshold",
            action="store",
            type=float,
            default=0.96,
            help="Mininum probability of splits to include (default=%(default)s)")

    args = parser.parse_args()
    if args.tree_output_path == "-":
        out = sys.stdout
    else:
        out = open(args.tree_output_path, "w")
    with open(args.bpp_results_path) as src:
        bpp_out = src.read()

    posterior_guide_tree_string = extract_bpp_output_posterior_guide_tree_string(bpp_out)
    collapsed_tree = calculate_bpp_full_species_tree(
            src_tree_string=posterior_guide_tree_string,
            minimum_species_probability_threshold=args.minimum_species_probability_threshold)
    collapsed_tree.write(file=out,
            schema="nexus")

if __name__ == '__main__':
    main()


def xmain():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "run_manifest_files",
            metavar="RUN-MANIFEST [RUN-MANIFEST [RUN-MANIFEST [...]]]",
            type=str,
            nargs="+",
            help="Path to BPP output file(s).")
    parser.add_argument(
            "-f", "--minimum-species-probability-threshold",
            type=float,
            default=0.95,
            help="Minimum probability for species node to be included in processed tree (default: %(default)s).")
    args = parser.parse_args()
    # if not args.bpp_output_file:
    #     src_path = "-"
    # else:
    #     src_path = args.bpp_output_file
    # if src_path == "-":
    #     s00.log("(reading from standard input)")
    #     src = sys.stdin
    # else:
    #     s00.log("Processing: {}".format(src_path))
    #     src = open(os.path.expandvars(os.path.expanduser(src_path)), "r")


    input_fields = [
            "speciation_initiation_from_orthospecies_rate",
            "speciation_initiation_from_incipient_species_rate",
            "speciation_completion_rate",
            "orthospecies_extinction_rate",
            "incipient_species_extinction_rate",
            "max_time",
            "max_extant_orthospecies",
            "num_extant_lineages",
            "num_extant_orthospecies",
            "source_tree_type",
            "population_size",
            "num_individuals_per_population",
            "total_number_of_individuals",
            "num_loci_per_individual",
            "mutation_rate_per_site",
            "num_source_tree_tips",
            "theta",
            "theta_prior_a",
            "theta_prior_b",
            "root_age",
            "tau_prior_a",
            "tau_prior_b",
            "source_tree_path",
            "results_filepath",
            "mcmc_filepath",
            "num_input_lineages",
            ]
    output_fields = ["num_estimated_orthospecies"] + input_fields
    out = sys.stdout
    sep = ","
    # out.write("{}\n".format(sep.join(output_fields)))
    num_manifest_entries = 0
    num_successes = 0
    fails = []
    output_rows = []
    for run_manifest_idx, run_manifest_path in enumerate(args.run_manifest_files):
        run_manifest_path = os.path.expandvars(os.path.expanduser(run_manifest_path))
        run_manifest_f = open(run_manifest_path, "r")
        run_manifest_dir = os.path.abspath(os.path.dirname(run_manifest_path))
        csv_reader = csv.DictReader(run_manifest_f)
        for input_row in csv_reader:
            num_manifest_entries += 1
            original_src_path = input_row["results_filepath"]
            full_src_path = os.path.join(run_manifest_dir, original_src_path)
            s00.log("Processing '{}': '{}'".format(run_manifest_path, original_src_path))
            try:
                src = open(full_src_path, "r")
                source_text = src.read()
                posterior_guide_tree_string = extract_bpp_output_posterior_guide_tree_string(source_text)
                collapsed_tree = calculate_bpp_full_species_tree(src_tree_string=posterior_guide_tree_string, minimum_species_probability_threshold=args.minimum_species_probability_threshold)
                num_tips = len(list(collapsed_tree.leaf_node_iter()))
                # s00.log("(number of tips = {})".format(num_tips))
                output_row = dict(input_row)
                output_row["num_estimated_orthospecies"] = num_tips
                output_rows.append(output_row)
                # out.write("{},{}\n".format(src_path, num_tips))
                # out_title = os.path.expandvars(os.path.expanduser(src_path)) + ".processed"
                # collapsed_tree.write(
                #         path=out_title + ".newick",
                #         schema="newick")
                num_successes += 1
            except (IOError, dendropy.DataError) as e:
                s00.log("Failed to process '{}': {}".format(full_src_path, e))
                fails.append( full_src_path )

    s00.log("{} out of {} results successfully processed".format(num_successes, num_manifest_entries))
    writer = csv.DictWriter(
            out,
            fieldnames=output_fields,
            restval="NA",
            delimiter=",",
            lineterminator=os.linesep,
            )
    writer.writeheader()
    writer.writerows(output_rows)

if __name__ == "__main__":
    main()

