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
from dendropy.model import reconcile
from dendropy.model import birthdeath
from dendropy.interop import seqgen
from dendropy.calculate import popgenstat
import dendropy

def _log(msg):
    sys.stderr.write("- {}\n".format(msg))

def _open_output_file_for_csv_writer(filepath, append=False):
    if filepath is None or filepath == "-":
        out = sys.stdout
    elif sys.version_info >= (3,0,0):
        out = open(filepath,
                "a" if append else "w",
                newline='')
    else:
        out = open(filepath, "ab" if append else "wb")
    return out

"""
Both parameters theta s and
tau s in the MSC model are measured by genetic distance, the expected number of
mutations/substitutions per site. For example, for the human species, theta H ~= 0.0006, which
means that two random sequences from the human population are different at ~0.06% of
sites, less than 1 difference per kb. A sensible diffuse prior is then ''thetaprior = 3 0.002''
with mean 0.001.
"""

BPP_TEMPLATE = """\

          seed =  -1

       seqfile = {chars_filepath}
      Imapfile = {imap_filepath}
       outfile = {out_filepath}
      mcmcfile = {mcmc_filepath}

* speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
  speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
          speciestree = 0        * species tree NNI/SPR
*        speciestree = 1  0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio


     speciesmodelprior = 1         * 0: uniform labeled histories; 1:uniform rooted trees

  species&tree = {num_populations}  {population_labels}
                     {num_individuals_per_population}
                 {bpp_guide_tree}

       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = {num_loci}    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = {theta_prior_a} {theta_prior_b}   # invgamma(a, b) for theta, mean = {theta_prior_mean}
      tauprior = {tau_prior_a}   {tau_prior_b}       # invgamma(a, b) for root tau & Dirichlet(a) for other tau's, mean = {tau_prior_mean}; root age (raw, unscaled) = {root_age}

      finetune =  1: 3 0.003 0.002 0.00002 0.005 0.9 0.001 0.001 # finetune for GBtj, GBspr, theta, tau, mix

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 4000
      sampfreq = 100
       nsample = 1000
"""

def try_to_coerce_to_float(v):
    try:
        return float(v)
    except ValueError:
        if v == "None":
            return "NA"
        else:
            return v

def generate_contained_trees(
        containing_tree,
        contained_taxon_namespace=None,
        population_size=1,
        num_subpopulation_lineages_per_population_fn=None,
        subpopulation_root_age=None,
        num_individuals_per_terminal_lineage=4,
        num_gene_trees=5,
        rng=None):
    assert len(containing_tree.taxon_namespace) > 0
    if num_subpopulation_lineages_per_population_fn is not None:
        pseudopopulation_tree = dendropy.Tree(containing_tree)
        pseudopopulation_tree.taxon_namespace = dendropy.TaxonNamespace()
        for parent_node in pseudopopulation_tree.postorder_node_iter():
            if not parent_node.is_leaf():
                parent_node.annotations["true_population_id"] = "null"
                continue
            parent_taxon = parent_node.taxon
            parent_node.taxon = None
            # parent_node.label = parent_taxon.label
            parent_node.annotations["true_population_id"] = parent_taxon.label
            num_subpops = num_subpopulation_lineages_per_population_fn()
            if num_subpops == 1:
                pseudopopulation_label = "{}.sub1".format(parent_taxon.label)
                parent_node.taxon = pseudopopulation_tree.taxon_namespace.require_taxon(label=pseudopopulation_label)
                parent_node.taxon.true_population_label = parent_taxon.label
                parent_node.edge.length += subpopulation_root_age
            else:
                subtree = birthdeath.birth_death_tree(
                        num_extant_tips=num_subpops,
                        birth_rate=0.1,
                        death_rate=0.0)
                # print(subtree.as_string("newick"))
                subtree.calc_node_ages()
                original_root_age = subtree.seed_node.age + subtree.seed_node.edge.length
                for subtree_node in subtree.postorder_node_iter():
                    subtree_node.age = (subtree_node.age / original_root_age) * subpopulation_root_age
                subtree.set_edge_lengths_from_node_ages(error_on_negative_edge_lengths=True)
                subtree_leaf_idx = 0
                for subtree_node in subtree.postorder_node_iter():
                    subtree_node.annotations["true_population_id"] = parent_taxon.label
                    if subtree_node.parent_node is subtree.seed_node:
                        parent_node.add_child(subtree_node)
                    if subtree_node.is_leaf():
                        subtree_leaf_idx += 1
                        pseudopopulation_label = "{}.sub{}".format(parent_taxon.label, subtree_leaf_idx)
                        subtree_node.taxon = pseudopopulation_tree.taxon_namespace.require_taxon(label=pseudopopulation_label)
                        subtree_node.taxon.true_population_label = parent_taxon.label
        containing_tree = pseudopopulation_tree
    else:
        for nd in containing_tree.postorder_node_iter():
            if nd.taxon is not None:
                nd.taxon.true_population_label = nd.taxon.label
                nd.annotations["true_population_id"] = nd.taxon.label
    if contained_taxon_namespace is None:
        contained_taxon_namespace = dendropy.TaxonNamespace()
    contained_to_containing_map = {}
    pop_samples = collections.Counter()
    for sp_idx, sp_tax in enumerate(containing_tree.taxon_namespace):
        for gidx in range(num_individuals_per_terminal_lineage):
            pop_samples[sp_tax.true_population_label] += 1
            sample_idx = pop_samples[sp_tax.true_population_label]
            glabel = "{sp}_i{ind}^{sp}_i{ind}".format(
                    sp=sp_tax.true_population_label,
                    ind=sample_idx)
            # glabel = "{sp}^{sp}_{ind}".format(sp=sp_tax.label, ind=gidx+1)
            g = contained_taxon_namespace.require_taxon(label=glabel)
            g.population_label = sp_tax.label
            contained_to_containing_map[g] = sp_tax
    ct = reconcile.ContainingTree(
            containing_tree=containing_tree,
            contained_taxon_namespace=contained_taxon_namespace,
            contained_to_containing_taxon_map=contained_to_containing_map)
    gene_trees = dendropy.TreeList(taxon_namespace=contained_taxon_namespace)
    for gtidx in range(num_gene_trees):
        gt = ct.simulate_contained_kingman(
                default_pop_size=population_size,
                rng=rng)
        gene_trees.append(gt)
    return containing_tree, gene_trees


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_trees",
            metavar="SOURCE_TREEFILE [SOURCE_TREEFILE [SOURCE_TREEFILE]]",
            nargs="+",
            help="Path to containing tree files. Specify '-' to read from standard input.")
    parser.add_argument("-f", "--input-format",
            default="nexus",
            dest="schema",
            help="Input trees format (default: $(default)s).")
    # parser.add_argument("--split-populations",
    #         help="Introduce splits between individuals within populations in *guide tree* passed to BPP."
    #              " The true tree, on which the sequences etc., remains the same. The purpose of this finer"
    #              " grained oversplitting is to allow BPP to collapse splits not supported by actual"
    #              " population structure.")
    parser.add_argument("-z", "--random-seed",
            type=int,
            default=None,
            help="Seed for random number generator engine.")
    parser.add_argument("-o", "--output-prefix",
            default="bpprun",
            help="Run title (default: '%(default)s')")
    data_options = parser.add_argument_group("Data Options")
    data_options.add_argument("--population-size",
            type=float,
            default=1e6,
            help="Population size [default: %(default)s].")
    data_options.add_argument(
            "--min-subpopulation-lineages-per-population",
            type=int,
            default=0,
            help="In guide tree passed to BPP, divide individuals from single populations into at least this number of multiple nominal populations (with minimal to no structure between them), allowing BPP to collapse these (default: %(default)s; i.e. do not create multiple subpopulation_lineages).")
    data_options.add_argument(
            "--max-subpopulation-lineages-per-population",
            type=int,
            default=0,
            help="In guide tree passed to BPP, divide individuals from single populations into at most this number of multiple nominal populations (with minimal to no structure between them), allowing BPP to collapse these (default: %(default)s; i.e. do not create multiple subpopulation_lineages).")
    data_options.add_argument(
            "--subpopulation-root-age",
            type=float,
            default=1e-10,
            help="Determines the amount of structuring in subpopulations; set to 0 for no structure at all [default: %(default)s].")
    data_options.add_argument("--num-individuals-per-terminal-lineage",
            type=int,
            default=4,
            help="Number of individuals sampled per terminal lineage within each population [default: %(default)s].")
    data_options.add_argument("--num-loci-per-individual",
            type=int,
            default=10,
            help="Number of loci sampled per individual [default: %(default)s].")
    data_options.add_argument("--num-characters-per-locus",
            type=int,
            default=1000,
            help="Number of characters sampled per locus [default: %(default)s].")
    data_options.add_argument("--mutation-rate-per-site",
            type=float,
            # default=0.00001,
            default=1e-8,
            help="Per-site mutation rate [default: %(default)s].")
    data_options.add_argument("--no-scale-tree-by-mutation-rate",
            action="store_true",
            help="Do not scale tree by mutation rate.")
    summary_options = parser.add_argument_group("Summary Options")
    summary_options.add_argument("-p", "--population-probability-threshold",
            action="store",
            type=float,
            default=0.95,
            help="Mininum probability of splits to include when summarizing tree (default=%(default)s)")
    job_options = parser.add_argument_group("Job Options")
    job_options.add_argument("--job-preamble",
            metavar="JOB-PREAMBLE-FILE",
            help="Path to file with job preamble (commands to insert in job file before main operations).")
    args = parser.parse_args()

    if args.random_seed is None:
        random_seed = random.randint(0, sys.maxsize-1)
    else:
        random_seed = args.random_seed

    rng = random.Random(random_seed)
    _log("Random seed: {}".format(random_seed))

    sg = seqgen.SeqGen()
    sg.seq_len = args.num_characters_per_locus
    sg.scale_branch_lens = args.mutation_rate_per_site

    if "-" in args.source_trees:
        filepaths = sys.stdin.read().split("\n")
        args.source_trees.remove("-")
    else:
        filepaths = []

    manifest_entries = []
    filepaths.extend(args.source_trees)
    for idx, filepath in enumerate(filepaths):
        job_title = "{}_{:05d}".format(args.output_prefix, idx+1)
        _log("{} of {}: {}: {}".format(idx+1, len(filepaths), job_title, filepath))
        source_tree = dendropy.Tree.get(
                path=filepath,
                schema=args.schema,
                extract_comment_metadata=True,
                preserve_underscores=True,
                )

        if args.min_subpopulation_lineages_per_population > args.max_subpopulation_lineages_per_population:
            sys.exit("Minimum number of subpopulation lineages must be greater than maximum number of subpopulation lineages")
        elif args.min_subpopulation_lineages_per_population < 0:
            sys.exit("Number of subpopulation lineages cannot be negative")
        elif args.min_subpopulation_lineages_per_population == args.max_subpopulation_lineages_per_population:
            if args.min_subpopulation_lineages_per_population <= 1:
                num_subpopulation_lineages_per_population_fn = lambda: 1
            else:
                num_subpopulation_lineages_per_population_fn = lambda: args.min_subpopulation_lineages_per_population
        else:
            num_subpopulation_lineages_per_population_fn = lambda:rng.randint(args.min_subpopulation_lineages_per_population, args.max_subpopulation_lineages_per_population)

        source_tree.calc_node_ages()
        population_tree, gene_trees = generate_contained_trees(
                containing_tree=source_tree,
                num_subpopulation_lineages_per_population_fn=num_subpopulation_lineages_per_population_fn,
                subpopulation_root_age=args.subpopulation_root_age,
                num_individuals_per_terminal_lineage=args.num_individuals_per_terminal_lineage,
                num_gene_trees=args.num_loci_per_individual,
                population_size=args.population_size,
                rng=rng,
                )
        gene_trees.write(
                path="{}.gene-trees.nex".format(job_title),
                schema="nexus")
        population_tree.write(
                path="{}.guide-tree.nex".format(job_title),
                schema="nexus")

        imap_filepath = "{}.input.imap.txt".format(job_title)
        f = open(imap_filepath, "w")
        for taxon in gene_trees.taxon_namespace:
            f.write("{}    {}\n".format(taxon.label.split("^")[1], taxon.population_label))
            # f.write("{}    {}\n".format(taxon.label.split("^")[0], taxon.population_label))
            # f.write("{}    {}\n".format(taxon.label, taxon.population_label))
        f.write("\n")

        d0 = sg.generate(gene_trees)
        chars_filepath = "{}.input.chars.txt".format(job_title)
        f = open(chars_filepath, "w")
        for cm_idx, cm in enumerate(d0.char_matrices):
            try:
                td = popgenstat.tajimas_d(cm)
            except ZeroDivisionError as e:
                td = "N/A"
            # td = 1
            wtheta = popgenstat.wattersons_theta(cm)
            sys.stderr.write("Locus {}: pi = {}, Watterson's theta per site = {}, Tajima's D = {}\n".format(
                cm_idx+1,
                popgenstat.nucleotide_diversity(cm),
                wtheta/args.num_characters_per_locus,
                td))
            cm.write(file=f, schema="phylip")
            f.write("\n")

        out_filepath = "{}.results.out.txt".format(job_title)
        mcmc_filepath = "{}.results.mcmc.txt".format(job_title)
        # Inverse Gamma Prior
        # IG(a,b), with mean given by b/(a-1)
        # So,
        #   thetaprior 3 0.002
        # has a mean of
        #   0.002/(3-1) = 0.001
        theta_prior_mean = args.population_size * 4 * args.mutation_rate_per_site
        theta_prior_a = 3.0
        theta_prior_b = theta_prior_mean * (theta_prior_a - 1)
        if args.no_scale_tree_by_mutation_rate:
            tau_prior_mean = population_tree.seed_node.age
        else:
            # tau_prior_mean = population_tree.seed_node.age * args.population_size * 4 * args.mutation_rate_per_site
            tau_prior_mean = population_tree.seed_node.age * args.mutation_rate_per_site * (1.0 / (args.num_loci_per_individual * args.num_characters_per_locus))
            # tau_prior_mean = population_tree.seed_node.age / 100000
        tau_prior_a = 3.0
        tau_prior_b = tau_prior_mean * (tau_prior_a - 1)

        num_populations = len(population_tree.taxon_namespace)
        population_labels = " ".join(t.label for t in population_tree.taxon_namespace)
        num_individuals_per_population = " ".join(str(args.num_individuals_per_terminal_lineage) for i in range(len(population_tree.taxon_namespace)))
        # num_individuals_per_population = " ".join(str(t.num_individuals_sampled) for t in population_tree.taxon_namespace)
        num_input_lineages = len(population_labels)

        bpp_guide_tree = population_tree.as_string(
                schema="newick",
                suppress_leaf_taxon_labels=False,
                suppress_leaf_node_labels=True,
                suppress_internal_taxon_labels=True,
                suppress_internal_node_labels=True,
                suppress_rooting=True,
                suppress_edge_lengths=True,
                unquoted_underscores=True,
                preserve_spaces=True,
                store_tree_weights=False,
                suppress_annotations=True,
                suppress_item_comments=True,
                )
        bpp_config = BPP_TEMPLATE.format(
                chars_filepath=chars_filepath,
                imap_filepath=imap_filepath,
                out_filepath=out_filepath,
                mcmc_filepath=mcmc_filepath,
                num_populations=num_populations,
                population_labels=population_labels,
                num_individuals_per_population=num_individuals_per_population,
                bpp_guide_tree=bpp_guide_tree,
                theta_prior_mean=theta_prior_mean,
                theta_prior_a=theta_prior_a,
                theta_prior_b=theta_prior_b,
                tau_prior_mean=tau_prior_mean,
                tau_prior_a=tau_prior_a,
                tau_prior_b=tau_prior_b,
                num_loci=args.num_loci_per_individual,
                root_age=population_tree.seed_node.age
                )
        bpp_ctl_filepath = "{}.input.bpp.ctl".format(job_title)
        f = open(bpp_ctl_filepath, "w")
        f.write(bpp_config)
        f.write("\n")

        jobf = open("{}.job".format(job_title), "w")
        jobf.write("#! /bin/bash\n\n")
        jobf.write("set -e -o pipefail\n")
        if args.job_preamble:
            with open(os.path.expanduser(os.path.expandvars(args.job_preamble))) as job_preamble_src:
                jobf.write(job_preamble_src.read())
        # jobf.write("#$ -cwd\n")
        # jobf.write("#$ -V\n")
        # jobf.write("#$ -S /bin/bash\n")
        # jobf.write("#$ -l h_vmem=12G\n")
        # jobf.write("#$ -l virtual_free=12G\n")
        jobf.write("bpp --cfile {}\n".format(bpp_ctl_filepath))
        jobf.write("spdw-extract-bpp-a10-tree.py -p {population_probability_threshold} -o {job_title}.results.summary {job_title}.results.out.txt {job_title}.guide-tree.nex\n".format(
            population_probability_threshold=args.population_probability_threshold,
            job_title=job_title,
            ))

if __name__ == "__main__":
    main()


