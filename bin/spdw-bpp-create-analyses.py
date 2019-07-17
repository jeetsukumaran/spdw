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
from dendropy.utility import textprocessing
from spdw import spdwlib
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
Some notes about the inverse-gamma distribution. Since BPP 3.4, both the θ and τ
parameters are assigned the inverse gamma priors rather than the gamma priors in version 3.3
or earlier. One difference is that the gamma is light-tailed while the inverse-gamma is heavy-
tailed, so that the inverse-gamma may be less influential than the gamma if your prior mean
is much too small. The inverse-gamma distribution IG( α , β ) has mean m = β /( α – 1) if α > 1
and variance s 2 = β 2 /[( α – 1) 2 ( α – 2)] if α > 2, while the coefficient of variation is s/m =
1 ( α − 2) . If little information is available about the parameters, you can use α = 3 for a
diffuse prior and then adjust the β so that the mean looks reasonable. Both parameters θ s and
τ s in the MSC model are measured by genetic distance, the expected number of
mutations/substitutions per site. For example, for the human species, θ H ≈ 0.0006, which
means that two random sequences from the human population are different at ~0.06% of
sites, less than 1 difference per kb. A sensible diffuse prior is then “ thetaprior = 3 0.002 ”,
with mean 0.001.
tauprior = 3 0.03 specifies the inverse-gamma prior IG( α , β ) for τ 0 , the divergence time
parameter for the root in the species tree. Other divergence times are generated from the
uniform Dirichlet distribution (Yang and Rannala, 2010: equation 2). In the example, the
mean is 0.03/(3 – 1) = 0.015 (which means 1.5% of sequence divergence between the root of
the species tree and the present time). If the mutation rate is 10 –9 mutations/site/year, this
distance will translate to a human-orangutan divergence time of 15MY .
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

    *  - Root age = {root_age}; N = {population_size}; mu = {mutation_rate_per_site};
    *  - Optimized Priors:
    *      - thetaprior = {optimized_theta_prior_a} {optimized_theta_prior_b} * Mean = {optimized_theta_prior_mean}
    *      - tauprior = {optimized_tau_prior_a} {optimized_tau_prior_b} * Mean = {optimized_tau_prior_mean};
    *  - Default Priors:
    *      - thetaprior = {default_theta_prior_a} {default_theta_prior_b}
    *      - tauprior = {default_tau_prior_a} {default_tau_prior_b}
    thetaprior = {theta_prior_a} {theta_prior_b}
      tauprior = {tau_prior_a}   {tau_prior_b}

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
            help="Path to containing tree files. Specify '-' to read from standard input. Tree branch lengths should be in units of N. Specify '--population-size' to scale appropriately.")
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
    analysis_options = parser.add_argument_group("Analysis Options")
    parser.add_argument("--default-priors",
            dest="is_use_informative_priors",
            action="store_false",
            help="Use default priors instead of optimizing for (known) true value.")
    data_options = parser.add_argument_group("Data Options")
    data_options.add_argument("--population-size",
            type=float,
            default=1e6,
            help="Population size [default: %(default)s]. Branch lengths will be scaled (multiplied) by this factor.")
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
            default=0,
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
    parser.add_argument("--save-alignments",
            action="store_true",
            default=False,
            help="Save alignment data.")
    parser.add_argument("--save-locus-profile",
            action="store_true",
            default=False,
            help="Save locus profile information.")
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
    input_sequences_d = {}
    for idx in range(1):
        input_sequences_d["t{}".format(idx+1)] = [rng.choice(["A","C","G","T"]) for x in range(sg.seq_len)]
    input_sequences = dendropy.DnaCharacterMatrix.from_dict(input_sequences_d)
    sg.ancestral_seq_idx = 1

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

        population_tree, gene_trees = generate_contained_trees(
                containing_tree=source_tree,
                num_subpopulation_lineages_per_population_fn=num_subpopulation_lineages_per_population_fn,
                subpopulation_root_age=args.subpopulation_root_age,
                num_individuals_per_terminal_lineage=args.num_individuals_per_terminal_lineage,
                num_gene_trees=args.num_loci_per_individual,
                population_size=args.population_size,
                rng=rng,
                )
        population_tree.calc_node_ages()
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

        field_names = (
                "Locus",
                "pi",
                "theta",
                "TajD",
                "tau_mean",
                "tau_max",
                "tau_min",
                )
        table_row_formatter = spdwlib.TableRowFormatter(
                field_names=field_names,
                field_lengths=[5] + [11] * 6,
                field_value_templates = ["{}"] + ["{:<11.8f}"] * 6,
                field_separator = "    ",
                )
        d0 = sg.generate(
                gene_trees,
                input_sequences=input_sequences,
                )
        chars_filepath = "{}.input.chars.txt".format(job_title)
        f = open(chars_filepath, "w")
        _log("{}".format(table_row_formatter.format_header_row()))
        locus_thetas = []
        root_tip_max_thetas = []
        root_tip_mean_thetas = []
        locus_profile_info = []
        for cm_idx, cm in enumerate(d0.char_matrices):
            try:
                td = popgenstat.tajimas_d(cm)
            except ZeroDivisionError as e:
                td = "N/A"
            root_tip_thetas = []
            nucleotide_diversity = popgenstat.nucleotide_diversity(cm)
            # locus_thetas.append(nucleotide_diversity)
            wtheta = popgenstat.wattersons_theta(cm) / args.num_characters_per_locus
            locus_thetas.append(wtheta)
            for seq_idx, seq in enumerate(cm.values()):
                pair_data = (input_sequences[0], seq)
                num_segregating_sites = popgenstat._num_segregating_sites(
                        char_sequences=pair_data,
                        state_alphabet=cm.default_state_alphabet,
                        ignore_uncertain=True)

                a1 = sum([1.0/i for i in range(1, len(pair_data))])
                root_tip_theta = float(num_segregating_sites) / a1
                root_tip_thetas.append(root_tip_theta/sg.seq_len)

                # root_tip_theta = num_segregating_sites / sg.seq_len
                # root_tip_thetas.append(root_tip_theta)

            root_tip_max_thetas.append(max(root_tip_thetas))
            root_tip_mean_thetas.append(sum(root_tip_thetas)/len(root_tip_thetas))
            field_values = [
                cm_idx+1,
                nucleotide_diversity,
                locus_thetas[-1],
                td,
                root_tip_mean_thetas[-1],
                root_tip_max_thetas[-1],
                min(root_tip_thetas),
            ]
            assert len(field_values) == len(field_names)
            profile = table_row_formatter.format_row(field_values=field_values)
            _log("{}".format(profile))
            locus_profile_info.append(dict(zip(field_names, field_values)))
            cm.write(file=f, schema="phylip")
            if args.save_alignments:
                cm.write(path="{}.locus-{:04d}.nex".format(job_title, cm_idx+1), schema="nexus")
            f.write("\n")
        if args.save_locus_profile:
            with open("{}.locus-profile.tsv".format(job_title), "w") as locus_profile_f:
                csv_writer = csv.DictWriter(
                        locus_profile_f,
                        delimiter="\t",
                        fieldnames=field_names,)
                csv_writer.writeheader()
                csv_writer.writerows(locus_profile_info)

        out_filepath = "{}.results.out.txt".format(job_title)
        mcmc_filepath = "{}.results.mcmc.txt".format(job_title)

        # Inverse Gamma Prior
        # IG(a,b), with mean given by b/(a-1)
        # So,
        #   thetaprior 3 0.002
        # has a mean of
        #   0.002/(3-1) = 0.001
        # optimized_theta_prior_mean = args.population_size * 4 * args.mutation_rate_per_site
        optimized_theta_prior_mean = sum(locus_thetas) / len(locus_thetas)
        optimized_theta_prior_a = 3.0
        optimized_theta_prior_b = optimized_theta_prior_mean * (optimized_theta_prior_a - 1)

        optimized_tau_prior_mean = sum(root_tip_mean_thetas) / len(root_tip_mean_thetas)
        optimized_tau_prior_a = 3.0
        optimized_tau_prior_b = optimized_tau_prior_mean * (optimized_tau_prior_a - 1)

        default_theta_prior_a = 3.0
        default_theta_prior_b = 0.002
        default_theta_prior_mean =  default_theta_prior_b / (default_theta_prior_a - 1)
        default_tau_prior_a = 3.0
        default_tau_prior_b = 0.03
        default_tau_prior_mean =  default_tau_prior_b / (default_tau_prior_a - 1)

        if args.is_use_informative_priors:
            theta_prior_a = optimized_theta_prior_a
            theta_prior_b = optimized_theta_prior_b
            tau_prior_a = optimized_tau_prior_a
            tau_prior_b = optimized_tau_prior_b
        else:
            theta_prior_a = default_theta_prior_a
            theta_prior_b = default_theta_prior_b
            tau_prior_a = default_tau_prior_a
            tau_prior_b = default_tau_prior_b
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
                theta_prior_a=theta_prior_a,
                theta_prior_b=theta_prior_b,
                tau_prior_a=tau_prior_a,
                tau_prior_b=tau_prior_b,
                num_loci=args.num_loci_per_individual,
                root_age=population_tree.seed_node.age,
                population_size=args.population_size,
                mutation_rate_per_site=args.mutation_rate_per_site,
                num_loci_per_individual=args.num_loci_per_individual,
                optimized_theta_prior_a=optimized_theta_prior_a,
                optimized_theta_prior_b=optimized_theta_prior_b,
                optimized_theta_prior_mean=optimized_theta_prior_mean,
                default_theta_prior_a=default_theta_prior_a,
                default_theta_prior_b=default_theta_prior_b,
                default_theta_prior_mean=default_theta_prior_mean,
                optimized_tau_prior_a=optimized_tau_prior_a,
                optimized_tau_prior_b=optimized_tau_prior_b,
                optimized_tau_prior_mean=optimized_tau_prior_mean,
                default_tau_prior_a=default_tau_prior_a,
                default_tau_prior_b=default_tau_prior_b,
                default_tau_prior_mean=default_tau_prior_mean,
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


