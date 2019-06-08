#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import os
import collections
import json
import sys
import argparse
import inspect
import itertools
import io
import dendropy

class Partition(object):

    @staticmethod
    def build_species_leafsets(d):
        if d is None:
            return None
        return frozenset(frozenset(s) for s in d)

    def __init__(self, **kwargs):
        self.part_id = kwargs.get("id", None)
        self.log_probability = kwargs.get("log_probability", None)
        self.probability = kwargs.get("probability", None)
        self.probability_given_constraints = kwargs.get("probability_given_constraints", None)
        self.log_probability_given_constraints = kwargs.get("log_probability_given_constraints", None)
        self.cumulative_probability = kwargs.get("cumulative_probability", None)
        self.cumulative_probability_given_constraints = kwargs.get("cumulative_probability_given_constraints", None)
        self.is_in_confidence_interval = kwargs.get("is_in_confidence_interval", None)
        self.species_leafsets = Partition.build_species_leafsets(kwargs.get("species_leafsets", None))

    def is_conspecific(self, lineage1, lineage2):
        for sp in self.species_leafsets:
            if lineage1 in sp and lineage2 in sp:
                return True
        return False

    def __len__(self):
        return len(self.species_leafsets)

class Evaluator(object):

    def __init__(self):
        pass

    def execute(self, args):
        self.config_path = args.configuration_filepath
        if not self.config_path:
            sys.exit("Path to configuration file needs to be specified")
        self.analysis_results_filepath = args.analysis_results_filepath
        if not self.analysis_results_filepath:
            sys.exit("Path to analysis results file needs to be specified")

        self.report_dest = sys.stdout
        self.load_data()
        perf_data_rows = []
        if args.is_evaluate_marginal:
            if args.lineage_tree_filepath is None:
                sys.exit("Require path to lineage tree filepath to analyze marginal probabilities")
            with open(os.path.expandvars(os.path.expanduser(args.lineage_tree_filepath))) as src:
                self.read_lineage_tree(src, schema="nexus")
            self.load_estimated_partitions()
            for lineage_pair in self.all_distinct_pairs_of_unconstrained_lineages():
                perf_data = collections.OrderedDict()
                self.store_basic_features(perf_data)
                perf_data.update(self.calc_lineage_pair_features(*lineage_pair))
                results = self.calc_marginal_probability_of_conspecificity(*lineage_pair)
                for key in [
                        "lineage_pair_is_true_conspecific",
                        "lineage_pair_conspecificity_marginal_probability",
                        "lineage_pair_conspecificity_marginal_probability_given_constraints",
                        "lineage_pair_nonconspecificity_marginal_probability",
                        "lineage_pair_nonconspecificity_marginal_probability_given_constraints",
                        ]:
                    perf_data[key] = results[key]
                perf_data_rows.append(perf_data)
        else:
            perf_data = collections.OrderedDict()
            self.store_basic_features(perf_data)
            self.standard_performance_assessment(perf_data)
            perf_data_rows.append(perf_data)
        assert perf_data_rows
        self.report(perf_data_rows)

    def load_data(self):
        with open(self.config_path) as src:
            self.config = json.load(src)
        with open(self.analysis_results_filepath) as src:
            self.estimation_results = json.load(src)
        self.set_true_partition(species_leafsets=self.config["test_info"]["true_species_leafsets"])

    def set_true_partition(self, species_leafsets):
        self.true_partition = Partition(species_leafsets=species_leafsets)

    def load_estimated_partitions(self):
        self.partitions = [Partition(**p) for p in self.estimation_results["partitions"]]
        return self.partitions

    def read_lineage_tree(self, src, schema="nexus"):
        self.set_lineage_tree(
                file=src,
                schema=schema,
                rooting="force-rooted")

    def set_lineage_tree(self, **kwargs):
        self.lineage_tree = dendropy.Tree.get(**kwargs)
        self.lineage_tree.encode_bipartitions()
        self.lineage_tree.calc_node_ages()
        # self.lineage_tree_label_node_map = {taxon.label:taxon for taxon in self.lineage_tree.taxon_namespace}
        self.lineage_tree_label_node_map = {nd.taxon.label:nd for nd in self.lineage_tree.leaf_node_iter()}
        self.phylogenetic_distance_matrix = self.lineage_tree.phylogenetic_distance_matrix(is_store_path_edges=False)

    def all_distinct_pairs_of_unconstrained_lineages(self):
        unconstrained_lineages = self.config["test_info"]["unconstrained_lineages"]
        for x in itertools.combinations(unconstrained_lineages, 2):
            yield x

    def is_span_root(self, lineage1, lineage2):
        n1 = self.lineage_tree_label_node_map[lineage1]
        n2 = self.lineage_tree_label_node_map[lineage2]
        assert n1 is not n2
        r_left, r_right = self.lineage_tree.seed_node._child_nodes
        if n1.bipartition.is_nested_within(r_left.bipartition):
            return n2.bipartition.is_nested_within(r_right.bipartition)
        else:
            return n2.bipartition.is_nested_within(r_left.bipartition)

    def calc_lineage_pair_features(self, lineage1, lineage2):
        result = {}
        result["lineage_pair_unnormalized_weighted_distance"] = self.phylogenetic_distance_matrix.distance(
                self.lineage_tree_label_node_map[lineage1].taxon,
                self.lineage_tree_label_node_map[lineage2].taxon,
                is_weighted_edge_distances=True,
                is_normalize_by_tree_size=False)
        result["lineage_pair_normalized_weighted_distance"] = self.phylogenetic_distance_matrix.distance(
                self.lineage_tree_label_node_map[lineage1].taxon,
                self.lineage_tree_label_node_map[lineage2].taxon,
                is_weighted_edge_distances=True,
                is_normalize_by_tree_size=True)
        result["lineage_pair_unnormalized_unweighted_distance"] = self.phylogenetic_distance_matrix.distance(
                self.lineage_tree_label_node_map[lineage1].taxon,
                self.lineage_tree_label_node_map[lineage2].taxon,
                is_weighted_edge_distances=False,
                is_normalize_by_tree_size=False)
        result["lineage_pair_normalized_unweighted_distance"] = self.phylogenetic_distance_matrix.distance(
                self.lineage_tree_label_node_map[lineage1].taxon,
                self.lineage_tree_label_node_map[lineage2].taxon,
                is_weighted_edge_distances=False,
                is_normalize_by_tree_size=True)
        mrca = self.phylogenetic_distance_matrix.mrca(
                self.lineage_tree_label_node_map[lineage1].taxon,
                self.lineage_tree_label_node_map[lineage2].taxon)
        result["lineage_pair_mrca_age"] = mrca.age
        result["lineage_pair_is_span_root"] = self.is_span_root(lineage1, lineage2)
        return result

    def calc_marginal_probability_of_conspecificity(self, lineage1, lineage2):
        results = {
            "lineage_pair_conspecificity_marginal_probability": 0.0,
            "lineage_pair_conspecificity_marginal_probability_given_constraints": 0.0,
            "lineage_pair_nonconspecificity_marginal_probability": 0.0,
            "lineage_pair_nonconspecificity_marginal_probability_given_constraints": 0.0,
            "lineage_pair_conspecific_partitions": [],
        }
        for partition in self.partitions:
            if partition.is_conspecific(lineage1, lineage2):
                results["lineage_pair_conspecific_partitions"].append(partition)
                results["lineage_pair_conspecificity_marginal_probability"] += partition.probability
                results["lineage_pair_conspecificity_marginal_probability_given_constraints"] += partition.probability_given_constraints
            else:
                results["lineage_pair_nonconspecificity_marginal_probability"] += partition.probability
                results["lineage_pair_nonconspecificity_marginal_probability_given_constraints"] += partition.probability_given_constraints
        results["lineage_pair_is_true_conspecific"] = self.true_partition.is_conspecific(lineage1, lineage2)
        return results

    def store_basic_features(self, perf_data):
        perf_data["batch_id"] = self.estimation_results["batch_id"]
        perf_data["root_age"] = self.estimation_results["root_age"]
        perf_data["num_tips"] = self.estimation_results["num_tips"]
        perf_data["total_num_partitions"] = self.estimation_results["num_partitions"]
        perf_data["true_speciation_completion_rate"] = self.config["test_info"]["true_speciation_completion_rate"]
        perf_data["true_num_species"] = len(self.true_partition)
        perf_data["num_constrained_species"] = self.config["test_info"]["species_partition_estimation_num_constrained_species"] # number of species defined (may not include all lineages)
        perf_data["num_constrained_lineages"] = self.config["test_info"]["species_partition_estimation_num_constrained_lineages"] # number of lineages assigned to species
        perf_data["num_unconstrained_lineages"] = self.config["test_info"]["species_partition_estimation_num_unconstrained_lineages"] # number of lineages not assigned to species
        perf_data["estimated_speciation_completion_rate"] = self.estimation_results["speciation_completion_rate"]
        perf_data["speciation_completion_rate_source"] = self.estimation_results["speciation_completion_rate_source"]

    def standard_performance_assessment(self, perf_data):
        perf_data["total_num_partitions_in_confidence_interval"] = self.estimation_results["num_partitions_in_confidence_interval"]
        perf_data["inferred_partition_num_species"] = len(self.estimation_results["partitions"][0]["species_leafsets"])
        # perf_data["inferred_partition_log_probability"] = self.estimation_results["partitions"][0]["log_probability"]
        perf_data["inferred_partition_probability"] = self.estimation_results["partitions"][0]["probability"]
        # perf_data["inferred_partition_log_probability_given_constraints"] = self.estimation_results["partitions"][0]["log_probability_given_constraints"]
        perf_data["inferred_partition_probability_given_constraints"] = self.estimation_results["partitions"][0]["probability_given_constraints"]
        for partition_idx, partition_info in enumerate(self.estimation_results["partitions"]):
            current_partition = Partition(**partition_info)
            if current_partition.species_leafsets == self.true_partition.species_leafsets:
                # perf_data["true_partition_log_probability"] = current_partition.log_probability
                perf_data["true_partition_probability"] = current_partition.probability
                perf_data["true_partition_cumulative_probability"] = current_partition.cumulative_probability
                # perf_data["true_partition_log_probability_given_constraints"] = current_partition.log_probability_given_constraints
                perf_data["true_partition_probability_given_constraints"] = current_partition.probability_given_constraints
                perf_data["true_partition_cumulative_probability_given_constraints"] = current_partition.cumulative_probability_given_constraints
                if partition_idx == 0:
                    perf_data["is_true_partition_preferred"] = True
                else:
                    perf_data["is_true_partition_preferred"] = False
                perf_data["is_true_partition_in_confidence_interval"] = current_partition.is_in_confidence_interval
                break
        else:
            raise ValueError("True partition not found in results")
        return perf_data

    def report(self, perf_data_rows):
        # json.dump(perf_data, sys.stdout, indent=4, separators=(',', ': '))
        delimiter = "\t"
        self.report_dest.write(delimiter.join(perf_data_rows[0].keys()))
        self.report_dest.write("\n")
        for perf_data in perf_data_rows:
            value_row = []
            for idx, v in enumerate(perf_data.values()):
                if isinstance(v, bool):
                    value_row.append(str(v).upper()) # for R
                else:
                    value_row.append(str(v))
            self.report_dest.write(delimiter.join(value_row))
            self.report_dest.write("\n")

class TestRunner(object):

    def __init__(self):
        self.test_log = lambda msg: sys.stdout.write("-[{}]: {}\n".format(inspect.stack()[1][3], msg))
        self.test_data_dir = os.path.join(os.path.abspath(__file__), "_test_data")

    def run_tests(self):
        self.test_is_conspecific()
        self.test_marginal_probability_of_conspecificity()
        self.test_is_span_root()
        self.test_lineage_pair_distances()
        self.test_all_distinct_pairs_of_unconstrained_lineages()

    def test_is_conspecific(self):
        d = {
            "species_leafsets": [
                    ["a", "b", "c"],
                    ["d", "e", "f"],
                    ["g"],
                    ["h"],
                ],
        }
        partition = Partition(**d)
        assert     partition.is_conspecific("a", "b")
        assert     partition.is_conspecific("a", "c")
        assert not partition.is_conspecific("a", "d")
        assert not partition.is_conspecific("a", "e")
        assert not partition.is_conspecific("a", "f")
        assert not partition.is_conspecific("a", "g")
        assert not partition.is_conspecific("a", "h")

        assert     partition.is_conspecific("b", "a")
        assert     partition.is_conspecific("b", "c")
        assert not partition.is_conspecific("b", "d")
        assert not partition.is_conspecific("b", "e")
        assert not partition.is_conspecific("b", "f")
        assert not partition.is_conspecific("b", "g")
        assert not partition.is_conspecific("b", "h")

        assert     partition.is_conspecific("c", "a")
        assert     partition.is_conspecific("c", "b")
        assert not partition.is_conspecific("c", "d")
        assert not partition.is_conspecific("c", "e")
        assert not partition.is_conspecific("c", "f")
        assert not partition.is_conspecific("c", "g")
        assert not partition.is_conspecific("c", "h")

        assert not partition.is_conspecific("d", "a")
        assert not partition.is_conspecific("d", "b")
        assert not partition.is_conspecific("d", "c")
        assert     partition.is_conspecific("d", "e")
        assert     partition.is_conspecific("d", "f")
        assert not partition.is_conspecific("d", "g")
        assert not partition.is_conspecific("d", "h")

        assert not partition.is_conspecific("e", "a")
        assert not partition.is_conspecific("e", "b")
        assert not partition.is_conspecific("e", "c")
        assert     partition.is_conspecific("e", "d")
        assert     partition.is_conspecific("e", "f")
        assert not partition.is_conspecific("e", "g")
        assert not partition.is_conspecific("e", "h")

        assert not partition.is_conspecific("f", "a")
        assert not partition.is_conspecific("f", "b")
        assert not partition.is_conspecific("f", "c")
        assert     partition.is_conspecific("f", "d")
        assert     partition.is_conspecific("f", "e")
        assert not partition.is_conspecific("f", "g")
        assert not partition.is_conspecific("f", "h")

        assert not partition.is_conspecific("g", "a")
        assert not partition.is_conspecific("g", "b")
        assert not partition.is_conspecific("g", "c")
        assert not partition.is_conspecific("g", "d")
        assert not partition.is_conspecific("g", "e")
        assert not partition.is_conspecific("g", "f")
        assert not partition.is_conspecific("g", "h")

        assert not partition.is_conspecific("h", "a")
        assert not partition.is_conspecific("h", "b")
        assert not partition.is_conspecific("h", "c")
        assert not partition.is_conspecific("h", "d")
        assert not partition.is_conspecific("h", "e")
        assert not partition.is_conspecific("h", "f")
        assert not partition.is_conspecific("h", "g")

        self.test_log("OK")

    def test_is_span_root(self):
        tree_src = io.StringIO("(((a:1.25, b:1.25):1.25, c:2.5):1.5, (d:2.25, (e:0.5,f:0.5):1.75):1.75):2.5;")
        ev = Evaluator()
        ev.read_lineage_tree(src=tree_src, schema="newick")
        expected = {
                "ab": False,
                "ac": False,
                "ad": True,
                "ae": True,
                "af": True,
                "ba": False,
                "bc": False,
                "bd": True,
                "be": True,
                "bf": True,
                "ca": False,
                "cb": False,
                "cd": True,
                "ce": True,
                "cf": True,
                "da": True,
                "db": True,
                "dc": True,
                "de": False,
                "df": False,
                "ea": True,
                "eb": True,
                "ec": True,
                "ed": False,
                "ef": False,
                "fa": True,
                "fb": True,
                "fc": True,
                "fd": False,
                "fe": False,
        }
        for k, val in expected.items():
            assert ev.is_span_root(k[0], k[1]) is val
        self.test_log("OK")

    def test_lineage_pair_distances(self):
        ev = Evaluator()
        tree_src = io.StringIO("(((a:1.25, b:1.25):1.25, c:2.5):1.5, (d:2.25, (e:0.5,f:0.5):1.75):1.75):2.5;")
        ev.read_lineage_tree(src=tree_src, schema="newick")
        expected = {
            frozenset({'a', 'b'}): {'lineage_pair_unnormalized_weighted_distance': 2.5, 'lineage_pair_normalized_weighted_distance': 0.14705882352941177,  'lineage_pair_unnormalized_unweighted_distance': 2, 'lineage_pair_normalized_unweighted_distance': 0.18181818181818182, 'lineage_pair_mrca_age': 1.25, 'lineage_pair_is_span_root': False},
            frozenset({'a', 'c'}): {'lineage_pair_unnormalized_weighted_distance': 5.0, 'lineage_pair_normalized_weighted_distance': 0.29411764705882354,  'lineage_pair_unnormalized_unweighted_distance': 3, 'lineage_pair_normalized_unweighted_distance': 0.2727272727272727,  'lineage_pair_mrca_age': 2.5,  'lineage_pair_is_span_root': False},
            frozenset({'a', 'd'}): {'lineage_pair_unnormalized_weighted_distance': 8.0, 'lineage_pair_normalized_weighted_distance': 0.47058823529411764,  'lineage_pair_unnormalized_unweighted_distance': 5, 'lineage_pair_normalized_unweighted_distance': 0.45454545454545453, 'lineage_pair_mrca_age': 4.0,  'lineage_pair_is_span_root': True},
            frozenset({'a', 'e'}): {'lineage_pair_unnormalized_weighted_distance': 8.0, 'lineage_pair_normalized_weighted_distance': 0.47058823529411764,  'lineage_pair_unnormalized_unweighted_distance': 6, 'lineage_pair_normalized_unweighted_distance': 0.5454545454545454,  'lineage_pair_mrca_age': 4.0,  'lineage_pair_is_span_root': True},
            frozenset({'a', 'f'}): {'lineage_pair_unnormalized_weighted_distance': 8.0, 'lineage_pair_normalized_weighted_distance': 0.47058823529411764,  'lineage_pair_unnormalized_unweighted_distance': 6, 'lineage_pair_normalized_unweighted_distance': 0.5454545454545454,  'lineage_pair_mrca_age': 4.0,  'lineage_pair_is_span_root': True},
            frozenset({'b', 'c'}): {'lineage_pair_unnormalized_weighted_distance': 5.0, 'lineage_pair_normalized_weighted_distance': 0.29411764705882354,  'lineage_pair_unnormalized_unweighted_distance': 3, 'lineage_pair_normalized_unweighted_distance': 0.2727272727272727,  'lineage_pair_mrca_age': 2.5,  'lineage_pair_is_span_root': False},
            frozenset({'b', 'd'}): {'lineage_pair_unnormalized_weighted_distance': 8.0, 'lineage_pair_normalized_weighted_distance': 0.47058823529411764,  'lineage_pair_unnormalized_unweighted_distance': 5, 'lineage_pair_normalized_unweighted_distance': 0.45454545454545453, 'lineage_pair_mrca_age': 4.0,  'lineage_pair_is_span_root': True},
            frozenset({'b', 'e'}): {'lineage_pair_unnormalized_weighted_distance': 8.0, 'lineage_pair_normalized_weighted_distance': 0.47058823529411764,  'lineage_pair_unnormalized_unweighted_distance': 6, 'lineage_pair_normalized_unweighted_distance': 0.5454545454545454,  'lineage_pair_mrca_age': 4.0,  'lineage_pair_is_span_root': True},
            frozenset({'b', 'f'}): {'lineage_pair_unnormalized_weighted_distance': 8.0, 'lineage_pair_normalized_weighted_distance': 0.47058823529411764,  'lineage_pair_unnormalized_unweighted_distance': 6, 'lineage_pair_normalized_unweighted_distance': 0.5454545454545454,  'lineage_pair_mrca_age': 4.0,  'lineage_pair_is_span_root': True},
            frozenset({'c', 'd'}): {'lineage_pair_unnormalized_weighted_distance': 8.0, 'lineage_pair_normalized_weighted_distance': 0.47058823529411764,  'lineage_pair_unnormalized_unweighted_distance': 4, 'lineage_pair_normalized_unweighted_distance': 0.36363636363636365, 'lineage_pair_mrca_age': 4.0,  'lineage_pair_is_span_root': True},
            frozenset({'c', 'e'}): {'lineage_pair_unnormalized_weighted_distance': 8.0, 'lineage_pair_normalized_weighted_distance': 0.47058823529411764,  'lineage_pair_unnormalized_unweighted_distance': 5, 'lineage_pair_normalized_unweighted_distance': 0.45454545454545453, 'lineage_pair_mrca_age': 4.0,  'lineage_pair_is_span_root': True},
            frozenset({'f', 'c'}): {'lineage_pair_unnormalized_weighted_distance': 8.0, 'lineage_pair_normalized_weighted_distance': 0.47058823529411764,  'lineage_pair_unnormalized_unweighted_distance': 5, 'lineage_pair_normalized_unweighted_distance': 0.45454545454545453, 'lineage_pair_mrca_age': 4.0,  'lineage_pair_is_span_root': True},
            frozenset({'e', 'd'}): {'lineage_pair_unnormalized_weighted_distance': 4.5, 'lineage_pair_normalized_weighted_distance': 0.2647058823529412,   'lineage_pair_unnormalized_unweighted_distance': 3, 'lineage_pair_normalized_unweighted_distance': 0.2727272727272727,  'lineage_pair_mrca_age': 2.25, 'lineage_pair_is_span_root': False},
            frozenset({'f', 'd'}): {'lineage_pair_unnormalized_weighted_distance': 4.5, 'lineage_pair_normalized_weighted_distance': 0.2647058823529412,   'lineage_pair_unnormalized_unweighted_distance': 3, 'lineage_pair_normalized_unweighted_distance': 0.2727272727272727,  'lineage_pair_mrca_age': 2.25, 'lineage_pair_is_span_root': False},
            frozenset({'f', 'e'}): {'lineage_pair_unnormalized_weighted_distance': 1.0, 'lineage_pair_normalized_weighted_distance': 0.058823529411764705, 'lineage_pair_unnormalized_unweighted_distance': 2, 'lineage_pair_normalized_unweighted_distance': 0.18181818181818182, 'lineage_pair_mrca_age': 0.5,  'lineage_pair_is_span_root': False},
        }
        for lineage1, lineage2 in itertools.combinations("abcdef", 2):
            d = ev.calc_lineage_pair_features(lineage1, lineage2)
            key = frozenset([lineage1,lineage2])
            # print("{}: {}".format(key, d))
            for field in d:
                assert expected[key][field] == d[field]
        self.test_log("OK")

    def test_marginal_probability_of_conspecificity(self):
        results_d = {
            "partitions": [
                { "id": 0, "conspecifics": set(["ab", "cd"]), "species_leafsets": [["a", "b"], ["c", "d"], ["e"]],  "probability":  2, "probability_given_constraints": 3, },   # ab, cd
                { "id": 1, "conspecifics": set(["ab",     ]), "species_leafsets": [["a", "b", "c",], ["d"], ["e"]], "probability":  4, "probability_given_constraints": 5, },   # ab
                { "id": 2, "conspecifics": set(["ab", "cd"]), "species_leafsets": [["a", "b", "c", "d"], ["e"]],    "probability":  6, "probability_given_constraints": 7, },   # ab, cd
                { "id": 3, "conspecifics": set(["ab", "cd"]), "species_leafsets": [["a", "b", "e"], ["c", "d"]],    "probability":  8, "probability_given_constraints": 9, },   # ab, cd
                { "id": 4, "conspecifics": set(["ab", "cd"]), "species_leafsets": [["a", "b", "c", "d"], ["e"]],    "probability": 10, "probability_given_constraints": 11, },  # ab, cd
                { "id": 5, "conspecifics": set([      "cd"]), "species_leafsets": [["a", "e"], ["c", "d"], ["b"]],  "probability": 12, "probability_given_constraints": 13, },  # cd
                { "id": 6, "conspecifics": set([      "cd"]), "species_leafsets": [["a"], ["b", "c", "d"], ["e"]],  "probability": 14, "probability_given_constraints": 15, },  # cd
                { "id": 7, "conspecifics": set([          ]), "species_leafsets": [["a", "e", "d"], ["b", "c"]],    "probability": 16, "probability_given_constraints": 17, },
                { "id": 8, "conspecifics": set([          ]), "species_leafsets": [["a", "c"], ["b", "d"], ["e"]],  "probability": 18, "probability_given_constraints": 19, },
                { "id": 9, "conspecifics": set(["ab"      ]), "species_leafsets": [["b", "d"], ["b", "c", "a"]],    "probability": 20, "probability_given_constraints": 21, },  # ab
            ]
        }
        id_partition_map = { p["id"]: p for p in results_d["partitions"] }
        pair_keys = [
                "ab",
                "cd",
                "ec",
                ]
        ev = Evaluator()
        ev.estimation_results = results_d
        # true conspecifics out of pair keys == "ab" only
        # ev.set_true_partition([["a", "b", "c",], ["d"], ["e"]])
        ev.set_true_partition(id_partition_map[1]["species_leafsets"])
        true_conspecifics = set(["ab",])
        ev.load_estimated_partitions()
        for pk in pair_keys:
            partition_ids = [pid for pid in id_partition_map if pk in id_partition_map[pid]["conspecifics"]]
            marginal_probability = sum([ id_partition_map[pid]["probability"] for pid in partition_ids ])
            marginal_probability_given_constraints = sum([ id_partition_map[pid]["probability_given_constraints"] for pid in partition_ids ])
            r = ev.calc_marginal_probability_of_conspecificity(pk[0], pk[1])
            obs_pids = [p.part_id for p in r["lineage_pair_conspecific_partitions"]]
            assert set(obs_pids) == set(partition_ids)
            assert abs(marginal_probability - r["lineage_pair_conspecificity_marginal_probability"]) <= 1e-6
            assert abs(marginal_probability_given_constraints - r["lineage_pair_conspecificity_marginal_probability_given_constraints"]) <= 1e-6
            assert (pk in true_conspecifics) is (r["lineage_pair_is_true_conspecific"])

        self.test_log("OK")

    def test_all_distinct_pairs_of_unconstrained_lineages(self):
        config_d = {
                "test_info": {
                    "unconstrained_lineages": [
                        "S1.25",
                        "S1.29",
                        "S1.31",
                        "S1.48",
                        "S2.39"
                    ],
                }
        }
        expected = set({
            frozenset({'S1.25', 'S1.29'}),
            frozenset({'S1.31', 'S1.25'}),
            frozenset({'S1.25', 'S1.48'}),
            frozenset({'S2.39', 'S1.25'}),
            frozenset({'S1.31', 'S1.29'}),
            frozenset({'S1.29', 'S1.48'}),
            frozenset({'S2.39', 'S1.29'}),
            frozenset({'S1.31', 'S1.48'}),
            frozenset({'S2.39', 'S1.31'}),
            frozenset({'S2.39', 'S1.48'}),
            })
        ev = Evaluator()
        ev.config = config_d
        observed = set()
        for p in ev.all_distinct_pairs_of_unconstrained_lineages():
            observed.add(frozenset(p))
        assert expected == observed
        self.test_log("OK")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "configuration_filepath",
            metavar="CONFIGURATION-FILEPATH",
            nargs="?",
            help="Path to analysis configuration file (JSON format).")
    parser.add_argument(
            "analysis_results_filepath",
            metavar="ANALYSIS-RESULTS-FILEPATH",
            nargs="?",
            help="Path to analysis results file (JSON format).")
    parser.add_argument(
            "-m", "--marginal",
            action="store_true",
            dest="is_evaluate_marginal",
            default=False,
            help="Evaluate marginal probabilities of conspecificity of lineage pairs.")
    parser.add_argument(
            "-t", "--lineage-tree-filepath",
            metavar="TREE-FILEPATH",
            default=None,
            help="Path to file with lineage tree.")
    parser.add_argument(
            "--test",
            action="store_true",
            dest="is_run_tests",
            default=False,
            help="Run tests.")
    args = parser.parse_args()
    if args.is_run_tests:
        test_runner = TestRunner()
        test_runner.run_tests()
    else:
        ev = Evaluator()
        ev.execute(args)

if __name__ == "__main__":
    main()
