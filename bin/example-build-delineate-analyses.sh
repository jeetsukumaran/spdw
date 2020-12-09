#! /bin/bash

set -e -o pipefail

python ../../bin/spdw-sim-protractedspeciation-trees.py -o pbd-trees --max-extant-lineages 12 --min-extant-orthospecies 5 -n 1
python ../../bin/spdw-delineate-create-analyses-from-trees.py --constrain-partitions random --num-unconstrained-leaves 6 pbd-trees001.lineage.tre
delineate-estimate partitions -t delineate-run.tree.nex -c delineate-run.tsv
python ../../bin/spdw-prune-delineate-result-trees-of-unconstrained-leaves.py --max 1 delineate-run.delimitation-results.trunc.trees
