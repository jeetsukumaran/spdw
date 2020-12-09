#! /bin/bash

set -e -o pipefail

p ../../bin/spdw-sim-protractedspeciation-trees.py -o pbd-trees-30tips-05species --max-extant-lineages 30 --min-extant-orthospecies 5 -n 1
p ../../bin/spdw-delineate-create-analyses-from-trees.py --constrain-partitions random --num-unconstrained-leaves 20 pbd-trees-30tips-05species001.lineage.tre
