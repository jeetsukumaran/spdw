#! /bin/bash

set -e -o pipefail

python bin/spdw_sim_coaltrees.py \
    -k 100 \
    --num-reps 100 \
    coal1
python bin/spdw_sim_bdstruct_coaltrees.py \
    -p50 \
    -k2 \
    --num-pop-tree-reps 1 \
    --num-reps-per-pop-tree 100 \
    bd1
python bin/spdw_plotcoaltimes.py *.trees
