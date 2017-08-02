#!/usr/bin/env bash
source activate tapir

nex_sub_path = ${1}
tapir_out_sub_path = ${2}
ref_tree_fn = ${3}

tapir_compute.py ${nex_sub_path} ${ref_tree_fn} \
--intervals=2-3,3-5,5-8,8-10,10-18 \
--times=2,3,5,8,10,17 \
--tree-format=newick \
--output ${tapir_out_sub_path}
