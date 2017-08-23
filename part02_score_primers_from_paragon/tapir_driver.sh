#!/usr/bin/env bash
source activate tapir

nex_sub_path=${1}
tapir_out_sub_path=${2}
ref_tree_fn=${3}
tapir_opts_string=${4}

tapir_compute.py ${nex_sub_path} ${ref_tree_fn} \
${tapir_opts_string}
--output ${tapir_out_sub_path}
