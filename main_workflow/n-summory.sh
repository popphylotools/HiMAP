#!/bin/bash

N=${1:-"0"}

"./step11_summary.py" \
--summary_fn "../data/output${N}/summory.csv" \
--primer3_path "../data/intermediate${N}/primer_design/" \
--unpadded_primer_product_path "../data/output${N}/primerProducts/" \
--json_path "../data/intermediate${N}/json/" \
--pi_score_path "../data/intermediate${N}/phylo_informativeness/sql/"
