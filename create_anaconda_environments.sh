#!/usr/bin/env bash
conda env create -f environment.yml
conda env create -f part01_find_ortho_exons/mafft_env.yml
conda env create -f part02_score_primers_from_paragon/tapir_env.yml

# before running python scripts, run `source activate orthoExon`
