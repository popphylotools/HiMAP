#!/usr/bin/env bash
conda env create -f environment.yml
conda env create -f mafft_env.yml
conda env create -f tapir_env.yml

# before running python scripts, run `source activate orthoExon`
