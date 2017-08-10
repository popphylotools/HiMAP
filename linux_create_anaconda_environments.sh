#!/usr/bin/env bash

# before running python scripts, run `source activate orthoExon`

conda env create -f environment.yml
conda env create -f mafft_env.yml
conda env create -f tapir_env.yml

# the tapir environment still needs hiphy2
source activate tapir

wget http://s3.faircloth-lab.org/packages/hyphy2.linux.gz
gunzip hyphy2.*.gz
chmod 0700 hyphy2.*
mv hyphy2.* $(dirname $(which python))/hyphy2

rm hyphy2*

source deactivate tapir