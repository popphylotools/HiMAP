#!/bin/bash
shopt -s expand_aliases
module load mafft
module load parallel
"/bin/ls" "../data/intermediate/template_scds_alignment/"*".fasta" | parallel "mafft --localpair --maxiterate 1000 {} > {}.aln"
