#!/bin/bash
shopt -s expand_aliases
module load mafft
module load parallel

rm "./data/intermediate/enhanced_scds_alignment/"*".fasta.aln"

"/bin/ls" "./data/intermediate/enhanced_scds_alignment/"*".fasta" | parallel "mafft --localpair --maxiterate 1000 {} > {}.aln"
