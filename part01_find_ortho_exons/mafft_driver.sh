#!/usr/bin/env bash

DIRECTORY=${1}

rm ${DIRECTORY}*".fasta.aln"

"/bin/ls" ${DIRECTORY}*".fasta" | parallel "mafft --localpair --maxiterate 1000 {} > {}.aln"
