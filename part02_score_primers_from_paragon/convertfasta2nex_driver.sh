#!/usr/bin/env bash

padded_primer_product_fn=${1}
nex_fn=${2}

perl ./convertfasta2nex.pl ${padded_primer_product_fn} > ${nex_fn}
