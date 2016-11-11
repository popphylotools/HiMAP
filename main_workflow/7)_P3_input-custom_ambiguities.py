#!/usr/bin/env python
# coding: utf-8

import os
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import shutil
import argparse
from config import full_species_list, species_list, transvestigated_species_set, primer3_path, \
    orthoCds_path

parser = argparse.ArgumentParser(description='This script creates primer3 input files')
parser.add_argument('--orthoCds_path', help='orthoCds_path', default=orthoCds_path)
parser.add_argument('--primer3_path',help='primer3_path', default=primer3_path)
parser.add_argument('-n','--ns_allowed',help="the number of n's allowed in primers", default='0')

args = parser.parse_args()

orthoCds_path = args.orthoCds_path
primer3_path = args.primer3_path

collapse_iupac = {
    ('-',): '-',
    ('A',): 'A',
    ('G',): 'G',
    ('C',): 'C',
    ('T',): 'T',
    ('C', 'T'): 'y',
    ('A', 'G'): 'r',
    ('A', 'T'): 'w',
    ('C', 'G'): 's',
    ('G', 'T'): 'k',
    ('A', 'C'): 'm',
    ('A', 'G', 'T'): 'd',
    ('A', 'C', 'G'): 'v',
    ('A', 'C', 'T'): 'h',
    ('C', 'G', 'T'): 'b',
    ('A', 'C', 'G', 'T'): 'n',
}

expand_iupac = {value.upper(): set(key) for key, value in collapse_iupac.items()}


def Consensus(aligned_seq_list):
    consensus = ""
    length = len(aligned_seq_list[0])
    for seq in aligned_seq_list:
        assert len(seq) == length
    for loc in range(length):
        base_set = set.union(*[expand_iupac[seq[loc].upper()] for seq in aligned_seq_list])
        if '-' in base_set:
            consensus += '-'
        else:
            consensus += collapse_iupac[tuple(sorted(base_set))]
    return consensus


# create handles for all .fasta files in fasta directory
fasta_fn = {name.split('.13spp.fasta')[0]: orthoCds_path + name for name in
            os.listdir(orthoCds_path) if
            ((".13spp.fasta" in name) and (".13spp.fasta.fai" not in name))}

# read and parse fasta files for each species
fasta = {}
for ortho in fasta_fn.keys():
    fasta[ortho] = {seq_record.id: seq_record
                    for seq_record in SeqIO.parse(fasta_fn[ortho],
                                                  "fasta", alphabet=IUPAC.ambiguous_dna)}

from Bio import motifs

fasta_degenerate = {}
for ortho in fasta:
    seq = Consensus([fasta[ortho][sp].upper().seq for sp in fasta[ortho].keys()])
    fasta_degenerate[ortho] = seq

# output
primer_product_size_range = '200-10000'
primer_thermodynamic_parameters_path = '/data0/opt/Primer3/primer3-2.3.6/src/primer3_config/'
primer_max_ns_accepted = args.ns_allowed
primer_liberal_base = '1'
shutil.rmtree(primer3_path)
os.makedirs(primer3_path, exist_ok=True)
for ortho in fasta.keys():
    filename = primer3_path + ortho + ".degenerate.p3"
    with open(filename, "w") as f:
        sequence_id = ortho
        sequence_template = fasta_degenerate[ortho]
        f.write(
            "SEQUENCE_ID={}\n"
            "SEQUENCE_TEMPLATE={}\n"
            "PRIMER_PRODUCT_SIZE_RANGE={}\n"
            "PRIMER_THERMODYNAMIC_PARAMETERS_PATH={}\n"
            "PRIMER_MAX_NS_ACCEPTED={}\n"
            "PRIMER_LIBERAL_BASE={}\n"
            "PRIMER_MIN_TM=53\n"
            "PRIMER_OPT_TM=56\n"
            "PRIMER_MAX_TM=60\n"
            "=".format(
                sequence_id,
                sequence_template,
                primer_product_size_range,
                primer_thermodynamic_parameters_path,
                primer_max_ns_accepted,
                primer_liberal_base))
