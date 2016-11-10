#!/usr/bin/env python
# coding: utf-8

import os
import gffutils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import json
from pyfaidx import Fasta
from Bio import SeqIO
import re

from config import full_species_list, species_list, transvestigated_species_set, \
    padded_primer_product_path, unpadded_primer_product_path, orthoCds_path, primer3_path, \
    json_path

# full_species_list = ['Bjar', 'Aobl', 'Bmin', 'Asus', 'Btry', 'Afra', 'Blat', 'Bzon', 'Bcor', 'Ccap', 'Bcur', 'Bole', 'Bdor']
# species_list = ["Bcur", "Bdor", "Bole", "Ccap"]
# transvestigated_species_set = {'Bcor', 'Blat', 'Bzon', 'Afra', 'Bmin', 'Bjar', 'Aobl'}
# padded_primer_product_path = "../intermediate/phylo_informativeness/fasta/"
# unpadded_primer_product_path = "../output/primerProducts/"
# orthoCds_path = "../output/orthoCds/"
# primer3_path = "../intermediate/primer_design/"
# json_path = "../intermediate/json/"

with open(json_path + "alternate_sp.json", 'r') as f:
    alternate_sp = json.load(f)

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

primer = {}
for p3_out_fn in (fn for fn in os.listdir(primer3_path) if ".p3.out" in fn):
    ortho = p3_out_fn.split('.degenerate.p3.out')[0]
    with open(primer3_path + p3_out_fn, 'r') as f:
        lines = f.readlines()
        lines = [line.strip().split('=') for line in lines]
        lines = {key: value for key, value in lines if key is not ''}
        if 'PRIMER_PAIR_NUM_RETURNED' not in lines.keys():
            continue
        if lines['PRIMER_PAIR_NUM_RETURNED'] is not '0':
            test = lines
            left, l_len = lines['PRIMER_LEFT_0'].split(',')
            right, r_len = lines['PRIMER_RIGHT_0'].split(',')
            start = int(left) + int(l_len)
            end = int(right) - int(r_len) + 1
            primer[ortho] = (start, end)

from copy import deepcopy

padded_fasta = {}
trimmed_fasta = {}
for ortho in fasta.keys():
    if ortho in primer.keys():
        start, end = primer[ortho]
    else:
        continue
    padding = {}
    for sp in full_species_list:
        if sp not in fasta[ortho].keys():
            for alt_sp in alternate_sp[sp]:
                if alt_sp in fasta[ortho].keys():
                    seq = fasta[ortho][alt_sp].seq[start:end]
                    des = "PADDING"
                    padding[sp] = SeqRecord(seq, id=sp, description=des)
                    break
    trimmed_fasta[ortho] = {sp: fasta[ortho][sp][start:end] for sp in fasta[ortho]}
    padded_fasta[ortho] = padding
    padded_fasta[ortho].update(trimmed_fasta[ortho])

for ortho in padded_fasta:
    for sp in padded_fasta[ortho]:
        padded_fasta[ortho][sp].description = ""

sp_order = {'Bcur': 1,
            'Bdor': 2,
            'Bole': 3,
            'Ccap': 4,
            'Bcor': 5,
            'Blat': 6,
            'Bzon': 7,
            'Afra': 8,
            'Bmin': 9,
            'Bjar': 10,
            'Aobl': 11,
            'Asus': 12,
            'Btry': 13}

# output fasta to pre_padding_species.json
os.makedirs(json_path, exist_ok=True)
filename = json_path + "pre_padding_species.json"
with open(filename, 'w') as f:
    json.dump({ortho: [sp for sp in trimmed_fasta[ortho]] for ortho in trimmed_fasta}, f)


shutil.rmtree(unpadded_primer_product_path)
os.makedirs(unpadded_primer_product_path, exist_ok=True)
for ortho in trimmed_fasta.keys():
    filename = unpadded_primer_product_path + ortho + ".13spp.fasta"
    with open(filename, "w") as f:
        for seqReq in sorted(trimmed_fasta[ortho].values(), key=lambda x: sp_order[x.id]):
            f.write(seqReq.format("fasta"))

shutil.rmtree(padded_primer_product_path)
os.makedirs(padded_primer_product_path, exist_ok=True)
for ortho in padded_fasta.keys():
    filename = padded_primer_product_path + ortho + ".13spp.fasta"
    with open(filename, "w") as f:
        for seqReq in sorted(padded_fasta[ortho].values(), key=lambda x: sp_order[x.id]):
            f.write(seqReq.format("fasta"))
