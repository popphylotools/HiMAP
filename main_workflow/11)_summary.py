#!/usr/bin/env python
# coding: utf-8

import os
import gffutils
from Bio.Alphabet import IUPAC
import json
from Bio import SeqIO
import sqlite3
import argparse
from config import full_species_list, species_list, transvestigated_species_set, summary_fn, \
    primer3_path, unpadded_primer_product_path, db_path, json_path, pi_score_path

parser = argparse.ArgumentParser(description='This script creates primer3 input files')
parser.add_argument('--summary_fn', help='summary_fn', default=summary_fn)
parser.add_argument('--primer3_path',help='primer3_path', default=primer3_path)
parser.add_argument('--unpadded_primer_product_path', help='unpadded_primer_product_path', default=unpadded_primer_product_path)
parser.add_argument('--db_path',help='db_path', default=db_path)
parser.add_argument('--json_path', help='json_path', default=json_path)
parser.add_argument('--pi_score_path',help='pi_score_path', default=pi_score_path)

args = parser.parse_args()

summary_fn = args.summary_fn
primer3_path = args.primer3_path
unpadded_primer_product_path = args.unpadded_primer_product_path
db_path = args.db_path
json_path = args.json_path
pi_score_path = args.pi_score_path

# create handles for all .fasta files in fasta directory
fasta_fn = {name.split('.13spp.fasta')[0]: unpadded_primer_product_path + name for name in
            os.listdir(unpadded_primer_product_path) if
            ((".13spp.fasta" in name) and (".13spp.fasta.fai" not in name))}

# read and parse fasta files for each species
fasta = {}
for ortho in fasta_fn.keys():
    fasta[ortho] = {seq_record.id: seq_record
                    for seq_record in SeqIO.parse(fasta_fn[ortho],
                                                  "fasta", alphabet=IUPAC.ambiguous_dna)}

primer = {}
for p3_out_fn in [fn for fn in os.listdir(primer3_path) if ".p3.out" in fn]:
    ortho = p3_out_fn.split('.degenerate.p3.out')[0]
    with open(primer3_path + p3_out_fn, 'r') as f:
        lines = f.readlines()
        lines = [line.strip().split('=') for line in lines]
        lines = {key: value for key, value in lines if key is not ''}
        if 'PRIMER_PAIR_NUM_RETURNED' not in lines.keys():
            continue
        if lines['PRIMER_PAIR_NUM_RETURNED'] is not '0':
            left, l_len = lines['PRIMER_LEFT_0'].split(',')
            right, r_len = lines['PRIMER_RIGHT_0'].split(',')
            primer[ortho] = (
                str(int(lines['PRIMER_PAIR_0_PRODUCT_SIZE']) - int(l_len) - int(r_len)),
                lines['PRIMER_LEFT_0_SEQUENCE'],
                lines['PRIMER_RIGHT_0_SEQUENCE'],
                lines['PRIMER_LEFT_0_TM'],
                lines['PRIMER_RIGHT_0_TM'])

conn = sqlite3.connect(pi_score_path)
name_score = conn.execute("select loci.locus, pi from loci, net where loci.id = net.id").fetchall()
name_score = {line[0].split(".13spp.fasta")[0]: line[1] for line in name_score}

# import pre_padding_species.json
with open(json_path + "pre_padding_species.json", 'r') as f:
    pre_padd_sp = json.load(f)

# import ortholog groups
with open(json_path + "groups.json", 'r') as f:
    parent_groups = json.load(f)

# create handles for all .db files in intermediate directory
gff_fn = {name.split('.gff.db')[0]: db_path + name for name in os.listdir(db_path) if
          ".gff.db" in name}
gff = {key: gffutils.FeatureDB(value) for key, value in gff_fn.items()}

data = []
for ortho in primer:
    for sp in pre_padd_sp[ortho]:
        if 'product' in gff[sp][parent_groups[ortho.split("_")[0]][sp]].attributes.keys():
            product = gff[sp][parent_groups[ortho.split("_")[0]][sp]]['product'][0]
        else:
            product = "N/A"
        score = name_score[ortho]
        data.append((ortho, str(score), sp, product, *primer[ortho]))

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

data = sorted(data, key=lambda x: (x[0], sp_order[x[2]]))

header = ['Exon_Name',
          'PI_Score',
          'Species',
          'Gene_Product',
          'Target_Sequence_Length',
          'PRIMER_LEFT_0_SEQUENCE',
          'PRIMER_RIGHT_0_SEQUENCE',
          'PRIMER_LEFT_0_TM',
          'PRIMER_RIGHT_0_TM']

import csv
filename = summary_fn
os.makedirs(os.path.dirname(filename), exist_ok=True)
with open(filename, "w") as f:
    f.write(",".join(header))
    for record in data:
        f.write("\n" + ",".join(record))
