# coding: utf-8

# In[ ]:

import json
import os

import gffutils
from Bio import SeqIO
from Bio.Alphabet import IUPAC

# In[ ]:

full_species_list = ['Bjar', 'Aobl', 'Bmin', 'Asus', 'Btry', 'Afra', 'Blat', 'Bzon', 'Bcor',
                     'Ccap', 'Bcur', 'Bole', 'Bdor']
species_list = ["Bcur", "Bdor", "Bole", "Ccap"]
transvestigated_species_set = {'Bcor', 'Blat', 'Bzon', 'Afra', 'Bmin', 'Bjar', 'Aobl'}

output_path = "../output/"
p3_out_path = "../intermediate/primer_design/"
primer_products_path = "../output/primerProducts/"
db_path = "../intermediate/gff_databases/"
json_path = "../intermediate/json/"

# In[ ]:

# create handles for all .fasta files in fasta directory
fasta_fn = {name.split('.13spp.fasta')[0]: primer_products_path + name for name in
            os.listdir(primer_products_path) if
            ((".13spp.fasta" in name) and (".13spp.fasta.fai" not in name))}

# read and parse fasta files for each species
fasta = {}
for ortho in fasta_fn.keys():
    fasta[ortho] = {seq_record.id: seq_record
                    for seq_record in SeqIO.parse(fasta_fn[ortho],
                                                  "fasta", alphabet=IUPAC.ambiguous_dna)}
# In[ ]:

primer = {}
for p3_out_fn in os.listdir(p3_out_path):
    ortho = p3_out_fn.split('.degenerate.p3.out')[0]
    with open(p3_out_path + p3_out_fn, 'r') as f:
        lines = f.readlines()
        lines = [line.strip().split('=') for line in lines]
        lines = {key: value for key, value in lines if key is not ''}
        if lines['PRIMER_PAIR_NUM_RETURNED'] is not '0':
            left, l_len = lines['PRIMER_LEFT_0'].split(',')
            right, r_len = lines['PRIMER_RIGHT_0'].split(',')
            primer[ortho] = (
                str(int(lines['PRIMER_PAIR_0_PRODUCT_SIZE']) - int(l_len) - int(r_len)),
                lines['PRIMER_LEFT_0_SEQUENCE'],
                lines['PRIMER_RIGHT_0_SEQUENCE'],
                lines['PRIMER_LEFT_0_TM'],
                lines['PRIMER_RIGHT_0_TM'])

# In[ ]:

with open("input/net_PI_avg_edited.txt", 'r') as f:
    name_score = [line.strip().split() for line in f.readlines()[1:]]
name_score = {line[0].split(".13spp.fasta")[0]: line[1] for line in name_score}

# In[ ]:

# import pre_padding_species.json
with open(json_path + "pre_padding_species.json", 'r') as f:
    pre_padd_sp = json.load(f)

# In[ ]:

# import ortholog groups
with open(json_path + "groups.json", 'r') as f:
    parent_groups = json.load(f)

# In[ ]:

# create handles for all .db files in intermediate directory
gff_fn = {name.split('.gff.db')[0]: db_path + name for name in os.listdir(db_path) if
          ".gff.db" in name}
gff = {key: gffutils.FeatureDB(value) for key, value in gff_fn.items()}

# In[ ]:

data = []
for ortho in primer:
    for sp in pre_padd_sp[ortho]:
        if 'product' in gff[sp][parent_groups[ortho.split("_")[0]][sp]].attributes.keys():
            product = gff[sp][parent_groups[ortho.split("_")[0]][sp]]['product'][0]
        else:
            product = "N/A"
        score = name_score[ortho]
        data.append((ortho, score, sp, product, *primer[ortho]))

# In[ ]:

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

# In[ ]:

data = sorted(data, key=lambda x: (x[0], sp_order[x[2]]))

# In[ ]:

header = ['Exon_Name',
          'PI_Score',
          'Species',
          'Gene_Product',
          'Target_Sequence_Length',
          'PRIMER_LEFT_0_SEQUENCE',
          'PRIMER_RIGHT_0_SEQUENCE',
          'PRIMER_LEFT_0_TM',
          'PRIMER_RIGHT_0_TM']

# In[ ]:

with open(output_path + "13spp_exon_primer_data.csv", "w") as f:
    f.write(",".join(header))
    for record in data:
        f.write("\n" + ",".join(record))
