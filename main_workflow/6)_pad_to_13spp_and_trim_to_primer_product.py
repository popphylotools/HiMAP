# coding: utf-8

# In[ ]:

import json
import os

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

# In[ ]:

full_species_list = ['Bjar', 'Aobl', 'Bmin', 'Asus', 'Btry', 'Afra', 'Blat', 'Bzon', 'Bcor',
                     'Ccap', 'Bcur', 'Bole', 'Bdor']
species_list = ["Bcur", "Bdor", "Bole", "Ccap"]
transvestigated_species_set = {'Bcor', 'Blat', 'Bzon', 'Afra', 'Bmin', 'Bjar', 'Aobl'}

output_padded_path = "../intermediate/phylo_informativeness/"
output_unpadded_path = "../output/primerProducts/"
aligned_fasta_path = "../output/orthoCds/"
p3_out_path = "../intermediate/primer_design/"
json_path = "../intermediate/json/"

# In[ ]:

with open(json_path + "alternate_sp.json", 'r') as f:
    alternate_sp = json.load(f)

# In[ ]:

# create handles for all .fasta files in fasta directory
fasta_fn = {name.split('.13spp.fasta')[0]: aligned_fasta_path + name for name in
            os.listdir(aligned_fasta_path) if
            ((".13spp.fasta" in name) and (".13spp.fasta.fai" not in name))}

# In[ ]:

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
            test = lines
            left, l_len = lines['PRIMER_LEFT_0'].split(',')
            right, r_len = lines['PRIMER_RIGHT_0'].split(',')
            start = int(left) + int(l_len)
            end = int(right) - int(r_len) + 1
            primer[ortho] = (start, end)

# In[ ]:

# trimmed_fasta = {ortho: {sp: fasta[ortho][sp][start:end] for sp in fasta[ortho]}}

# In[ ]:

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
                    des = fasta[ortho][alt_sp].description
                    des = "PADDING " + des
                    padding[sp] = SeqRecord(seq, id=sp, description=des)
                    break
    trimmed_fasta[ortho] = {sp: fasta[ortho][sp][start:end] for sp in fasta[ortho]}
    padded_fasta[ortho] = padding
    padded_fasta[ortho].update(trimmed_fasta[ortho])

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

# output fasta to pre_padding_species.json
with open(json_path + "pre_padding_species.json", 'w') as f:
    json.dump({ortho: [sp for sp in trimmed_fasta[ortho]] for ortho in trimmed_fasta}, f)

# In[ ]:

for ortho in trimmed_fasta.keys():
    with open(output_unpadded_path + ortho + ".13spp.fasta", "w") as f:
        for seqReq in sorted(trimmed_fasta[ortho].values(), key=lambda x: sp_order[x.id]):
            f.write(seqReq.format("fasta"))

# In[ ]:

for ortho in padded_fasta.keys():
    with open(output_padded_path + ortho + ".13spp.fasta", "w") as f:
        for seqReq in sorted(padded_fasta[ortho].values(), key=lambda x: sp_order[x.id]):
            f.write(seqReq.format("fasta"))
