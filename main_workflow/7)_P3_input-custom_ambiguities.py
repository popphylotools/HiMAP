#!/usr/bin/env python
# coding: utf-8

# In[2]:

import os
import gffutils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from pyfaidx import Fasta
from Bio import SeqIO
import re
from pprint import pprint as pp


# In[6]:

full_species_list = ['Bjar', 'Aobl', 'Bmin', 'Asus', 'Btry', 'Afra', 'Blat', 'Bzon', 'Bcor', 'Ccap', 'Bcur', 'Bole', 'Bdor']
species_list = ["Bcur", "Bdor", "Bole", "Ccap"]
transvestigated_species_set = {'Bcor', 'Blat', 'Bzon', 'Afra', 'Bmin', 'Bjar', 'Aobl'}

output_path = "../intermediate/primer_design/"
aligned_fasta_path = "../output/orthoCds/"


# In[4]:

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

expand_iupac = {value.upper():set(key) for key,value in collapse_iupac.items()}

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
    return(consensus)


# In[8]:

# create handles for all .fasta files in fasta directory
fasta_fn = {name.split('.13spp.fasta')[0]: aligned_fasta_path + name for name in os.listdir(aligned_fasta_path) if
         ((".13spp.fasta" in name) and (".13spp.fasta.fai" not in name))}


# In[9]:

# read and parse fasta files for each species
fasta = {}
for ortho in fasta_fn.keys():
    fasta[ortho] = {seq_record.id : seq_record 
                                      for seq_record in SeqIO.parse(fasta_fn[ortho],
                                                                    "fasta", alphabet=IUPAC.ambiguous_dna)}


# In[10]:

from Bio import motifs
fasta_degenerate = {}
for ortho in fasta:
    seq = Consensus([fasta[ortho][sp].upper().seq for sp in fasta[ortho].keys()])
    fasta_degenerate[ortho] = seq


# In[12]:

# output
primer_product_size_range = '200-10000'
primer_thermodynamic_parameters_path = '/data0/opt/Primer3/primer3-2.3.6/src/primer3_config/'
primer_max_ns_accepted = '1'
primer_liberal_base = '1'
for ortho in fasta.keys():
    with open(output_path + ortho + ".degenerate.p3", "w") as f:
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

