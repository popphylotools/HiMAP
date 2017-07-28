import json
import os
import sqlite3
import subprocess
from collections import Counter
from functools import reduce
from operator import mul

import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import IUPAC

# utility data structures and functions
full_species_list = ['Bcur', 'Bdor', 'Bole', 'Ccap', 'Bcor', 'Blat', 'Bzon', 'Afra', 'Bmin', 'Bjar', 'Aobl', 'Asus',
                     'Btry']

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

possibilities = {key: len(value) for key, value in expand_iupac.items()}


def amb_count(seq):
    return sum([value for key, value in Counter(seq).items() if key in "yrwskmdvhbn-"])


def amb_combinations(seq):
    return reduce(mul, [possibilities[nuc.upper()] for nuc in seq])


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


def amb_count(seq):
    return sum([value for key, value in Counter(seq).items() if key in "yrwskmdvhbn-"])

# load fasta data

# fasta path
unpadded_primer_product_path = "./orthoCds/"

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

# load alternate_sp
with open("alternate_sp.json", 'r') as f:
    alternate_sp = json.load(f)

# load paragon data
with open("./PGC106.selectedAmp.tab") as f:
    lines = f.readlines()
    df = pd.DataFrame([line.strip().split() for line in lines[1:]], columns=lines[0].strip().split())

# parse ints
df.amp_start = df.amp_start.astype(int)
df.amp_end = df.amp_end.astype(int)
df.ampInsert_start = df.ampInsert_start.astype(int)
df.ampInsert_end = df.ampInsert_end.astype(int)
df.amp_len = df.amp_len.astype(int)
df.index = df.loci_ID

# add columns
df["species"] = df.loci_ID.apply(lambda ortho: [req.id for req in fasta[ortho].values()])
df["consensus"] = df.loci_ID.apply(lambda ortho: Consensus([req.seq for req in fasta[ortho].values()]))
df["left_primer"] = df.apply(lambda req: req["consensus"][req["amp_start"]:req["ampInsert_start"]], axis=1)
df["insert"] = df.apply(lambda req: req["consensus"][req["ampInsert_start"]:req["ampInsert_end"]], axis=1)
df["right_primer"] = df.apply(lambda req: req["consensus"][req["ampInsert_end"]:req["amp_end"]], axis=1)

# set column names lowercase
df.columns = map(str.lower, df.columns)

# output unpadded primer products
for ortho in df.loci_id:
    with open("primerProducts/" + ortho + ".fasta", 'w') as f:
        start = df.loc[ortho]["ampinsert_start"]
        end = df.loc[ortho]["ampinsert_end"]
        for sp in full_species_list:
            if sp in df.loc[ortho]["species"]:
                seqReq = fasta[ortho][sp][start:end]
                f.write(seqReq.format("fasta"))

# output padded primer products
for ortho in df.loci_id:
    with open("P3/fasta/" + ortho + ".padded.fasta", 'w') as f:
        start = df.loc[ortho]["ampinsert_start"]
        end = df.loc[ortho]["ampinsert_end"]
        for sp in full_species_list:
            if sp in df.loc[ortho]["species"]:
                seqReq = fasta[ortho][sp][start:end]
            else:
                for alt_sp in alternate_sp[sp]:
                    if alt_sp in fasta[ortho].keys():
                        seqReq = fasta[ortho][alt_sp][start:end]
                        break
            seqReq.description = ""
            seqReq.id = sp
            seqReq.name = sp
            f.write(seqReq.format("fasta"))

# calculate phylogenetic informativeness
subprocess.call(["./pf.sh", "./P3"])

# grab pi scores from sql database
pi_score_path = "./P3/sql/"
name_score = []
for fn in os.listdir(pi_score_path):
    conn = sqlite3.connect(pi_score_path + fn)
    name_score += conn.execute(
        "select loci.locus, avg(pi) from loci, discrete where loci.id = discrete.id group by loci.locus").fetchall()
name_score_df = pd.DataFrame([(line[0].split(".padded.fasta")[0], line[1]) for line in name_score],
                             columns=["ortho", "score"])
name_score_df.index = name_score_df.ortho
name_score_df = name_score_df.drop("ortho", axis=1)

name_score_df.score = name_score_df.score.astype(float)
df = pd.concat([df, name_score_df], axis=1)

# create more columns
df["insert_length"] = df["insert"].apply(len)
df["total_ambiguities(l+r)"] = df.apply(lambda row: amb_count(row['left_primer']) + amb_count(row['right_primer']),
                                        axis=1)
df["left_ambiguities"] = df.apply(lambda row: amb_count(row['left_primer']), axis=1)
df["right_ambiguities"] = df.apply(lambda row: amb_count(row['right_primer']), axis=1)
df["max_ambiguities(l|r)"] = df.apply(lambda row: max(amb_count(row['left_primer']), amb_count(row['right_primer'])),
                                      axis=1)
df["total_combinations(l+r)"] = df.apply(
    lambda row: amb_combinations(row['left_primer']) + amb_combinations(row['right_primer']), axis=1)
df["left_combinations"] = df.apply(lambda row: amb_combinations(row['left_primer']), axis=1)
df["right_combinations"] = df.apply(lambda row: amb_combinations(row['right_primer']), axis=1)

# add count of unique sequences found in data which match primers
primer_set = []
for ortho in df.loci_id:

    l_primer_set = set()
    r_primer_set = set()

    l_start = df.loc[ortho]["amp_start"]
    l_end = df.loc[ortho]["ampinsert_start"]
    r_start = df.loc[ortho]["ampinsert_end"]
    r_end = df.loc[ortho]["amp_end"]

    for sp in df.loc[ortho]["species"]:
        l_primer_set.add(str(fasta[ortho][sp][l_start:l_end].seq))
        r_primer_set.add(str(fasta[ortho][sp][r_start:r_end].seq))

    primer_set.append((ortho, len(l_primer_set), len(r_primer_set), len(l_primer_set) + len(r_primer_set)))

primer_set_df = pd.DataFrame(primer_set,
                             columns=["ortho", "l_primer_versions", "r_primer_versions", "total_primer_versions"])
primer_set_df.index = primer_set_df.ortho
del primer_set_df["ortho"]

df = pd.concat([df, primer_set_df], axis=1)

# sort
df.sort_values(by='loci_id', ascending=True, inplace=True)

# output
df.to_csv("summary.csv", index=False,
          columns=["loci_id", "amplicon_id", "score",
                   'total_ambiguities(l+r)', 'max_ambiguities(l|r)', 'left_ambiguities', 'right_ambiguities',
                   "total_primer_versions", "l_primer_versions", "r_primer_versions",
                   'total_combinations(l+r)', 'left_combinations', 'right_combinations',
                   "insert_length", "ampinsert_start", "ampinsert_end",
                   "amp_len", "amp_start", "amp_end",
                   "left_primer", "right_primer",
                   "comment", "species"])
