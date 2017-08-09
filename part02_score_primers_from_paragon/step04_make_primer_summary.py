#!/usr/bin/env python

import json
import multiprocessing as mp
import os
import shutil
import sqlite3
import subprocess
from collections import Counter
from functools import reduce
from multiprocessing.pool import ThreadPool
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


def consensus(aligned_seq_list):
    con = ""
    length = len(aligned_seq_list[0])
    for seq in aligned_seq_list:
        assert len(seq) == length
    for loc in range(length):
        base_set = set.union(*[expand_iupac[seq[loc].upper()] for seq in aligned_seq_list])
        if '-' in base_set:
            con += '-'
        else:
            con += collapse_iupac[tuple(sorted(base_set))]
    return con


def create_padded_primer_products(orthoExon_path, padded_primer_product_path, unpadded_primer_product_path,
                                  alternate_sp_fn,
                                  paragon_fn):
    # create handles for all .fasta files in fasta directory
    fasta_fn = {name.split('.full.fasta')[0]: orthoExon_path + name for name in
                os.listdir(orthoExon_path) if
                ((".full.fasta" in name) and (".full.fasta.fai" not in name))}

    # read and parse fasta files for each species
    fasta = {}
    for ortho in fasta_fn.keys():
        fasta[ortho] = {seq_record.id: seq_record
                        for seq_record in SeqIO.parse(fasta_fn[ortho],
                                                      "fasta", alphabet=IUPAC.ambiguous_dna)}

    # load alternate_sp
    with open(alternate_sp_fn, 'r') as f:
        alternate_sp = json.load(f)

    # load paragon data
    with open(paragon_fn) as f:
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
    df["consensus"] = df.loci_ID.apply(lambda ortho: consensus([req.seq for req in fasta[ortho].values()]))
    df["left_primer"] = df.apply(lambda req: req["consensus"][req["amp_start"]:req["ampInsert_start"]], axis=1)
    df["insert"] = df.apply(lambda req: req["consensus"][req["ampInsert_start"]:req["ampInsert_end"]], axis=1)
    df["right_primer"] = df.apply(lambda req: req["consensus"][req["ampInsert_end"]:req["amp_end"]], axis=1)

    # set column names lowercase
    df.columns = map(str.lower, df.columns)

    # output unpadded primer products
    shutil.rmtree(unpadded_primer_product_path, ignore_errors=True)
    os.makedirs(unpadded_primer_product_path, exist_ok=True)
    for ortho in df.loci_id:
        with open(unpadded_primer_product_path + ortho + ".fasta", 'w') as f:
            start = df.loc[ortho]["ampinsert_start"]
            end = df.loc[ortho]["ampinsert_end"]
            for sp in full_species_list:
                if sp in df.loc[ortho]["species"]:
                    seqReq = fasta[ortho][sp][start:end]
                    f.write(seqReq.format("fasta"))

    # output padded primer products
    shutil.rmtree(padded_primer_product_path, ignore_errors=True)
    os.makedirs(padded_primer_product_path, exist_ok=True)
    for ortho in df.loci_id:
        with open(padded_primer_product_path + ortho + ".padded.fasta", 'w') as f:
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

    return df, fasta


def convertfasta2nex(padded_primer_product_fn, nex_fn):
    p = subprocess.Popen(["./convertfasta2nex_driver.sh", padded_primer_product_fn, nex_fn],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return out, err


def tapir_driver(nex_sub_path, tapir_out_sub_path, ref_tree_fn):
    p = subprocess.Popen(["./tapir_driver.sh", nex_sub_path, tapir_out_sub_path, ref_tree_fn],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    with open(tapir_out_sub_path + "tapir.out", 'wb') as f:
        f.write(out)
    with open(tapir_out_sub_path + "tapir.err", 'wb') as f:
        f.write(err)
    return out, err


def phylogenetic_informativeness(padded_primer_product_path, nex_path, tapir_out_path, pi_score_path, ref_tree_fn):
    cpu_count = mp.cpu_count()

    # remove and recreate nex and tapir ouput directories
    shutil.rmtree(nex_path, ignore_errors=True)
    shutil.rmtree(tapir_out_path, ignore_errors=True)
    shutil.rmtree(pi_score_path, ignore_errors=True)
    os.makedirs(pi_score_path, exist_ok=True)

    # create tapir input and output directories for each core
    for i in range(cpu_count):
        os.makedirs(nex_path + str(i), exist_ok=True)
        os.makedirs(tapir_out_path + str(i), exist_ok=True)

    # for each fasta file, make a nex (split into cpu_count subdirectories)
    i = 0
    convertfasta2nex_input = []
    for file in os.listdir(padded_primer_product_path):
        group = str(i % cpu_count)
        padded_primer_product_fn = padded_primer_product_path + file
        nex_fn = nex_path + group + "/" + file + ".nex"
        convertfasta2nex_input.append((padded_primer_product_fn, nex_fn))
        i += 1

    pool = ThreadPool(cpu_count)
    pool.starmap(convertfasta2nex, convertfasta2nex_input)
    pool.close()
    pool.join()

    # for each core, process a subdirectory
    tapir_driver_input = []
    for i in range(cpu_count):
        nex_sub_path = nex_path + str(i) + "/"
        tapir_out_sub_path = tapir_out_path + str(i) + "/"
        tapir_driver_input.append((nex_sub_path, tapir_out_sub_path, ref_tree_fn))

    pool = ThreadPool(cpu_count)
    pool.starmap(tapir_driver, tapir_driver_input)
    pool.close()
    pool.join()

    # collect sql to a directory
    sql_list = []
    for i in range(cpu_count):
        tapir_out_sub_path = tapir_out_path + str(i) + "/"
        sql_list.extend([(tapir_out_sub_path + fn, pi_score_path + str(i) + '_' + fn )for fn in os.listdir(tapir_out_sub_path) if ".sqlite" in fn])
    for old_fn, new_fn in sql_list:
        shutil.move(old_fn, new_fn)


def summarize(df, fasta, pi_score_path, summary_fn):
    # grab pi scores from sql database
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
    df["max_ambiguities(l|r)"] = df.apply(
        lambda row: max(amb_count(row['left_primer']), amb_count(row['right_primer'])),
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
    df.to_csv(summary_fn, index=False,
              columns=["loci_id", "amplicon_id", "score",
                       'total_ambiguities(l+r)', 'max_ambiguities(l|r)', 'left_ambiguities', 'right_ambiguities',
                       "total_primer_versions", "l_primer_versions", "r_primer_versions",
                       'total_combinations(l+r)', 'left_combinations', 'right_combinations',
                       "insert_length", "ampinsert_start", "ampinsert_end",
                       "amp_len", "amp_start", "amp_end",
                       "left_primer", "right_primer",
                       "comment", "species"])


if __name__ == '__main__':
    import argparse
    import pytoml

    parser = argparse.ArgumentParser()
    parser.add_argument('--configPath', help='configPath', default='../config.toml')
    args = parser.parse_args()

    # load config file
    with open(args.configPath) as toml_data:
        config = pytoml.load(toml_data)

    working_df, working_fasta = create_padded_primer_products(config['orthoExon_path'],
                                                              config['padded_primer_product_path'],
                                                              config['unpadded_primer_product_path'],
                                                              config['alternate_sp_fn'], config['paragon_fn'])

    phylogenetic_informativeness(config['padded_primer_product_path'], config['nex_path'], config['tapir_out_path'],
                                 config['pi_score_path'], config['ref_tree_fn'])

    summarize(working_df, working_fasta, config['pi_score_path'], config['summary_fn'])
