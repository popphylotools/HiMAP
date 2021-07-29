#!/usr/bin/env python

import collections
import itertools
import multiprocessing as mp
import os

import pandas as pd
from Bio import SeqIO
from Bio import motifs
from Bio.Data import IUPACData
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def process_sample(sample_dir, called_sequences_path):
    """for each individual (folder)
    for each target (file)
    combine reads (entries) and make a call
    """
    sample_nm = sample_dir.split("/")[-2]
    fastq_fn = {fn.split(".fastq")[0]: sample_dir + fn for fn in os.listdir(sample_dir) if
                (".fastq" in fn)}

    mots = dict()
    good_count = dict()
    trashed_count = dict()
    trashed_percent = dict()
    for sample_ortho in fastq_fn:
        reads = [seq_record.seq for seq_record in
                 SeqIO.parse(fastq_fn[sample_ortho], "fastq", alphabet=IUPAC.ambiguous_dna)]
        common_length = collections.Counter((len(seq) for seq in reads)).most_common()[0][0]

        good_reads = list(filter(lambda x: len(x) == common_length, reads))

        good_count[sample_ortho] = len(good_reads)
        trashed_count[sample_ortho] = len(reads) - len(good_reads)
        trashed_percent[sample_ortho] = trashed_count[sample_ortho] * 100. / (
            trashed_count[sample_ortho] + good_count[sample_ortho])

        mots[sample_ortho] = motifs.create(good_reads, alphabet=IUPACData.ambiguous_dna_letters)
        _seq = mots[sample_ortho].degenerate_consensus

        # write fasta for each sample_ortho
        os.makedirs(called_sequences_path + sample_nm, exist_ok=True)
        with open(called_sequences_path + sample_nm + '/' + sample_ortho + "_" + str(
                len(good_reads)) + "-reads.consensus.fasta", "w") as f:
            seqReq = SeqRecord(_seq, id=sample_ortho,
                               description="good_{} trashed_{}".format(good_count[sample_ortho],
                                                                       trashed_count[sample_ortho]))
            f.write(seqReq.format("fasta"))

    # write csv for each sample
    consensus = [(sample_ortho, mots[sample_ortho].degenerate_consensus) for sample_ortho in mots]
    df = pd.DataFrame(consensus, columns=["sample_ortho", "seq"])
    df["good_count"] = df["sample_ortho"].apply(lambda x: good_count[x]).astype(int)
    df["trashed_count"] = df["sample_ortho"].apply(lambda x: trashed_count[x]).astype(int)
    df["trashed_percent"] = df["sample_ortho"].apply(lambda x: trashed_percent[x]).astype(float)
    df["seq_len"] = df["seq"].apply(len)
    df["seq_str"] = df["seq"].apply(str)
    df[["sample_ortho", "good_count", "trashed_count", "trashed_percent", "seq_len", "seq_str"]].to_csv(
        called_sequences_path + sample_nm + ".csv", index=False)


def summerize(sample_dirs, called_sequences_path, filtered_summary_sequences_path, min_length, max_len_dev):
    frames = [pd.read_csv(called_sequences_path + sample_nm + ".csv", index_col=0, parse_dates=True) for sample_nm in sample_dirs.keys()]
    df = pd.concat(frames)
    df["sample_ortho"] = df.index
    df = df.loc[df['sample_ortho'].str.contains('.')]  # gets rid of "unknown" pools??
    df["sample"] = df["sample_ortho"].apply(lambda x: x.split('.',1)[0])
    df["ortho"] = df["sample_ortho"].apply(lambda x: x.split('.',1)[1])
    df.drop("sample_ortho", axis=1, inplace=True)

    # split("N * 40")[0]
    # if contains string of 40 N's, split there and take first half only
    df = df.loc[df["seq_str"].notnull()]
    df.loc[df["seq_str"].str.contains("N" * 40), "seq_str"] = df.loc[
        df["seq_str"].str.contains("N" * 40), "seq_str"].apply(
        lambda x: x.split("N" * 40)[0])

    # exclude if len < min_len
    df["seq_len"] = df["seq_str"].apply(len)
    df = df.loc[df["seq_len"] >= min_length]

    # exclude if len_dev > max_len_dev
    avg_len_df = df[["ortho", "seq_len"]].groupby("ortho").mean()
    df["len_deviation_pre_cut"] = df.apply(lambda x: abs(x["seq_len"] - avg_len_df.loc[x["ortho"]]), axis=1)
    df = df.loc[df["len_deviation_pre_cut"] <= max_len_dev]

    avg_len_df = df[["ortho", "seq_len"]].groupby("ortho").mean()
    df["len_deviation_post_cut"] = df.apply(lambda x: abs(x["seq_len"] - avg_len_df.loc[x["ortho"]]), axis=1)

    df.sort_values(by=['ortho', 'sample'], ascending=[True, True], inplace=True)
    df.reindex(["sample", "ortho", "good_count", "trashed_count", "trashed_percent", "len_deviation_pre_cut",
                          "len_deviation_post_cut", "seq_len", "seq_str"], axis=1)

    # write to csv
    os.makedirs(filtered_summary_sequences_path, exist_ok=True)
    df.to_csv(filtered_summary_sequences_path + "summary.csv", index=False)

    # create summary fasta file per ortho
    df["fasta_str"] = df.apply(
        lambda rec: SeqRecord(Seq(rec["seq_str"]), id=rec["sample"],
                              description="{} matched_{} trashed_{}".format(rec["ortho"], rec["good_count"],
                                                                            rec["trashed_count"])).format("fasta"),
        axis=1)

    os.makedirs(filtered_summary_sequences_path + "summary_fastas", exist_ok=True)
    for ortho in df["ortho"].unique():
        with open(filtered_summary_sequences_path + "summary_fastas/" + ortho + ".fasta", "w") as f:
            for fasta_str in df.loc[df["ortho"] == ortho, "fasta_str"]:
                f.write(fasta_str)


if __name__ == '__main__':
    import argparse
    import pytoml

    parser = argparse.ArgumentParser()
    parser.add_argument('--configPath', help='configPath', default='../config.toml')
    args = parser.parse_args()

    # load config file
    with open(args.configPath) as toml_data:
        config = pytoml.load(toml_data)

    sample_dirs = {val: config['sequences_path'] + val + "/" for val in
                   os.listdir(config['sequences_path']) if
                   os.path.isdir(config['sequences_path'] + val)}

    with mp.Pool(min(len(sample_dirs), mp.cpu_count())) as p:
        p.starmap(process_sample, zip(sample_dirs.values(), itertools.repeat(config['called_sequences_path'])))

    summerize(sample_dirs, config['called_sequences_path'], config['filtered_summary_sequences_path'],
              config['min_length'], config['max_len_dev'])
