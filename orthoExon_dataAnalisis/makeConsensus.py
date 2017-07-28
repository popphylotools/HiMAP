#!/usr/bin/env python

import collections
import multiprocessing as mp
import os

import Bio
import pandas as pd


def process_sample(sample_dir):
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
                 Bio.SeqIO.parse(fastq_fn[sample_ortho], "fastq", alphabet=Bio.IUPAC.ambiguous_dna)]
        common_length = collections.Counter([len(seq) for seq in reads]).most_common()[0][0]

        good_reads = list(filter(lambda x: len(x) == common_length, reads))

        good_count[sample_ortho] = len(good_reads)
        trashed_count[sample_ortho] = len(reads) - len(good_reads)
        trashed_percent[sample_ortho] = trashed_count[sample_ortho] * 100. / (
            trashed_count[sample_ortho] + good_count[sample_ortho])

        mots[sample_ortho] = Bio.motifs.create(good_reads)
        _seq = mots[sample_ortho].degenerate_consensus

        # write fasta for each sample_ortho
        with open(sample_dir + sample_ortho + "_" + str(len(good_reads)) + "-reads.consensus.fasta", "w") as f:
            seqReq = Bio.SeqRecord(_seq, id=sample_ortho,
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
        sample_dir + sample_nm + ".csv", index=False)


def summerize(sample_dirs):
    frames = [pd.DataFrame.from_csv(s_dir + sample_nm + ".csv") for sample_nm, s_dir in sample_dirs.items()]
    df = pd.concat(frames)
    df["sample_ortho"] = df.index
    df = df.loc[df['sample_ortho'].str.contains('orth')]  # gets rid of "unknown" pools
    df["sample"] = df["sample_ortho"].apply(lambda x: x.split('orth')[0][:-1])
    df["ortho"] = df["sample_ortho"].apply(lambda x: 'orth' + x.split('orth')[1])
    df.drop("sample_ortho", axis=1, inplace=True)

    # split("N * 40")[0]
    df = df.loc[df["seq_str"].notnull()]
    df.loc[df["seq_str"].str.contains("N" * 40), "seq_str"] = df.loc[
        df["seq_str"].str.contains("N" * 40), "seq_str"].apply(
        lambda x: x.split("N" * 40)[0])

    # cut len < 65
    df["seq_len"] = df["seq_str"].apply(len)
    df = df.loc[df["seq_len"] >= 65]

    # cut len_dev > 20
    avg_len_df = df[["ortho", "seq_len"]].groupby("ortho").mean()
    df["len_deviation_pre_cut"] = df.apply(lambda x: abs(x["seq_len"] - avg_len_df.ix[x["ortho"]]), axis=1)
    df = df.loc[df["len_deviation_pre_cut"] <= 20]

    avg_len_df = df[["ortho", "seq_len"]].groupby("ortho").mean()
    df["len_deviation_post_cut"] = df.apply(lambda x: abs(x["seq_len"] - avg_len_df.ix[x["ortho"]]), axis=1)

    df.sort_values(by=['ortho', 'sample'], ascending=[True, True], inplace=True)
    df = df.reindex_axis(["sample", "ortho", "good_count", "trashed_count", "trashed_percent", "len_deviation_pre_cut",
                          "len_deviation_post_cut", "seq_len", "seq_str"], axis=1)

    # write to csv
    df.to_csv("summary.csv", index=False)

    # create summary fasta file per ortho
    df["fasta_str"] = df.apply(
        lambda rec: Bio.SeqRecord(Bio.Seq(rec["seq_str"], alphabet=Bio.IUPAC.ambiguous_dna), id=rec["sample"],
                                  description="{} matched_{} trashed_{}".format(rec["ortho"], rec["good_count"],
                                                                                rec["trashed_count"])).format("fasta"),
        axis=1)

    os.makedirs("summary_fasta", exist_ok=True)
    for ortho in df["ortho"].unique():
        with open("summary_fasta/" + ortho + ".fasta", "w") as f:
            for fasta_str in df.loc[df["ortho"] == ortho, "fasta_str"]:
                f.write(fasta_str)


if __name__ == '__main__':
    sample_dirs = {val: "./data/" + val + "/" for val in os.listdir("./data") if os.path.isdir("./data/" + val)}

    with mp.Pool(min(len(sample_dirs), mp.cpu_count())) as p:
        p.map(process_sample, sample_dirs.values())

    summerize(sample_dirs)
