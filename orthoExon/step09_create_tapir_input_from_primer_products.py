#!/usr/bin/env python

import json
import os
import shutil

import config
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def pad_to_full_sp_count_and_trim_to_primer_product(alternate_sp_fn, padded_primer_product_path, unpadded_primer_product_path,
                                                    primer3_path, orthoCds_path, enhanced_species_list):
    with open(alternate_sp_fn, 'r') as f:
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
                primer[ortho] = []
                for variation in range(int(lines['PRIMER_PAIR_NUM_RETURNED'])):
                    left, l_len = lines['PRIMER_LEFT_{}'.format(variation)].split(',')
                    right, r_len = lines['PRIMER_RIGHT_{}'.format(variation)].split(',')
                    start = int(left) + int(l_len)
                    end = int(right) - int(r_len) + 1
                    primer[ortho].append((start, end))

    trimmed_fasta = {}
    padded_fasta = {}
    for ortho in fasta.keys():
        if ortho not in primer.keys():
            continue
        trimmed_fasta[ortho] = []
        padded_fasta[ortho] = []
        for variation in range(len(primer[ortho])):
            start, end = primer[ortho][variation]
            padding = {}
            for sp in enhanced_species_list:
                if sp not in fasta[ortho].keys():
                    for alt_sp in alternate_sp[sp]:
                        if alt_sp in fasta[ortho].keys():
                            seq = fasta[ortho][alt_sp].seq[start:end]
                            des = "PADDING"
                            padding[sp] = SeqRecord(seq, id=sp, description=des)
                            break
            trimmed_fasta[ortho].append({sp: fasta[ortho][sp][start:end] for sp in fasta[ortho]})
            padded_fasta[ortho].append(padding)
            padded_fasta[ortho][variation].update(trimmed_fasta[ortho][variation])

    for ortho in padded_fasta:
        for variation in range(len(padded_fasta[ortho])):
            for sp in padded_fasta[ortho][variation]:
                padded_fasta[ortho][variation][sp].description = ""

    # output fasta to pre_padding_species.json
    os.makedirs(json_path, exist_ok=True)
    filename = json_path + "pre_padding_species.json"
    with open(filename, 'w') as f:
        json.dump({ortho: [sp for sp in trimmed_fasta[ortho][0]] for ortho in trimmed_fasta}, f)

    shutil.rmtree(unpadded_primer_product_path, ignore_errors=True)
    os.makedirs(unpadded_primer_product_path, exist_ok=True)
    for ortho in trimmed_fasta.keys():
        for variation in range(len(trimmed_fasta[ortho])):
            filename = unpadded_primer_product_path + ortho + "_V" + str(variation) + ".13spp.fasta"
            with open(filename, "w") as f:
                for seqReq in sorted(trimmed_fasta[ortho][variation].values(), key=lambda x: config.sp_order[x.id]):
                    f.write(seqReq.format("fasta"))

    shutil.rmtree(padded_primer_product_path, ignore_errors=True)
    os.makedirs(padded_primer_product_path, exist_ok=True)
    for ortho in padded_fasta.keys():
        for variation in range(len(padded_fasta[ortho])):
            filename = padded_primer_product_path + ortho + "_V" + str(variation) + ".13spp.fasta"
            with open(filename, "w") as f:
                for seqReq in sorted(padded_fasta[ortho][variation].values(), key=lambda x: config.sp_order[x.id]):
                    f.write(seqReq.format("fasta"))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='This script creates primer3 input files')

    parser.add_argument('--orthoCds_path', help='orthoCds_path', default=config.orthoCds_path)
    parser.add_argument('--primer3_path', help='primer3_path', default=config.primer3_path)
    parser.add_argument('--padded_primer_product_path', help='padded_primer_product_path',
                        default=config.padded_primer_product_path)
    parser.add_argument('--unpadded_primer_product_path', help='unpadded_primer_product_path',
                        default=config.unpadded_primer_product_path)
    parser.add_argument('--alternate_sp_fn', help='alternate_sp.json path', default=config.alternate_sp_fn)

    args = parser.parse_args()

    pad_to_full_sp_count_and_trim_to_primer_product(args.alternate_sp_fn, args.padded_primer_product_path,
                                                    args.unpadded_primer_product_path, args.primer3_path,
                                                    args.orthoCds_path, config.enhanced_species_list)
