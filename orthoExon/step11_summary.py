#!/usr/bin/env python

import json
import os
import sqlite3

import orthoExon.config as config
import gffutils
import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import IUPAC


def summary(json_path, db_path, primer3_path, unpadded_primer_product_path, pi_score_path, summary_fn):
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
                primer[ortho] = []
                for variation in range(int(lines['PRIMER_PAIR_NUM_RETURNED'])):
                    left, l_len = lines['PRIMER_LEFT_{}'.format(variation)].split(',')
                    right, r_len = lines['PRIMER_RIGHT_{}'.format(variation)].split(',')
                    primer[ortho].append((
                        str(int(lines['PRIMER_PAIR_{}_PRODUCT_SIZE'.format(variation)]) - int(l_len) - int(r_len)),
                        lines['PRIMER_LEFT_{}_SEQUENCE'.format(variation)],
                        lines['PRIMER_RIGHT_{}_SEQUENCE'.format(variation)],
                        lines['PRIMER_LEFT_{}_TM'.format(variation)],
                        lines['PRIMER_RIGHT_{}_TM'.format(variation)]))

    # grab pi scores from sql database
    name_score = []
    for fn in os.listdir(pi_score_path):
        conn = sqlite3.connect(pi_score_path + fn)
        name_score += conn.execute(
            "select loci.locus, avg(pi) from loci, discrete where loci.id = discrete.id group by loci.locus").fetchall()
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
    for ortho_plus in name_score:
        ortho = ortho_plus[:-3]
        variation = int(ortho_plus[-1:])
        for sp in pre_padd_sp[ortho]:
            if 'product' in gff[sp][parent_groups[ortho.split("_")[0]][sp]].attributes.keys():
                product = gff[sp][parent_groups[ortho.split("_")[0]][sp]]['product'][0]
            else:
                product = "N/A"
            score = name_score[ortho_plus]
            data.append((ortho_plus, str(score), sp, product, *(primer[ortho][variation])))

    data = sorted(data, key=lambda x: (x[0], config.sp_order[x[2]]))

    header = ['Exon_Name',
              'PI_Score',
              'Species',
              'Gene_Product',
              'Target_Sequence_Length',
              'PRIMER_LEFT_0_SEQUENCE',
              'PRIMER_RIGHT_0_SEQUENCE',
              'PRIMER_LEFT_0_TM',
              'PRIMER_RIGHT_0_TM']

    df = pd.DataFrame(data, columns=header)

    df["species"] = df.Exon_Name.apply(lambda ortho: [req.id for req in fasta[ortho].values()])

    filename = summary_fn
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    df.to_csv(summary_fn, index=False)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='This script creates primer3 input files')

    parser.add_argument('--summary_fn', help='summary_fn', default=config.summary_fn)
    parser.add_argument('--primer3_path', help='primer3_path', default=config.primer3_path)
    parser.add_argument('--unpadded_primer_product_path', help='unpadded_primer_product_path',
                        default=config.unpadded_primer_product_path)
    parser.add_argument('--db_path', help='db_path', default=config.db_path)
    parser.add_argument('--json_path', help='json_path', default=config.json_path)
    parser.add_argument('--pi_score_path', help='pi_score_path', default=config.pi_score_path)

    args = parser.parse_args()

    summary(args.json_path, args.db_path, args.primer3_path, args.unpadded_primer_product_path, args.pi_score_path,
            args.summary_fn)
