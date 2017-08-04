#!/usr/bin/env python

import json
import multiprocessing as mp
import os
import shutil
from multiprocessing.pool import ThreadPool

import gffutils
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta


def create_padded_cds(template_species_list, fasta_path, template_alignment_path, db_path, json_path, n_count):
    # create handles for all .db files in intermediate directory
    gff = {name.split('.gff.db')[0]: name for name in os.listdir(db_path) if ".gff.db" in name}
    gff = {key: gffutils.FeatureDB(db_path + value) for key, value in gff.items()}

    # create handles for all .fasta files in fasta directory
    fasta = {name.split('.nt.fasta')[0]: name for name in os.listdir(fasta_path) if
             ((".nt.fasta" in name) and (".nt.fasta.fai" not in name))}
    fasta = {key: Fasta(fasta_path + value) for key, value in fasta.items()}

    # import ortholog groups
    with open(json_path + "groups.json", 'r') as f:
        parent_groups = json.load(f)

    # concatenate cds's for each species,ortho and output a fasta for each ortho
    nnn = Seq('n' * n_count, IUPAC.ambiguous_dna)
    shutil.rmtree(template_alignment_path, ignore_errors=True)
    os.makedirs(template_alignment_path, exist_ok=True)
    for ortho in parent_groups:
        filename = template_alignment_path + ortho + "template.fasta"
        with open(filename, "w") as f:
            for sp in template_species_list:
                parent = gff[sp][parent_groups[ortho][sp]]
                strand = parent.strand
                cds_list = gff[sp].children(parent, featuretype="CDS", order_by="start")
                cat_seq = Seq("", IUPAC.ambiguous_dna)
                for i, cds in enumerate(cds_list):
                    if i > 0:
                        cat_seq += nnn
                    cat_seq += Seq(str(cds.sequence(fasta=fasta[sp], use_strand=False)),
                                   IUPAC.ambiguous_dna)
                if strand == '-':
                    cat_seq = cat_seq.reverse_complement()
                seqReq = SeqRecord(cat_seq, id=sp, description=parent.id)
                f.write(seqReq.format("fasta"))


def mafft_driver_file(file):
    p = subprocess.Popen(["./mafft_driver.sh", file, file + ".aln"],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return out, err


def mafft_driver_path(path):
    # remove old alignments
    rm_files = [path + file for file in os.listdir(path) if ".aln" in file]
    for file in rm_files:
        os.remove(file)

    # call maft on each fasta
    files = [path + file for file in os.listdir(path) if ".fasta" in file]
    pool = ThreadPool(mp.cpu_count())
    pool.map(mafft_driver_file, files)
    pool.close()
    pool.join()

if __name__ == "__main__":
    import subprocess
    import argparse
    import pytoml

    parser = argparse.ArgumentParser()
    parser.add_argument('--configPath', help='configPath', default='../config.toml')
    args = parser.parse_args()

    # load config file
    with open(args.configPath) as toml_data:
        config = pytoml.load(toml_data)

    create_padded_cds(config['template_species_list'], config['fasta_path'], config['template_alignment_path'],
                      config['db_path'], config['json_path'], config['n_count'])

    mafft_driver_path(config['template_alignment_path'])
