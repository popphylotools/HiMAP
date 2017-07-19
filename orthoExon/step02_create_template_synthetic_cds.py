#!/usr/bin/env python

import json
import os
import shutil

import gffutils
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta

from config import n_count


def create_template_synthetic_cds(template_species_list, fasta_path, template_species_alignment_path, db_path,
                                  json_path):
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

    # concatinate cds's for each species,ortho and output a fasta for each ortho
    nnn = Seq("".join([n for n in range(n_count)]), IUPAC.ambiguous_dna)
    shutil.rmtree(template_species_alignment_path, ignore_errors=True)
    os.makedirs(template_species_alignment_path, exist_ok=True)
    for ortho in parent_groups:
        filename = template_species_alignment_path + ortho + ".fasta"
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


if __name__ == "__main__":
    from config import template_species_list, fasta_path, template_species_alignment_path, db_path, json_path

    create_template_synthetic_cds(template_species_list, fasta_path, template_species_alignment_path, db_path,
                                  json_path)
