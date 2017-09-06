#!/usr/bin/env python

import gffutils
import json
import logging
import multiprocessing as mp
import os
import re
import shutil
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing.pool import ThreadPool
from pyfaidx import Fasta


# define functions to parse coordinates of cds's from concatenated aligned fasta w/ n's and -'s

def findBreakpoints(seq, n_count):
    breakpoints = []
    loc = 0
    regex = re.compile(r"n+[-+n]*")
    while True:
        # loc = seq.find(nnn, loc)
        match = regex.search(seq, loc)
        if not match:
            break
        if len(match.group().replace('-', '')) >= n_count:
            breakpoints.append(match.span())
        loc = match.end()
    return breakpoints


def findExonCoords(seq, n_count):
    breakpoints = findBreakpoints(seq, n_count)
    length = len(seq)

    if len(breakpoints) == 0:
        return [(0, length)]

    if len(breakpoints) == 1:
        bp = breakpoints[0]
        return [(0, bp[0]), (bp[1], length)]

    elif len(breakpoints) > 0:
        exonCoords = [(0, breakpoints[0][0])]

        for i in range(len(breakpoints) + 1)[1:-1]:  # all intermediate exons
            ex_start = breakpoints[i - 1][1]
            ex_end = breakpoints[i][0]
            exonCoords.append((ex_start, ex_end))

        exonCoords.append((breakpoints[-1][1], length))  # last exon
        return exonCoords


def gapPercent(seq):
    seq = str(seq)
    gappedLen = len(seq)
    gapCount = seq.count('-')
    return (100.0 * gapCount) / gappedLen


def longestGap(seq):
    seq = str(seq)
    gap_regex = re.compile(r"-+")
    gap_list = gap_regex.findall(seq)
    if gap_list:
        return sorted([len(gap) for gap in gap_list], reverse=True)[0]
    else:
        return 0


def create_raw_exons(fasta_path, enhanced_alignment_path, template_alignment_path, db_path, json_path,
                     template_species_list, max_gap_percent, max_gap_length,
                     min_exon_length, n_count):
    # create handles for all .db files in intermediate directory
    gff_fn = {name.split('.gff.db')[0]: db_path + name for name in os.listdir(db_path) if
              ".gff.db" in name}
    gff = {key: gffutils.FeatureDB(value) for key, value in gff_fn.items()}

    # create handles for all .fasta files in fasta directory
    fasta_fn = {name.split('.')[0]: fasta_path + name for name in os.listdir(fasta_path) if
                ((".fasta" in name) and (".fai" not in name))}
    fasta = {key: Fasta(value) for key, value in fasta_fn.items()}

    # import ortholog groups
    with open(json_path + "groups.json", 'r') as f:
        parent_groups = json.load(f)

    # create handles for all .fasta files in aligned template fasta directory
    aligned_fasta_fn = {name.split('.template.fasta')[0]: template_alignment_path + name for name in
                        os.listdir(template_alignment_path) if
                        ((".fasta.aln" in name) and (".fasta.aln.fai" not in name))}

    # read and parse fasta files for each species
    aligned_fasta = {}
    for ortho in aligned_fasta_fn.keys():
        aligned_fasta[ortho] = {seq_record.id: seq_record
                                for seq_record in SeqIO.parse(aligned_fasta_fn[ortho],
                                                              "fasta", alphabet=IUPAC.ambiguous_dna)}

    # parse coords from aligned fasta's
    coords = {}  # coords[ortho][sp] = [coord, ]
    for ortho in aligned_fasta:
        coords[ortho] = {}
        for sp in template_species_list:
            coords[ortho][sp] = findExonCoords(str(aligned_fasta[ortho][sp].seq), n_count)

    # Filter aligned exons
    ortho_coords = {}
    coord_index_dict = {}
    for ortho in coords:
        ortho_coords[ortho] = {}
        coord_index_dict[ortho] = {}
        for sp in template_species_list:
            for i in range(len(coords[ortho][sp])):
                coord = coords[ortho][sp][i]

                # filter for length
                start, end = coord
                length = end - start
                # if not min_exon_length <= length:
                if not min_exon_length <= length:
                    continue

                # filter for gap percent
                seq = str(aligned_fasta[ortho][sp].seq[start:end])
                if gapPercent(seq) > max_gap_percent:
                    continue

                # filter for gap length
                if longestGap(seq) > max_gap_length:
                    continue

                # populate coord_index_dict so I can find unique identifing info with ortho,sp,coord later
                if coord not in coord_index_dict[ortho].keys():
                    coord_index_dict[ortho][coord] = {}
                coord_index_dict[ortho][coord][sp] = i

                # prep to filter for species membership of ortho
                if coord not in ortho_coords[ortho].keys():
                    ortho_coords[ortho][coord] = set()
                ortho_coords[ortho][coord].add(sp)

    # set of coords per ortho which were represented in all species
    universal_ortho_coords = {}
    for ortho in ortho_coords:
        for coord in ortho_coords[ortho]:
            sp_set = ortho_coords[ortho][coord]
            if len(sp_set) == len(template_species_list):
                if ortho not in universal_ortho_coords.keys():
                    universal_ortho_coords[ortho] = set()
                universal_ortho_coords[ortho].add(coord)

    fasta_prep = {}
    for ortho in universal_ortho_coords:
        fasta_prep[ortho] = {}
        cds_list = {}
        for coord in universal_ortho_coords[ortho]:
            fasta_prep[ortho][coord] = {}
            for sp in template_species_list:
                if sp not in cds_list.keys():
                    parent = gff[sp][parent_groups[ortho][sp]]
                    cds_list[sp] = [cds for cds in
                                    gff[sp].children(parent, featuretype="CDS", order_by="start")]

                index = coord_index_dict[ortho][coord][sp]
                cds = cds_list[sp][index]
                fasta_prep[ortho][coord][sp] = cds

    nnn = Seq('n' * n_count, IUPAC.ambiguous_dna)
    shutil.rmtree(enhanced_alignment_path, ignore_errors=True)
    os.makedirs(enhanced_alignment_path, exist_ok=True)
    for ortho in fasta_prep:
        for coord in fasta_prep[ortho]:
            filename = enhanced_alignment_path + ortho + "_" + str(coord[0]) + "-" + str(coord[1]) + ".full.fasta"
            with open(filename, "w") as f:
                for sp in template_species_list:
                    cds = fasta_prep[ortho][coord][sp]
                    start, end = coord
                    seq = nnn + aligned_fasta[ortho][sp].seq[start:end] + nnn
                    seqReq = SeqRecord(seq, id=sp, description=cds.id)
                    f.write(seqReq.format("fasta"))

                for sp in sorted(parent_groups[ortho]):

                    if sp in template_species_list:
                        continue

                    parent = gff[sp][parent_groups[ortho][sp]]
                    strand = parent.strand
                    cds_list = [cds for cds in gff[sp].children(parent, featuretype="CDS", order_by="start")]
                    cat_seq = Seq("", IUPAC.ambiguous_dna)

                    for cds in cds_list:
                        try:
                            cat_seq += Seq(str(cds.sequence(fasta=fasta[sp], use_strand=False)), IUPAC.ambiguous_dna)
                        except ValueError as e:
                            if "imply a diffent length than sequence" in str(e):
                                cat_seq += Seq(str(fasta[sp][cds.chrom]), IUPAC.ambiguous_dna)
                                logging.debug(
                                    "coordinates from gff don't fall within with scaffold, grabbing entire scaffold")
                            else:
                                raise

                    if strand == '-':
                        cat_seq = cat_seq.reverse_complement()
                    seqRec = SeqRecord(cat_seq, id=sp, description=parent.id)
                    f.write(seqRec.format("fasta"))


def mafft_driver_file(file):
    p = subprocess.Popen(["./mafft_driver.sh", file, file + ".aln"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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


if __name__ == '__main__':
    import subprocess
    import argparse
    import pytoml

    parser = argparse.ArgumentParser()
    parser.add_argument('--configPath', help='configPath', default='../config.toml')
    args = parser.parse_args()

    # load config file
    with open(args.configPath) as toml_data:
        config = pytoml.load(toml_data)

    create_raw_exons(config['fasta_path'], config['enhanced_alignment_path'], config['template_alignment_path'],
                     config['db_path'],
                     config['json_path'], config['template_species_list'],
                     config['max_gap_percent'], config['max_gap_length'], config['min_exon_length'], config['n_count'])

    mafft_driver_path(config['enhanced_alignment_path'])
