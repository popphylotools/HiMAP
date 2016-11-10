#!/usr/bin/env python
# coding: utf-8

import os
import gffutils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import json
from pyfaidx import Fasta
from Bio import SeqIO
import re
from config import template_species_list, transvestigated_species_set, fasta_path, orthoCds_path, \
    fullset_alignment_path, db_path, json_path

# globals
# template_species_list = ["Bcur", "Bdor", "Bole", "Ccap"]
# transvestigated_species_set = {'Bcor', 'Blat', 'Bzon', 'Afra', 'Bmin', 'Bjar', 'Aobl'}
# fasta_path = "../input/fasta/"
# orthoCds_path = "../output/orthoCds/"
# fullset_alignment_path = "../intermediate/13spp_alignment/"
# db_path = "../intermediate/gff_databases/"
# json_path = "../intermediate/json/"

# gap filter parameters
max_gap_percent = 0
max_gap_length = 0
# cds length filter parameters
min_cds_length = 100
# max_cds_length = 600

# create handles for all .db files in intermediate directory
gff_fn = {name.split('.gff.db')[0]: db_path + name for name in os.listdir(db_path) if
          ".gff.db" in name}
gff = {key: gffutils.FeatureDB(value) for key, value in gff_fn.items()}

# create handles for all .fasta files in fasta directory
fasta_fn = {name.split('.nt.fasta')[0]: fasta_path + name for name in os.listdir(fasta_path) if
            ((".nt.fasta" in name) and (".nt.fasta.fai" not in name))}
fasta = {key: Fasta(value) for key, value in fasta_fn.items()}

# import ortholog groups
with open(json_path + "groups.json", 'r') as f:
    parent_groups = json.load(f)

# create handles for all .fasta files in aligned_13spp_fasta directory
aligned_fasta_fn = {name.split('.13spp')[0]: fullset_alignment_path + name for name in
                    os.listdir(fullset_alignment_path) if
                    ((".fasta.aln" in name) and (".fasta.aln.fai" not in name))}

# define functions to parse coordinates of cds's from concatinated aligned fasta w/ n's and -'s
nnn = 50


def findBreakpoints(seq):
    breakpoints = []
    loc = 0
    regex = re.compile(r"n+[-+n+]*")
    while True:
        # loc = seq.find(nnn, loc)
        match = regex.search(seq, loc)
        if not match:
            break
        if len(match.group().replace('-', '')) >= nnn:
            breakpoints.append(match.span())
        loc = match.end()
    return breakpoints


def findExonCoords(seq):
    breakpoints = findBreakpoints(seq)
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


# read and parse fasta files for each species
aligned_fasta = {}
for ortho in aligned_fasta_fn.keys():
    aligned_fasta[ortho] = {seq_record.id: seq_record
                            for seq_record in SeqIO.parse(aligned_fasta_fn[ortho],
                                                          "fasta", alphabet=IUPAC.ambiguous_dna)}

# parse coords from template species in aligned fasta's and trash entries w/ all gaps
coords = {}  # coords[ortho][sp] = [coord, ]
for ortho in aligned_fasta:
    coords[ortho] = {}
    for sp in template_species_list:
        seq = str(aligned_fasta[ortho][sp].seq)
        temp_coords = findExonCoords(str(aligned_fasta[ortho][sp].seq))
        for start, end in temp_coords:
            cds = seq[start:end]
            if len(cds) != cds.count('-'):
                if sp not in coords[ortho]:
                    coords[ortho][sp] = (start, end)
                elif type(coords[ortho][sp]) is list:
                    coords[ortho][sp].append((start, end))
                else:
                    temp = coords[ortho][sp]
                    coords[ortho][sp] = [temp, (start, end)]

# sanity check for multiple non gap template cds's per ortho,sp
for ortho in coords:
    for sp in coords[ortho]:
        if type(coords[ortho][sp]) is list:
            print("error, multiple non-gap template cds's for {},{}: {}".format(ortho, sp,
                                                                                coords[ortho][sp]))

# Filter aligned exons
ortho_coords = {}
for ortho in coords:
    ortho_coords[ortho] = {}
    for sp in coords[ortho]:
        coord = coords[ortho][sp]

        # filter for length
        start, end = coord
        length = end - start
        # if not min_cds_length <= length <= max_cds_length:
        if not min_cds_length <= length:
            continue

        # filter for gap percent
        seq = str(aligned_fasta[ortho][sp].seq[start:end])
        if gapPercent(seq) > max_gap_percent:
            continue

        # filter for gap length
        if longestGap(seq) > max_gap_length:
            continue

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
        else:
            print("warning, {} {} has only {}".format(ortho, coord, sp_set))

# fasta prep
fasta_prep = {}
for ortho in universal_ortho_coords:
    fasta_prep[ortho] = []
    for coord in universal_ortho_coords[ortho]:
        temp_sp_list = []
        for sp in sorted(aligned_fasta[ortho]):
            start, end = coord
            seq = aligned_fasta[ortho][sp].seq[start:end]
            des = aligned_fasta[ortho][sp].description
            seqReq = SeqRecord(seq, id=sp, description=des)
            if sp in template_species_list:
                fasta_prep[ortho].append(seqReq)
            else:
                temp_sp_list.append(seqReq)

        fasta_prep[ortho].extend(temp_sp_list)

for ortho in fasta_prep:
    fasta_prep[ortho] = [seqReq for seqReq in fasta_prep[ortho] if
                         (gapPercent(seqReq.seq) <= max_gap_percent) and (
                             longestGap(seqReq.seq) <= max_gap_length)]

fasta_prep = {ortho: seq_list for ortho, seq_list in fasta_prep.items() if len(seq_list) >= 8}

# fasta output
for ortho in fasta_prep:
    filename = orthoCds_output + ortho + ".13spp.fasta"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "w") as f:
        for seqReq in fasta_prep[ortho]:
            f.write(seqReq.format("fasta"))
