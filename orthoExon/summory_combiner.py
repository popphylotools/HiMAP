#!/usr/bin/env python

import csv
import shutil
from collections import Counter
from functools import reduce
from operator import mul

import config

possibilities = {key: len(value) for key, value in config.expand_iupac.items()}


def amb_count(seq):
    return sum([value for key, value in Counter(seq).items() if key in "yrwskmdvhbn-"])


def amb_combinations(seq):
    return reduce(mul, [possibilities[nuc.upper()] for nuc in seq])


def summory_combiner(filename_list):
    data = []
    for i in range(len(filename_list)):
        with open(filename_list[i]) as f:
            lines = f.readlines()
            lines = [line.strip().split(",") for line in lines[1:]]
            lines = {line[0]: (line[0][:-3],
                               int(line[0][-1:]),
                               i,
                               amb_count(line[5] + line[6]),
                               amb_count(line[5]),
                               amb_count(line[6]),
                               max(amb_count(line[5]), amb_count(line[6])),
                               amb_combinations(line[5]) + amb_combinations(line[6]),
                               amb_combinations(line[5]),
                               amb_combinations(line[6]),
                               float(line[1]),
                               int(line[4]),
                               line[5],
                               line[6],
                               float(line[7]),
                               float(line[8])) for line in lines}

            data += [rec for rec in lines.values()]

    data = sorted(data, key=lambda x: x[0])

    shutil.rmtree('combined_summory.csv', ignore_errors=True)
    with open('combined_summory.csv', "w") as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
        header = ('CDS_Ortholog',
                  'Primer_Version',
                  'Ambiguities_Allowed(p3_input)',
                  'Total_Ambiguities(L+R)',
                  'Left_Ambiguities',
                  'Right_Ambiguities',
                  'Max_Ambiguities(L|R)',
                  'Total_Combinations(L+R)',
                  'Left_Combinations',
                  'Right_Combinations',
                  'PI_Score',
                  'Target_Sequence_Length',
                  'PRIMER_LEFT_0_SEQUENCE',
                  'PRIMER_RIGHT_0_SEQUENCE',
                  'PRIMER_LEFT_0_TM',
                  'PRIMER_RIGHT_0_TM')
        writer.writerow(header)
        writer.writerows(data)


if __name__ == '__main__':
    filename_list = ["../output{}/summory.csv".format(i) for i in range(3)]

    summory_combiner(filename_list)
