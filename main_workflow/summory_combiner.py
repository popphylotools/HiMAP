#!/usr/bin/env python

# coding: utf-8

# In[1]:

import csv
from collections import Counter


# In[2]:

fn = ["../output{}/summory.csv".format(i) for i in range(3)]


# In[3]:

header = ('Exon_Name',
 'Ambiguities_Allowed',
 'Total_Ambiguities(L+R)',
 'Left_Ambiguities',
 'Right_Ambiguities(L+R)',
 'Max_Ambiguities(L|R)',
 'PI_Score',
 'Target_Sequence_Length',
 'PRIMER_LEFT_0_SEQUENCE',
 'PRIMER_RIGHT_0_SEQUENCE',
 'PRIMER_LEFT_0_TM',
 'PRIMER_RIGHT_0_TM')


# In[4]:

def amb_count(seq):
    return sum([value for key,value in Counter(seq).items() if key in "yrwskmdvhbn-"])


# In[5]:

data = []
for i in range(len(fn)):
    with open(fn[i]) as f:
        lines = f.readlines()
        lines = [line.strip().split(",") for line in lines[1:]]
        lines = {line[0]:(line[0], str(i), str(amb_count(line[6] + line[7])), str(amb_count(line[6])), str(amb_count(line[7])), str(amb_count(line[6]) + amb_count(line[7])), line[1], line[4], line[5], line[6], line[7], line[8]) for line in lines}
        data = data + [rec for rec in lines.values()]


# In[6]:

data = sorted(data, key=lambda x: x[0])


# In[7]:

ofile  = open('combined_summory.csv', "w")
writer = csv.writer(ofile, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)

writer.writerow(header)
for row in data:
    writer.writerow(row)

