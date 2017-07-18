#!/bin/bash

# move nex files to nex directories
I=0
for file in $(/bin/ls ./fasta/*.nex)
do
    mv ${file} ./nex$(($I % 32))/
    I=$(($I + 1))
done
