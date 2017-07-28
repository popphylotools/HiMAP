#!/bin/bash
shopt -s expand_aliases
module load tapir
module load parallel

DIRECTORY=${1:-"../intermediate/phylo_informativeness"}

/bin/ls $DIRECTORY/fasta/*.fasta | parallel 'perl ./convertfasta2nex.pl {} > {}.nex'

# remove and recreate nex directories
rm -r $DIRECTORY/nex*

for i in {0..31}
do
    mkdir $DIRECTORY/nex$i/
done

# move nex files to nex directories
I=0
for file in $(/bin/ls $DIRECTORY/fasta/*.nex)
do
    mv $file $DIRECTORY/nex$(($I % 32))/
    I=$(($I + 1))
done

# remove and recreate tapir ouput directories
rm -r $DIRECTORY/tapir_out*

for i in {0..31}
do
    mkdir $DIRECTORY/tapir_out$i/
done

#tapir_compute.py $DIRECTORY/nex/ ../input/959genes.phy.contree_AgeRooted_wholeNumbs.tre --intervals=2-3,3-5,5-8,8-10,10-18 --times=2,3,5,8,10,17 --tree-format=newick --output $DIRECTORY/tapir_out

seq 0 31 | parallel tapir_compute.py $DIRECTORY/nex{}/ ./959genes.phy.contree_AgeRooted_wholeNumbs.tre --intervals=2-3,3-5,5-8,8-10,10-18 --times=2,3,5,8,10,17 --tree-format=newick --output $DIRECTORY/tapir_out{}

# collect output
rm -r $DIRECTORY/sql
mkdir $DIRECTORY/sql

for i in {0..31}
do
    mv $DIRECTORY/tapir_out$i/*.sqlite $DIRECTORY/sql/$i.sqlite
done

