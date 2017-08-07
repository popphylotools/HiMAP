orthoExon
=========

Here you will find the bioinformatic locus selection and data processing pipelines used in Dupuis et al. 2017(link to paper eventually) for generating phylogenomic datasets using highly multiplexed amplicon sequencing. The code is divided into three main sections: 

Part01 is the main locus selection pipeline, which takes as input ortholog prediction from a variety of genomic and transcriptomic resources (some of which with relatively trustworthy structural annotations: "high quality annotations"). Exon/intron boundaries from the "high quality annotations" are used to predict exon/intron boundaries across all data, and several filtering steps identify conserved exons across data inputs.

Part02 generates summary information used to manually filter and select amplicons to include in a final amplicon panel, based on amplicon length, primer characteristics, phylogenetic informativeness, etc. In Dupuis et al. 2017, the input to this step is predicted primers generated from the Paragon Genomics CleanPlex custom amplicon design service.

Part03 is for post-sequencing data processing. This takes as input demultiplexed, adapter-trimmed and FLASh-joined FASTQ files, and calls consensus sequences based on read length distributions. It outputs aligned multi-FASTA formatted files (1 per gene), that can be used directly for gene-tree analysis or concatenated for an "all data" approach.

For a detailed tutuorial in re-creating the results of Dupuis et al. 2017, see here(still working on this).


installation 
------------

 - install Anaconda
 - create Anaconda environments with:
  - `conda env create -f create_anaconda_environments.sh`
 - activate env with:
  - linux/osx:
   - `source activate orthoExon`
 - deactivate env with
  - linux/osx:
   - `source deactivate orthoExon`
 - if needed, remove with:
  - `conda env remove --name orthoExon`

data
----
imput and intermediate data avalible at:
http://67.52.95.73/~woods26/
