# HiMAP: Highly Multiplexed Amplicon-based Phylogenomics

## Description

Here you will find the bioinformatic locus selection and data processing pipelines used in Dupuis et al. 2017(link to paper eventually) for generating phylogenomic datasets using the HiMAP approach. The code is divided into three main sections: 

**Part01** is the main locus selection pipeline, which takes as input ortholog prediction from a variety of genomic and transcriptomic resources (some of which with relatively trustworthy structural annotations: "high quality annotations"). Exon/intron boundaries from the "high quality annotations" are used to predict exon/intron boundaries across all data, and several filtering steps identify conserved exons across data inputs. Dupuis et al. 2017 uses [OrthoMCL](http://orthomcl.org/orthomcl/) for ortholog prediction, and this code is built to ingest OrthoMCL input format, however other ortholog prediction approaches can also be used.

**Part02** generates summary information used to manually filter and select amplicons to include in a final amplicon panel, based on amplicon length, primer characteristics, phylogenetic informativeness, etc. In Dupuis et al. 2017, the input to this step is predicted primers generated from the Paragon Genomics CleanPlex custom amplicon design service.

**Part03** is for post-sequencing data processing. This takes as input demultiplexed, adapter-trimmed and FLASh-joined FASTQ files, and calls consensus sequences based on read length distributions. It outputs aligned multi-FASTA formatted files (1 per gene), that can be used directly for gene-tree analysis or concatenated for an "all data" approach.

For a detailed tutuorial in re-creating the results of Dupuis et al. 2017, see **Usage** section below.


## Installation

#### Install Anaconda distribution of python

[Anaconda download page](https://www.continuum.io/downloads)

#### Create Anaconda environments

osx:
```
./osx_create_anaconda_environments.sh
```

linux:
```
./linux_create_anaconda_environments.sh
```

#### Download this git repo
```
curl "https://codeload.github.com/popphylotools/HiMAP/zip/master" -o "HiMAP-master.zip"
unzip HiMAP-master.zip
cd HiMAP-master
```


#### Download Data

Input, intermediate, and output data for each subsection (part) avalible at:
```
curl "http://67.52.95.73/~woods26/HiMAP_data.zip" -o "HiMAP_data.zip"
unzip HiMAP_data.zip
mv data data.bak
mv HiMAP_data data
```

HiMAP_data.zip can be unziped and used to replace the data directory from this git repo.

## Usage: Quick and Dirty
Configure data paths by editing `config.toml`.
As configured, it will run with the example data linked above.

Old output of each step is deleted before new output files are written, therefore all `intermediate` and `output` directories are optional when starting from the begining of any subsection (part).

Activate main HiMAP env with:
```
source activate HiMAP
```

Enter the appropriate directory for the desired subsection:
```
cd part01_find_ortho_exons
```

Run the desired script:
```
./step01_create_padded_cds_and_align.py
```

## Usage: Detailed
#### config.toml
The file `config.toml` contains various parameters and data paths for running all three parts of this pipeline. Notably, this includes 2 lists of species names (spelled as in the output of the ortholog prediction): one of all species, and one of species with "high quality annotations". The latter are used to predict exon/intron boundaries for the remaining species. This config file also contains paths to the various inputs/intermediates/outputs for each of the three parts. Only the inputs directories are required to run the scripts (intermediate & output dirs will be generated automatically). If using all steps of this pipeline, we suggest using a data structure such as the following and refer to this directory structure below:

```
.
|-- data
|   |-- part01
|       |-- input
|           |-- [ortholog prediction results].txt
|           |-- gff
|               |-- [input gff files]
|           |-- fasta
|               |-- [input fasta files]
|   |-- part02
|       |-- input
|   |-- part03
|       |-- input
```

#### Part01
Part01 takes as input predicted single copy orthologs (and their corresponding genomic/transcriptomic data files), and outputs conserved exons. Ortholog prediction results must conform to OrthoMCL input format:
```
orth2333: Afra|m.10026 Aobl|m.12854 Asus|M.10026_R0 Bcor|m.49852 Bcur|XP_011176775.1 Bdor|XP_011207020.1 Bjar|m.26921 Blat|m.41225 Bmin|m.16652 Bole|XP_014097011.1 Btry|M.41225_R0 Bzon|m.21935 Ccap|XP_004523616.1
orth2334: Afra|m.10105 Aobl|m.12162 Asus|M.10104_R0 Bcor|m.20328 Bcur|XP_011189361.1 Bdor|XP_011202495.1 Bjar|m.4341 Blat|m.39156 Bmin|m.5536 Bole|XP_014100885.1 Btry|M.39156_R0 Bzon|m.47619 Ccap|XP_004521443.1
orth2335: Afra|m.10131 Aobl|m.7135 Asus|M.10131_R0 Bcor|m.509 Bcur|XP_011191044.1 Bdor|XP_011202275.1 Bjar|m.10381 Blat|m.18775 Bmin|m.8724 Bole|XP_014102377.1 Btry|M.18775_R0 Bzon|m.34590 Ccap|XP_004518888.1
```
In this example, each ortholog is on a row, followed by the species ("Afra", "Aobl", "Asus", "Bcor", etc.) and their ortholog accession (after the |). FASTA (in peptide and nucleotide space) and GFF files for the original data sources are also required in additional directories. Part01 is run through 5 scripts, that can be run in direct succession:
```
./step00a_create_gff_db.py
./step00b_create_parent_groups_json.py
./step01_create_padded_cds_and_align.py
./step02_create_raw_exons_and_align.py
./step03_create_filtered_exons.py
```
Details of the individual scripts can be found in the part01 [README](https://github.com/popphylotools/HiMAP/tree/master/part01_find_ortho_exons).







## Uninstall

If needed, remove Anaconda environments
```
conda env remove --name HiMAP
conda env remove --name mafft
conda env remove --name tapir
```

## Credits

 * [Forest Bremer](https://github.com/Woods26):</br>
   primary programmer and github maintainer
 * [Julian Dupuis](https://github.com/jrdupuis):</br>
   primary author of companion paper and github documentation contributer
