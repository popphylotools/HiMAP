# orthoExon

## Description

Here you will find the bioinformatic locus selection and data processing pipelines used in Dupuis et al. 2017(link to paper eventually) for generating phylogenomic datasets using highly multiplexed amplicon sequencing. The code is divided into three main sections: 

**Part01** is the main locus selection pipeline, which takes as input ortholog prediction from a variety of genomic and transcriptomic resources (some of which with relatively trustworthy structural annotations: "high quality annotations"). Exon/intron boundaries from the "high quality annotations" are used to predict exon/intron boundaries across all data, and several filtering steps identify conserved exons across data inputs. Dupuis et al. 2017 uses [OrthoMCL](http://orthomcl.org/orthomcl/) for ortholog prediction.

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
curl "https://codeload.github.com/popphylotools/orthoExon/zip/master" -o "orthoExon-master.zip"
unzip orthoExon-master.zip
cd orthoExon-master
```


#### Download Data

Input, intermediate, and output data for each subsection (part) avalible at:
```
curl "http://67.52.95.73/~woods26/orthoExon_data.zip" -o "orthoExon_data.zip"
unzip orthoExon_data.zip
mv data data.bak
mv orthoExon_data data
```

orthoExon_data.zip can be unziped and used to replace the data directory from this git repo.

## Usage: Quick and Dirty
Configure data paths by editing `config.toml`.
As configured, it will run with the example data linked above.

Old output of each step is deleted before new output files are written, therefore all `intermediate` and `output` directories are optional when starting from the begining of any subsection (part).

Activate main orthoExon env with:
```
source activate orthoExon
```

Enter the appropriate directory for the desired subsection:
```
cd part01
```

Run the desired script:
```
./step01_create_padded_cds_and_align.py
```

## Usage: Detailed
#### config.toml
The file `config.toml` contains various parameters and data paths for running all three parts of this pipeline. Notably, this includes 2 lists of species names (which match the output of the ortholog prediction): one of all species, and one of species with "high quality annotations". The latter are used to predict exon/intron boundaries for the remaining species. This config file also contains paths to the various inputs/intermediates/outputs for each of the three parts. Only the inputs directories are required to run the scripts (intermediate & output dirs will be generated automatically). If using all steps of this pipeline, we suggest using a data structure such as this:

```
.
+-- _data
|   +-- _part01
|       +-- _input
|           +-- [ortholog prediction results].txt
|           +-- _gff
|               +-- [input gff files]
|           +-- _fasta
|               +-- [input fasta files]
|   +-- _part02
|       +-- _input
|   +-- _part03
|       +-- _input
```

#### Part01
Part01 takes as input predicted single copy orthologs (and their corresponding genomic/transcriptomic data files), and outputs conserved exons. [UNDER CONSTRUCTION]



## Uninstall

If needed, remove Anaconda environments
```
conda env remove --name orthoExon
conda env remove --name mafft
conda env remove --name tapir
```

## Credits

 * [Forest Bremer](https://github.com/Woods26):</br>
   primary programmer and github maintainer
 * [Julian Dupuis](https://github.com/jrdupuis):</br>
   primary author of companion paper and github documentation contributer
