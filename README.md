# HiMAP: Highly Multiplexed Amplicon-based Phylogenomics

## Description

Here you will find the bioinformatic locus selection and data processing pipelines used in Dupuis et al. 2017(link to paper eventually) for generating phylogenomic datasets using the HiMAP approach. The code is divided into three main sections: 

**Part01** is the main locus selection pipeline, which takes as input ortholog prediction from a variety of genomic and transcriptomic resources (some of which with relatively trustworthy structural annotations: "high quality annotations"). Exon/intron boundaries from the "high quality annotations" are used to predict exon/intron boundaries across all data, and several filtering steps identify conserved exons across data inputs. Dupuis et al. 2017 uses [OrthoMCL](http://orthomcl.org/orthomcl/) for ortholog prediction, and this code is built to ingest OrthoMCL input format, however other ortholog prediction approaches can also be used.

**Part02** generates summary information used to manually filter and select amplicons to include in a final amplicon panel, based on amplicon length, primer characteristics, phylogenetic informativeness, etc. In Dupuis et al. 2017, the input to this step is predicted primers generated from the Paragon Genomics CleanPlex custom amplicon design service.

**Part03** is for post-sequencing data processing. This takes as input demultiplexed, adapter-trimmed and FLASh-joined FASTQ files, and calls consensus sequences based on read length distributions. It outputs aligned multi-FASTA formatted files (1 per gene), that can be used directly for gene-tree analysis or concatenated for an "all data" approach.

For a detailed tutuorial in re-creating the results of Dupuis et al. 2017, see **Usage** section below. Note, the terminology from Dupuis et al. 2017 is used here ("raw exons", "filtered exons", etc.), so refer to the paper for those details.


## Installation

#### Install Anaconda distribution of python

[Anaconda download page](https://www.continuum.io/downloads)

#### Download this git repo
```
curl "https://codeload.github.com/popphylotools/HiMAP/zip/master" -o "HiMAP-master.zip"
unzip HiMAP-master.zip
cd HiMAP-master
```

#### Create Anaconda environments

osx:
```
./osx_create_anaconda_environments.sh
```

linux:
```
./linux_create_anaconda_environments.sh
```

#### Download Data

Input, intermediate, and output data for each subsection (part) avalible at:
```
curl "http://67.52.95.73/~woods26/HiMAP_data.zip" -o "HiMAP_data.zip"
unzip HiMAP_data.zip
mv data data.bak
mv HiMAP_data data
```

HiMAP_data.zip can be unzipped and used to replace the data directory from this git repo.

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
...
```
In this example, each ortholog is on a row, followed by the species ("Afra", "Aobl", "Asus", "Bcor", etc.) and their ortholog accession (after the |). FASTA (in peptide and nucleotide space) and GFF files for the original data sources are also required in additional directories, FASTA headers need to match the Sequence ID of their respecitve GFF files (first column of GFF). Note that Dupuis et al. (2017) conducted 

Part01 is run through 5 scripts, that can be run in direct succession:
```
./step00a_create_gff_db.py
./step00b_create_parent_groups_json.py
./step01_create_padded_cds_and_align.py
./step02_create_raw_exons_and_align.py
./step03_create_filtered_exons.py
```
Details of the individual scripts can be found in the part01 [README](https://github.com/popphylotools/HiMAP/tree/master/part01_find_ortho_exons). The output of part01 is multi-FASTA formatted conserved exons, and these files are named with the ortholog ID followed by position coordinates for that exon: `orth4955_1027-1296`. These coordinates refer to the positions in the [WHICH IS THIS?].





#### Part02
Part02 ingests putative primer information, calculates phylogenetic informativeness and other amplicon summary information, and outputs a summary.csv file containing all summary information. The goal of this process is to generate a file that can then be manually viewed and sorted either on the command line or in Excel or another spreadsheet editing platform. Dupuis et al. (2017) uses Paragon Genomics CleanPlex custom amplicon design service, and thus follow the output format used therein:
```
amplicon_ID	loci_ID	amp_start	amp_end	ampInsert_start	ampInsert_end	amp_len	comment
SET3976_3049	orth4955_1027-1296	9	258	39	232	250	NoDegenerate
SET1560_252046	orth3516_2121-3096	372	605	394	583	234	NoDegenerate
SET4061_2820	orth5003_8684-8996	15	236	38	213	222	NoDegenerate
SET10422_95779	orth9368_146-542	155	373	181	346	219	NoDegenerate
SET9146_445206	orth8353_0-1374	641	858	663	833	218	NoDegenerate
...
```
This format can be easily created from other primer generating software. If creating this file from another primer generation approach, `loci_ID` must match the amplicon file names, and `amplicon_ID` can be an arbitrary unique code. `amp_start/end` and `ampInsert_start/end` essentially give coordinates of the beginning and end of each primer set, with `amp` referring to the entire amplicon and `ampInsert` referring to the amplicon excluding the primers.

Part02 calculates phylogenetic informativeness (PI) for each amplicon (specifically the `ampInsert`, excluding the primers) using [TAPIR](https://github.com/faircloth-lab/tapir). This requires a dated reference tree for the species included in the original ortholog prediction, in NEWICK format. Installation of TAPIR is automatic when using the Anaconda Environments outlined above. By default, this script averages the PI values at specified times along the reference tree (the `--times` parameter in TAPIR's `compute.py`). These times need to be specified by the user, and can be edited in `part02_score_primers_from_paragon/tapir_driver.sh` found [here](https://github.com/popphylotools/HiMAP/tree/master/part02_score_primers_from_paragon/tapir_driver.sh). The source code can also be edited to output different estimates generated by TAPIR.

Part02 is ran through a single script:
```
./step04_make_primer_summary.py
```
And generates an output `summary.csv` that is formatted like this:
```
loci_id,amplicon_id,score,total_ambiguities(l+r),max_ambiguities(l|r),left_ambiguities,right_ambiguities,total_primer_versions,l_primer_versions,r_primer_versions,total_combinations(l+r),left_combinations,right_combinations,insert_length,ampinsert_start,ampinsert_end,amp_len,amp_start,amp_end,left_primer,right_primer,comment,species
orth10136_1834-2121,SET77_21882,0.5639935823860325,2,1,1,1,4,2,2,4,2,2,97,88,185,140,69,208,GyTGTCCCATCTGTGACAA,CAGTCGCACCTAACkCAACATTT,DegenerateNeeded,"['Bcur', 'Bdor', 'Bole', 'Ccap', 'Afra', 'Asus', 'Bcor', 'Bzon']"
orth10223_1638-1932,SET120_9638,0.3953473750685413,2,1,1,1,5,2,3,5,2,3,120,51,171,155,34,188,AGTAyGCTCGCACCGCC,ATGACCACbGTCTGGGC,DegenerateNeeded,"['Bcur', 'Bdor', 'Bole', 'Ccap', 'Bcor', 'Blat', 'Bmin', 'Btry']"
orth10262_1395-1632,SET130_19936,0.38906508363165077,2,1,1,1,5,2,3,5,2,3,96,93,189,137,71,207,GGTCATCGTTTyGTGAATGTGG,TGGATGATGTGAAGGCbG,DegenerateNeeded,"['Bcur', 'Bdor', 'Bole', 'Ccap', 'Afra', 'Bcor', 'Blat', 'Bzon']"
orth10262_574-1105,SET131_138012,0.8773904790654025,2,1,1,1,4,2,2,4,2,2,212,213,425,251,196,446,CGGAGGTGCGyAGCTTT,GATTTGyCATCGGGATGGGGT,DegenerateNeeded,"['Bcur', 'Bdor', 'Bole', 'Ccap', 'Afra', 'Bcor', 'Blat', 'Bzon']"
orth10297_2517-2710,SET155_248,0.38186248775817155,2,1,1,1,4,2,2,4,2,2,126,34,160,172,10,181,TGCGAAACTCACCAAAGGAwATTG,CGGTkGAAAGTTTACCCAGCG,DegenerateNeeded,"['Bcur', 'Bdor', 'Bole', 'Ccap', 'Afra', 'Blat', 'Bmin', 'Btry']"
orth10315_1281-1727,SET164_65033,0.4744932094613006,3,2,2,1,5,3,2,6,4,2,161,146,307,206,129,334,CmGCCGATCGTGATGAr,CATATTGAAGAAGAAACCCATGCAyGA,DegenerateNeeded,"['Bcur', 'Bdor', 'Bole', 'Ccap', 'Bcor', 'Blat', 'Btry', 'Bzon']"
orth10339_505-906,SET187_90989,0.11472834404492394,1,1,1,0,3,2,1,3,2,1,158,153,311,199,134,332,GTTAGCGACGGGTATGCTm,CTATGTTGCGCTCAGTTTGGT,NoDegenerate,"['Bcur', 'Bdor', 'Bole', 'Ccap', 'Bcor', 'Blat', 'Btry', 'Bzon']"
orth10362_263-601,SET193_30850,0.10833524405700566,0,0,0,0,2,1,1,2,1,1,113,86,199,159,63,221,CTTCTGCGATGATGTGTCAATCC,CTCATCACTCATATCGGAGCGT,NoDegenerate,"['Bcur', 'Bdor', 'Bole', 'Ccap', 'Bcor', 'Blat', 'Btry', 'Bzon']"
orth10375_1062-1241,SET202_5619,0.14771179397003442,2,1,1,1,4,2,2,4,2,2,86,59,145,142,26,167,ACTGCCATTAAATGTTGTTAATGAGAAAATwTA,GyATGAGGGTGTACTTGTTTCG,DegenerateNeeded,"['Bcur', 'Bdor', 'Bole', 'Ccap', 'Bcor', 'Blat', 'Btry', 'Bzon']"
...
```

Many of the columns in `summary.csv` match those in the putative primer information input (e.g. `loci_id`, `amplicon_id`). `score` refers to PI score, and [INCLUDE MORE INFO ABOUT PRIMER COLUMNS].

Following the general HiMAP approach, this output is then used to select final amplicons based on primer characteristics (e.g. number of degenerate bases), phylogenetic informativeness, etc. See Dupuis et al. (2017) for an example of this filtering process.


#### Part03
Part03 handles the post-sequencing data processing, and calling consensus sequences....[UNDER CONSTRUCTION]



MAKE SURE TO INCLUDE name splitting by orth, so sample names coming out of cutadapt (second time) needs to be [ind]orth[orth#]. Code splits by "orth", and 








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
