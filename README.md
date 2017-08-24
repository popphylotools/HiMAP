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

Input, intermediate, and output data for each subsection (part) avalible at: (to be replaced with Dryad accession)
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

Run the desired script, e.g.:
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
In this example, each ortholog is on a row, followed by the species ("Afra", "Aobl", "Asus", "Bcor", etc.) and their ortholog accession (after the |). FASTA (in peptide and nucleotide space) and GFF files for the original data sources are also required in additional directories, and FASTA headers need to match the Sequence ID of their respecitve GFF files (first column of GFF). 

Part01 is run through 5 scripts, that can be run in direct succession:
```
./step00a_create_gff_db.py
./step00b_create_parent_groups_json.py
./step01_create_padded_cds_and_align.py
./step02_create_raw_exons_and_align.py
./step03_create_filtered_exons.py
```
Details of the individual scripts can be found in the part01 [README](https://github.com/popphylotools/HiMAP/tree/master/part01_find_ortho_exons). The output of part01 is multi-FASTA formatted conserved exons, and these files are named with the ortholog ID followed by position coordinates for that exon: `orth4955_1027-1296`. These coordinates refer to the positions in the "padded exons" (before orthologs are split up into multiple exons). Note, that part03 (which calls the final consensus sequences, post-sequencing) requires that amplicon names begin with "orth", so this part of the name should be preserved.

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

Part02 calculates phylogenetic informativeness (PI) for each amplicon (specifically the `ampInsert`, excluding the primers) using [TAPIR](https://github.com/faircloth-lab/tapir). This requires a dated reference tree for the species included in the original ortholog prediction, in NEWICK format. Installation of TAPIR is automatic when using the Anaconda Environments outlined above. By default, this script averages the PI values at specified times along the reference tree (the `--times` parameter in TAPIR's `compute.py`). These times need to be specified by the user, and can be edited in `config.toml`. The source code can also be edited to output different estimates generated by TAPIR, if needed.

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

Many of the columns in `summary.csv` match those in the putative primer information input (e.g. `loci_id`, `amplicon_id`, `amp_len`, `insert_length`). `score` refers to PI score, and columns dealing with `ambiguities`, `primer_versions`, and `combinations` are various statistics of the ambiguities present in given primers, and how many individual primer sequences are required given those ambiguities (e.g. a `R` in a primer would lead to 2 "versions" or sequences for that primer, and two potential combinations with a non-degerate partner primer).

Following the general HiMAP approach, this output is then used to select final amplicons based on primer characteristics (e.g. number of degenerate bases), phylogenetic informativeness, etc. See Dupuis et al. (2017) for an example of this filtering process.


#### Part03
Part03 handles the post-sequencing data processing, and calling consensus sequences. The input for this script is adapter-trimmed and demultiplexed (by individual) FASTQ files. We suggest using [cutadapt](http://cutadapt.readthedocs.io/en/stable/index.html) and [FLASh](https://ccb.jhu.edu/software/FLASH/) for these steps, and provide details of the filtering done for Dupuis et al. (2017) below:

First, reads need to be demultiplexed by individual. If sequencing is done on an Illumina MiSeq or HiSeq connected to BaseSpace, this demultiplexing may be done on BaseSpace. Alternatively, cutadapt can be used to demultiplex based on individual-specific barcodes.

Next, cutadapt can be used to remove any additional Illumina adapters. We will assume that the raw FASTQ files are in `./RawData/`, and cutadapt can be run like:
```
for x in `cat individuals`; do echo "$x" |  cutadapt -a agatcggaagagcacacgtctgaa -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o AdapterTrimmed/"$x"_R1_adaptertrimmed.fastq -p AdapterTrimmed/"$x"_R2_adaptertrimmed.fastq  RawData/"$x"_R1_combined.fastq.gz RawData/"$x"_R2_combined.fastq.gz > CutAdapt_logs/log."$x".log ; done
```

Then, FLASh can be used to join the paired reads for all data:
```
for x in `cat individuals` ; do flash AdapterTrimmed/"$x"_R1_adaptertrimmed.fastq AdapterTrimmed/"$x"_R2_adaptertrimmed.fastq -o Flash/"$x"_flash.fastq | tee Flash/"$x".log ; done
```

Finally, cutadapt can be used again, but this time to demultiplex each individual, FLASh-joined FASTQ file by amplicon. Dupuis et al. (2017) pooled 384 individuals and 878 amplicons into single sequencing lanes, so we used a job array on a cluster to speed this process up. The job file looked something like this:
```
#!/bin/sh
#$-S /bin/sh
#$-o Demultiplexed_CutAdapt_Logs/$JOB_ID.out 
#$-e Demultiplexed_CutAdapt_Logs/$JOB_ID.err
#$-q all.q
#$-pe orte 1
#$-cwd
#$-V
#$-t 1-384
INFILE=`awk "NR==$SGE_TASK_ID" individuals`
echo "$INFILE"
module load cutadapt
cutadapt -O 10 -a orth10028_611-1041=TGCCCATCGCCTCCGAy...GAGGTGTACTTGGTGGGCG \
-a orth10034_825-2024=GGCACCACATTCTCACAm...TGTATCAACTGCAGGCGAG \
-a orth10118_500-1241=TGTGGCTACTCGTGTCGw...CCACATGAATGTAAAGTATGTGGACG \
... .... ...
-a orth9938_2137-2366=AGCsAAGGTGCAACAAGTCT...GGCGATCGTCGGGATCA \
-a orth9941_1692-2530=AAGAGAGCAACCCACCTr...AGCCATGGAACTCGCCAA \
-o Demultiplexed_CutAdapt/"$INFILE"/"$INFILE".{name}.fastq Flash/"$INFILE"_flash.fastq.extendedFrags.fastq | tee Demultiplexed_CutAdapt_Logs/"$INFILE".log
```
Here, the `-a` options specify the amplicon primer pairs, so this job file would contain 878 `-a` options. We also split the output into multiple directories, one per individual, which is required for the part03 script that calls consensus sequences. Splitting the files into multiple directories is also a good idea, as 384 individuals x 878 amplicons = 337,152 files at the end of this step (which could stall bash commands or general command line manipulation). This process could also be run through a for loop, with modified variables in the input/output paths of the cutadapt command.

Following these steps, the sequencing data is now demultiplexed by individual and amplicon, and the data structure should be one main directory (`Demultiplexed_CutAdapt`) containing a single directory for each individual; each individual directory contains one FASTQ file per amplicon. Part03 expects this data structure, and expects individual FASTQ file names to be structured as, e.g. `Btau_Nepal_2983_1.orth9941_1692-2530.fastq` where `Btau_Nepal_2983_1` is the individual identifier, and `orth9941_1692-2530` is the amplicon identifier. With this data structure and naming scheme, part03 can be ran through a single script:
```
./step05_make_consensus.py 
```

This script finds the most prevalent read length for each individual per amplicon and calls a degenerate consensus sequence based on all reads of that length, using the rules of [Cavener 1987 Nucleic Acids Research 15:1353–1361](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/15.4.1353) (via [Bio.motifs](http://biopython.org/DIST/docs/tutorial/Tutorial.html) in BioPython). By default, a minimum of five reads per consensus is required, any consensus reads <65 bp are removed, and if an individuals' consensus sequence length deviates >20 bp from the mean consensus sequence length for that amplicon, it is removed. These default values can be edited in `config.toml`. 

The output of part03 includes individual FASTA files for each consensus sequence, and a single multi-FASTA per amplicon containing all individuals' consensus sequences. This latter format is easily used directly to generate gene-trees, or concatenated using something like [catfasta2phyml.pl](https://github.com/nylander/catfasta2phyml) for concatenated phylogenetic analyses.

Note, the length deviation filter of part03 can throw away potentially good consensus sequences in the case that an amplicon primer set was amplifying multiple products of much different length. For example, if half of the consensus sequences are 300 bp and the other half are 100 bp, the average would be 200 bp, and >20 bp different than every consensus sequence (thus throwing away all sequences). It is a good idea to check over the summary .csv files created by part03 for each individual. If an amplicon's consensus sequences are found in these files, but are not being written into the overall summary.csv or multi-FASTA files, this may be the culprit. The 20 bp threshold in Dupuis et al. (2017) was based on a natural break in the data to remove very short sequences that were obviously not full amplicon sequences. In this scenario, a simple `cat` of all of the singel FASTA files for each individual can create an amplicon-specific mutli-FASTA.

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
