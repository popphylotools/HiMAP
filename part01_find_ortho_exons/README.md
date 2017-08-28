#### step00a
This script reads in the .gff files and produces sqlite databases for use with gffutils

#### step00b
This script uses the ortho -> sp -> accession relationships from the filtered orthomcl output
to create a json file with ortho -> sp -> cds_parent_database_id relationships which are used later.

#### step01
This script takes the "high-quality data sources" (well-annotated genomes) and creates template fastas of gene orthologs
consisting of CDS's concatenated with fixed length strings of 'N's in place of introns.
It then runs MAFFT on these fasta's to align them.

#### step02
This script takes the aligned padded orthologs from the privious step and splits them on the exon boundries
where conserved. It filters them based on parameters in the config file, then outputs fastas with the template
orthologus exons and the gene sequences from the remaining species. These fasta's are then aligned with mafft.

#### step03
This script takes the alignments from the privious step and splits them on the exon boundries of the template
species if conserved. The resulting orthologus exons are then filtered based on parameters in the config file
and output to fasta files.
