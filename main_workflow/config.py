species_list = ["Bcur", "Bdor", "Bole", "Ccap"]
full_species_list = ['Bjar', 'Aobl', 'Bmin', 'Asus', 'Btry', 'Afra', 'Blat', 'Bzon', 'Bcor',
                     'Ccap', 'Bcur', 'Bole', 'Bdor']
transvestigated_species_set = {'Bcor', 'Blat', 'Bzon', 'Afra', 'Bmin', 'Bjar', 'Aobl'}

groups_fn = "../input/groups_filtered_6181genes.txt"
fasta_path = "../input/fasta/"

db_path = "../intermediate/gff_databases/"
json_path = "../intermediate/json/"

subset_alignment_path = "../intermediate/4spp_alignment/"
fullset_alignment_path = "../intermediate/13spp_alignment/"

primer3_path = "../intermediate/primer_design/"

padded_primer_product_path = "../intermediate/phylo_informativeness/fasta/"
pi_score_path = "../intermediate/phylo_informativeness/tapir_out/phylogenetic-informativeness.sqlite"

orthoCds_path = "../output/orthoCds/"
unpadded_primer_product_path = "../output/primerProducts/"
summary_fn = "../output/summory.csv"
