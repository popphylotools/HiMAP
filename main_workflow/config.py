template_species_list = ["Bcur", "Bdor", "Bole", "Ccap"]
full_species_list = ['Bjar', 'Aobl', 'Bmin', 'Asus', 'Btry', 'Afra', 'Blat', 'Bzon', 'Bcor',
                     'Ccap', 'Bcur', 'Bole', 'Bdor']
transvestigated_species_set = {'Bcor', 'Blat', 'Bzon', 'Afra', 'Bmin', 'Bjar', 'Aobl'}

n_count = 50  # number of n's used to create alignment boarders

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

collapse_iupac = {
    ('-',): '-',
    ('A',): 'A',
    ('G',): 'G',
    ('C',): 'C',
    ('T',): 'T',
    ('C', 'T'): 'y',
    ('A', 'G'): 'r',
    ('A', 'T'): 'w',
    ('C', 'G'): 's',
    ('G', 'T'): 'k',
    ('A', 'C'): 'm',
    ('A', 'G', 'T'): 'd',
    ('A', 'C', 'G'): 'v',
    ('A', 'C', 'T'): 'h',
    ('C', 'G', 'T'): 'b',
    ('A', 'C', 'G', 'T'): 'n',
}

expand_iupac = {value.upper(): set(key) for key, value in collapse_iupac.items()}

sp_order = {'Bcur': 1,
            'Bdor': 2,
            'Bole': 3,
            'Ccap': 4,
            'Bcor': 5,
            'Blat': 6,
            'Bzon': 7,
            'Afra': 8,
            'Bmin': 9,
            'Bjar': 10,
            'Aobl': 11,
            'Asus': 12,
            'Btry': 13}
