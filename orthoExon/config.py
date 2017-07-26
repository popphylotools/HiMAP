# species with good annotations
template_species_list = ["Bcur", "Bdor", "Bole", "Ccap"]

# full species list
enhanced_species_list = ['Bjar', 'Aobl', 'Bmin', 'Asus', 'Btry', 'Afra', 'Blat', 'Bzon', 'Bcor', 'Ccap', 'Bcur', 'Bole',
                         'Bdor']

# order in which species should appear in results
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

# set of species who's data came from transvestigator (needs to be parsed slightly different)
transvestigated_species_set = {'Bcor', 'Blat', 'Bzon', 'Afra', 'Bmin', 'Bjar', 'Aobl'}

# number of n's used to create alignment boarders
n_count = 50

# cds gap filter parameters
max_gap_percent = 0
max_gap_length = 0

# cds length filter parameters
min_cds_length = 100

# path to the groups file
groups_fn = "../data/input/groups_filtered_6181genes.txt"

# path to gff folder
gff_path = "../data/input/gff/"

# path to fasta folder
fasta_path = "../data/input/fasta/"

# path where databases are created
db_path = "../data/intermediate/gff_databases/"

# path where json files are created
json_path = "../data/intermediate/json/"

# path where alignment with template species is created
template_alignment_path = "../data/intermediate/template_alignment/"

# path where alignment with all species is created
enhanced_alignment_path = "../data/intermediate/supplemental_alignment/"

#
fasta_output_path = "../data/output/fasta/"

#
summary_fn = "../data/output/summory.csv"

#
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

#
expand_iupac = {value.upper(): set(key) for key, value in collapse_iupac.items()}
