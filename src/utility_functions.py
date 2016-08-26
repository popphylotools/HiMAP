import os
import gffutils


def create_db(db_path, gff_path, species):
    gff_name = species + ".gff"
    created = False
    if not os.path.isfile(
            os.path.join(os.path.abspath(os.path.curdir), db_path + gff_name + ".db")):
        fn = gffutils.example_filename(
            os.path.join(os.path.abspath(os.path.curdir), gff_path + gff_name))
        db = gffutils.create_db(fn,
                                dbfn=db_path + gff_name + ".db",
                                force=True,
                                merge_strategy='merge',
                                id_spec=['ID', 'Name'])
        created = True
    return (species, created)


def connect_db(db_path, species):
    gff_name = species + ".gff"
    if os.path.isfile(os.path.join(os.path.abspath(os.path.curdir), db_path + gff_name + ".db")):
        db = gffutils.FeatureDB(
            os.path.join(os.path.abspath(os.path.curdir), db_path + gff_name + ".db"))
    return (species, db)


def get_exons(db_path, species):
    species, db = connect_db(db_path, species)
    return (species, [exon for exon in db.features_of_type('exon')])


def get_ortho_groups(ortho_group_path):
    # import ortholog groups
    with open(ortho_group_path, 'r') as f:
        groups_raw = f.readlines()
        groups_raw = [line.strip() for line in groups_raw]
    groups = dict()
    for line in groups_raw:
        ortho, data = line.split(':')
        ortho = ortho.strip()
        data = data.strip().split()
        data = {elem.split("|")[0]: elem.split("|")[1] for elem in data}
        groups[ortho] = data
    return (groups)
