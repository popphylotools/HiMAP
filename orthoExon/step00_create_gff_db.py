#!/usr/bin/env python

import itertools
import os
import time
from multiprocessing import Pool

import gffutils


def create_db(db_path, gff_path, species):
    created = False
    os.makedirs(db_path, exist_ok=True)
    gff_name = gff_path + species + ".gff"
    db_name = db_path + species + ".gff.db"
    if not os.path.isfile(db_name):
        gffutils.create_db(gff_name,
                           dbfn=db_name,
                           force=True,
                           merge_strategy='merge',
                           id_spec=['ID', 'Name'])
        created = True
    return species, created


def connect_db(db_path, species):
    gff_name = species + ".gff"
    if os.path.isfile(db_path + gff_name + ".db"):
        db = gffutils.FeatureDB(
            db_path + gff_name + ".db")
        return species, db
    else:
        return False


def get_exons(db_path, species):
    species, db = connect_db(db_path, species)
    return species, [exon for exon in db.features_of_type('exon')]


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
    return groups


def create_and_populate_dbs(groups_fn, gff_path, db_path, template_species_list):
    start = time.clock()
    print("created?\n" +
          "--------")
    with Pool(len(template_species_list)) as p:
        results = {sp: db for sp, db in p.starmap(create_db, zip(itertools.repeat(db_path),
                                                                 itertools.repeat(gff_path),
                                                                 template_species_list))}
    for sp, status in results.items():
        print("{}: {}".format(sp, status))
    end = time.clock()
    print("time: {}".format(end - start))

    start = time.clock()
    gff_dbs = {sp: db for sp, db in [connect_db(db_path, sp) for sp in template_species_list]}
    exons = {}
    for sp in template_species_list:
        exons[sp] = gff_dbs[sp].count_features_of_type('exon')
        print("{} exons: {}".format(sp, exons[sp]))
    end = time.clock()
    print("time: {}".format(end - start))

    print(get_ortho_groups(groups_fn))


if __name__ == '__main__':
    import config

    create_and_populate_dbs(config.groups_fn, config.gff_path, config.db_path, config.template_species_list)
