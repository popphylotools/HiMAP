#!/usr/bin/env python

import itertools
import os
import time
import multiprocessing as mp

import gffutils


def create_db(db_path, gff_path, species):
    created = False
    os.makedirs(db_path, exist_ok=True)
    gff_name = gff_path + species + ".gff"
    db_name = db_path + species + ".gff.db"
    if not os.path.isfile(db_name):
        gffutils.create_db(data=gff_name,
                           dbfn=db_name,
                           force=True,
                           merge_strategy='merge',
                           id_spec=['ID', 'Name'])
        created = True
    return species, created


def connect_db(db_path, species):
    gff_name = species + ".gff"
    if os.path.isfile(db_path + gff_name + ".db"):
        db = gffutils.FeatureDB(db_path + gff_name + ".db")
        return species, db
    else:
        return False


def get_exons(db_path, species):
    species, db = connect_db(db_path, species)
    return species, [exon for exon in db.features_of_type('exon')]


def create_gff_db(gff_path, db_path, full_species_list):
    start = time.clock()
    print("created?\n" +
          "--------")
    with mp.Pool(min(len(full_species_list), mp.cpu_count())) as p:
        results = {sp: db for sp, db in p.starmap(create_db, zip(itertools.repeat(db_path),
                                                                 itertools.repeat(gff_path),
                                                                 full_species_list))}
    for sp, status in results.items():
        print("{}: {}".format(sp, status))
    end = time.clock()
    print("time: {}".format(end - start))

    start = time.clock()
    gff_dbs = {sp: db for sp, db in [connect_db(db_path, sp) for sp in full_species_list]}
    exons = {}
    for sp in full_species_list:
        exons[sp] = gff_dbs[sp].count_features_of_type('CDS')
        print("{} CDS's: {}".format(sp, exons[sp]))
    end = time.clock()
    print("time: {}".format(end - start))


if __name__ == '__main__':
    import argparse
    import pytoml

    parser = argparse.ArgumentParser()
    parser.add_argument('--configPath', help='configPath', default='../config.toml')
    args = parser.parse_args()

    # load config file
    with open(args.configPath) as toml_data:
        config = pytoml.load(toml_data)

    create_gff_db(config['gff_path'], config['db_path'], config['full_species_list'])
