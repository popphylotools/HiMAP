#!/usr/bin/env python

import os
import gffutils

db_path = "./intermediate/"
db = {name.split('_')[0]: name for name in os.listdir(db_path) if ".db" in name}
print(db)
db = {key: gffutils.FeatureDB(db_path + value) for key, value in db.items()}
db_features = {name: sorted([i for i in handle.featuretypes()]) for name, handle in
               sorted(db.items())}

for name, features in db_features.items():
    print(name)
    for feature in features:
        print('\t' + feature)
    print()
