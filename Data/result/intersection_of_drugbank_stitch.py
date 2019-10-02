#! /usr/bin/env python

import pandas as pd
import sys

if __name__ == '__main__':
    db_label = pd.read_csv(sys.argv[1])
    cids = pd.read_csv(sys.argv[2])
    set_stitch = set(cids['CIDs'].str[4:].astype(int))
    set_drugbank = set(db_label[~db_label['PubChem Compound ID'].isnull()]['PubChem Compound ID'].astype(int))
    common = set_stitch & set_drugbank
    pd.DataFrame({'cids': list(common)}).to_csv(sys.argv[3], index = False)
