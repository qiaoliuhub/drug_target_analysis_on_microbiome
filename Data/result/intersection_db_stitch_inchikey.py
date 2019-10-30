#! /usr/bin/env python

import pandas as pd
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("drugbank_df", help = 'drugbank screened drugs information')
    parser.add_argument("pubchem_df", help = 'stitch screened drugs information')
    parser.add_argument("output_file", help='file to save filtered drugs')
    parser.add_argument("--header", action='store_true', help="if set to false, don't output file header")
    args = parser.parse_args()
    db_info = pd.read_csv(args.drugbank_df).drop_duplicates()
    pc_info = pd.read_csv(args.pubchem_df).drop_duplicates()
    db_info['join_inchikey'] = db_info['InChIKey'].str.split("-").str[0]
    pc_info['join_inchikey'] = pc_info['InChIKey'].str.split("-").str[0]
    pc_info = pc_info.drop_duplicates(subset=['join_inchikey'])
    db_info = db_info[~db_info['join_inchikey'].isnull()]
    pc_info = pc_info[~pc_info['join_inchikey'].isnull()]
    new_df = db_info.merge(pc_info, on = 'join_inchikey')
    new_df.to_csv(args.output_file, header = args.header, index = False, mode='a+')

    print("Out of {} compounds, {} drugs are found".format(len(db_info), len(new_df)))
