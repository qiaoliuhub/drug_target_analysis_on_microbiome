#! /usr/bin/env python

import pandas as pd
from drugbank_search.drugbankid_search import get_number_of_drugs
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("input_df", help="search results with pubchemid on drugbank database", type=str)
    parser.add_argument('column_name', help="column name to be count")
    args = parser.parse_args()
    df = pd.read_csv(args.input_df)

    new_df = df.drop_duplicates()

    print("Out of {} compounds, {} drugs are found".format(len(new_df), get_number_of_drugs(new_df, args.column_name)))