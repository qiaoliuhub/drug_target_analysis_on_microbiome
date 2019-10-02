#! /usr/bin/env python

try:
    from urllib.error import HTTPError
    from urllib.parse import quote, urlencode
    from urllib.request import urlopen
except ImportError:
    from urllib import urlencode
    from urllib2 import quote, urlopen, HTTPError

import argparse
import time
import pandas as pd
import json
import xml.etree.ElementTree as ET

# use this url to save cid list in eutils server
CIDS_LISTKEY_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/cids/JSON?list_return=listkey"

def xml2df(root):

    cols = 'ID,Name,PharmaActions,InChIKey'
    df = pd.DataFrame(columns=cols.split(","))
    for i, drug in enumerate(root):

        print(i)

        info = {}
        id = drug.find("./*[@Name='CID']")
        info['ID'] = id.text if id is not None else None
        name = drug.find("./*[@Name='MeSHHeadingList']/*[@Name='string']")
        info['Name'] = name.text if name is not None else None
        info['PharmaActions'] = ";".join(x.text for x in drug.findall("./*[@Name='PharmActionList']/*[@Name='string']"))
        inchikey = drug.find("./*[@Name='InChIKey']")
        info['InChIKey'] = inchikey.text if inchikey is not None else None
        df.loc[i] = pd.Series(info)

    return df


if __name__ == "__main__":

    argparser = argparse.ArgumentParser()
    argparser.add_argument('cid_ls_df', help = 'the data frame saving cids to be searched')
    argparser.add_argument('result_df', help = 'directory to save final results')
    args = argparser.parse_args()

    cids = pd.read_csv(args.cid_ls_df)
    cids_list = set(cids['CIDs'].str[4:].astype(int).astype(str))

    # construct apiurl and add cids into POST body
    identifiers = ','.join(str(x) for x in cids_list)
    namespace = 'cid'
    post_body = urlencode([(namespace, identifiers)]).encode('utf8')
    try:
        response = urlopen(CIDS_LISTKEY_API, post_body)
    except HTTPError as e:
        print("Fail to retrieve messages, caused by {}".format(str(e)))
        raise

    # Construct esummary retrieve url
    lsit_key_result = json.loads(response.read())
    esummary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi" \
               "?db=pccompound&WebEnv=" + lsit_key_result['IdentifierList']['EntrezWebEnv']\
               + "&query_key=" + str(lsit_key_result['IdentifierList']['EntrezQueryKey'])

    try:
        time.sleep(5)
        response2 = urlopen(esummary)
    except HTTPError as e:
        print("Fail to retrieve messages, caused by {}".format(str(e)))
        raise

    # Parsing the downloaded esummary xml string
    root = ET.fromstring(response2.read().decode('utf-8'))
    drug_df = xml2df(root)
    drug_df.to_csv(args.result_df)



