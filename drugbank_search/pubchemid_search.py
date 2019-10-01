import sys
import pandas as pd
import xml.etree.ElementTree as ET
from drugbank_search.drugbankid_search import pubchemid2infos
if __name__ == "__main__":

    with open(sys.argv[1]) as drugs_file:
        drugs = drugs_file.readline().split(',')

    drugs = list(map(str.strip, drugs))
    parse_tree = ET.parse(sys.argv[2])
    root = parse_tree.getroot()

    cols = 'DrugBank ID,Name,CAS Number,Drug Groups,InChIKey,InChI,SMILES,Formula,KEGG Compound ID,KEGG Drug ID,PubChem Compound ID,PubChem Substance ID,ChEBI ID,ChEMBL ID,HET ID,ChemSpider ID,BindingDB ID'
    drugs_df = pd.DataFrame(columns=cols.split(","))
    for i, drug in enumerate(drugs):
        drugs_df.loc[i] = pubchemid2infos(root, drug)

    unfound = drugs_df[drugs_df['Name'].isnull()]['DrugBank ID']
    with open(sys.argv[4], 'w+') as unfound_file:
        unfound_file.write(','.join(list(unfound)))
    drugs_df.to_csv(sys.argv[3], index=False)