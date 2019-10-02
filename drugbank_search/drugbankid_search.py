#!/usr/bin/env python

import xml.etree.ElementTree as ET
import sys
import pandas as pd


# def dbID2property(root, dbID, db_property):
#     ns = {'db_ns': 'http://www.drugbank.ca'}
#     for drug in root.findall(".//*[db_ns:drugbank-id='" + dbID + "']/*[@primary='true']/..", ns):
#
#         for drug_property in drug.findall(".//*[db_ns:kind='" + db_property + "']", ns):
#             return drug_property.find('db_ns:value', ns).text


def dbID2property(db_property, drug):
    ns = {'db_ns': 'http://www.drugbank.ca'}
    info = drug.find("./*/db_ns:property[db_ns:kind='" + db_property + "']", ns)
    if info is not None:
        return info.find('db_ns:value', ns).text
    else:
        return None


def dbID2externalid(resource, drug):
    ns = {'db_ns': 'http://www.drugbank.ca'}
    eid = drug.find("./*/db_ns:external-identifier[db_ns:resource='" + resource + "']", ns)
    if eid is not None:
        return eid.find('db_ns:identifier', ns).text
    else:
        return None


def dbID2infos(root, dbID):
    drug_info = {'DrugBank ID': dbID}
    ns = {'db_ns': 'http://www.drugbank.ca'}
    print(dbID)
    for drug in root.findall(".//*[db_ns:drugbank-id='" + dbID + "']/*[@primary='true']/..", ns):
        name = drug.find('./db_ns:name', ns)
        drug_info['Name'] = name.text if name is not None else None
        cas = drug.find('./db_ns:cas-number', ns)
        drug_info['CAS Number'] = cas.text if cas is not None else None
        drug_info['Drug Groups'] = '; '.join([x.text for x in drug.findall('./db_ns:groups/db_ns:group', ns)])
        drug_info['InChIKey'] = dbID2property('InChIKey', drug)
        drug_info['InChI'] = dbID2property('InChI', drug)
        drug_info['SMILES'] = dbID2property('SMILES', drug)
        drug_info['Formula'] = dbID2property('Molecular Formula', drug)
        drug_info['KEGG Compound ID'] = dbID2externalid('KEGG Compound', drug)
        drug_info['KEGG Drug ID'] = dbID2externalid('KEGG Drug', drug)
        drug_info['PubChem Compound ID'] = dbID2externalid('PubChem Compound', drug)
        drug_info['PubChem Substance ID'] = dbID2externalid('PubChem Substance', drug)
        drug_info['ChEBI ID'] = dbID2externalid('ChEBI', drug)
        drug_info['ChEMBL ID'] = dbID2externalid('ChEMBL', drug)
        drug_info['HET ID'] = dbID2externalid('PDB', drug)
        drug_info['ChemSpider ID'] = dbID2externalid('ChemSpider', drug)
        drug_info['BindingDB ID'] = dbID2externalid('BindingDB', drug)

    return pd.Series(drug_info)


def pubchemid2infos(root, pubchemid):
    
    drug_info = {'PubChem Compound ID': pubchemid}
    ns = {'db_ns': 'http://www.drugbank.ca'}
    print(pubchemid)
    for drug in root.findall(".//*[db_ns:resource='PubChem Compound'][db_ns:identifier='" + pubchemid + "']/../..", ns):
        name = drug.find('./db_ns:name', ns)
        drug_info['Name'] = name.text if name is not None else None
        cas = drug.find('./db_ns:cas-number', ns)
        drug_info['CAS Number'] = cas.text if cas is not None else None
        drug_info['Drug Groups'] = '; '.join([x.text for x in drug.findall('./db_ns:groups/db_ns:group', ns)])
        drug_info['InChIKey'] = dbID2property('InChIKey', drug)
        drug_info['InChI'] = dbID2property('InChI', drug)
        drug_info['SMILES'] = dbID2property('SMILES', drug)
        drug_info['Formula'] = dbID2property('Molecular Formula', drug)
        drug_info['KEGG Compound ID'] = dbID2externalid('KEGG Compound', drug)
        drug_info['KEGG Drug ID'] = dbID2externalid('KEGG Drug', drug)
        drug_info['DrugBank ID'] = drug.find("./db_ns:drugbank-id[@primary='true']", ns)
        drug_info['PubChem Substance ID'] = dbID2externalid('PubChem Substance', drug)
        drug_info['ChEBI ID'] = dbID2externalid('ChEBI', drug)
        drug_info['ChEMBL ID'] = dbID2externalid('ChEMBL', drug)
        drug_info['HET ID'] = dbID2externalid('PDB', drug)
        drug_info['ChemSpider ID'] = dbID2externalid('ChemSpider', drug)
        drug_info['BindingDB ID'] = dbID2externalid('BindingDB', drug)

    return pd.Series(drug_info)

def get_number_of_drugs(df, col):

    return len(df[~df[col].isnull()])


if __name__ == "__main__":

    with open(sys.argv[1]) as drugs_file:
        drugs = drugs_file.readline().split(',')

    drugs = list(map(str.strip, drugs))
    parse_tree = ET.parse(sys.argv[2])
    root = parse_tree.getroot()

    cols = 'DrugBank ID,Name,CAS Number,Drug Groups,InChIKey,InChI,SMILES,Formula,KEGG Compound ID,KEGG Drug ID,PubChem Compound ID,PubChem Substance ID,ChEBI ID,ChEMBL ID,HET ID,ChemSpider ID,BindingDB ID'
    drugs_df = pd.DataFrame(columns=cols.split(","))
    for i, drug in enumerate(drugs):
        drugs_df.loc[i] = dbID2infos(root, drug)

    unfound = drugs_df[drugs_df['Name'].isnull()]['DrugBank ID']
    with open(sys.argv[4], 'w+') as unfound_file:
        unfound_file.write(','.join(list(unfound)))
    drugs_df.to_csv(sys.argv[3], index=False)
