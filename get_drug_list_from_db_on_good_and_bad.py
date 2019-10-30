#!/usr/bin/env python


import pandas as pd
import os
import sys

beneficial_path = 'Data/new_ori_beneficial.csv'
pathogen_path = 'Data/new_ori_pathogen.csv'

beneficial = pd.read_csv(beneficial_path)
pathogen = pd.read_csv(pathogen_path)

beneficial_names = set(beneficial['name'])
pathogen_names = set(pathogen['name'])

input_prefix = '/data/saturn/a/qliu/blastdrugbank'
input_postfix = 'local/db/repository/ncbi/dacc_reference_genomes/20141006/drugs_2'
cur_folders = ['hmp_ref_' + str(i) for i in range(8)]

target_bad_drugs = set()
target_good_drugs = set()

for cur_folder in cur_folders: 
    micro_folders = os.path.join(input_prefix, cur_folder, input_postfix)

    for micro_file in os.listdir(micro_folders):
        if micro_file in beneficial_names:
            drugfile = os.path.join(micro_folders, micro_file)
            drugs = set(pd.read_csv(drugfile, header = None)[0])
            target_good_drugs = target_good_drugs.union(drugs)
        elif micro_file in pathogen_names:
            drugfile = os.path.join(micro_folders, micro_file)
            drugs = set(pd.read_csv(drugfile, header = None)[0])
            target_bad_drugs = target_bad_drugs.union(drugs)

sel_drugs = target_bad_drugs - target_good_drugs
output_file = sys.argv[1]
target_good_drug_file = sys.argv[2]
target_bad_drug_file = sys.argv[3]

with open(output_file, 'w+') as output:
    output.write(",".join(list(sel_drugs)))

with open(target_good_drug_file, 'w+') as target_good_drug_output:
    target_good_drug_output.write(",".join(list(target_good_drugs)))

with open(target_bad_drug_file, 'w+') as target_bad_drug_output:
    target_bad_drug_output.write(",".join(list(target_bad_drugs)))

