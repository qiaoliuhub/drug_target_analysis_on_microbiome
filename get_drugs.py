#!/usr/bin/env python

import pandas as pd
from shutil import copyfile
import os
import re
import sys

def process_file(old_file):
    with open(old_file) as microbe_file:
        first_line = microbe_file.readline()
        if not len(first_line):
            return
        folder_1, filename_1 = os.path.split(old_file)
        microbe = " ".join(first_line.split("[")[1].split()[:2])
    print(old_file)
    blast_file_part = ".".join(old_file.split("/")[-1].split(".")[:2])
    folder, filename = os.path.split(os.path.split(old_file)[0])
    src = os.path.join(folder, "blast_results", blast_file_part + ".pep.xml")
    if not os.path.exists(src):
        return
    pattern = 'DB\d+'

    drugs_list = []
    
    with open(src) as src_file:
        for line in src_file:
            result = re.findall(pattern, line)
            drugs_list += result

    if not os.path.exists(os.path.join(folder,"drugs_2")):
        os.mkdir(os.path.join(folder,"drugs_2"))

    target = os.path.join(folder,"drugs_2",microbe)

    if os.path.exists(target):	
        with open(target, 'a+') as target_file:
            target_file.write(",\n")

    with open(target, 'a+') as target_file:
        target_file.write(",\n".join(list(set(drugs_list))))

if __name__ == '__main__':
    
    old_files = []
    for root, dirs, files in os.walk(sys.argv[1]):
        for old_file in files:
            if old_file.endswith(".fsa"):
                cur_file = os.path.join(root, old_file)
                old_files.append(cur_file)
                
    print("get all files")
    print(str(old_files))
    for old_file in old_files:
        process_file(old_file)
