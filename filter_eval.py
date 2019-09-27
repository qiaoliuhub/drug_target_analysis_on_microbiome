import pandas as pd
import os
import re
from apscheduler.schedulers.background import BackgroundScheduler, BlockingScheduler
from apscheduler.executors.pool import ThreadPoolExecutor
import argparse
import logging
import config
import xml.etree.ElementTree as ET
import fcntl
import time

formatter = logging.Formatter(fmt='%(asctime)s %(levelname)s %(name)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S')
filehandler = logging.FileHandler(config.logfile, mode='w+')
filehandler.setFormatter(fmt=formatter)
logger = logging.getLogger("Filter Homologs")
logger.addHandler(filehandler)
logger.setLevel(logging.DEBUG)
sche_logger = logging.getLogger('apscheduler.executors.default')
sche_logger.addHandler(filehandler)
sche_logger.setLevel(logging.INFO)

def process_file(old_file):

    ### Read in all the xml files and only select the balst results with evalue <10e-60
    ### build file with name as microbiome in drug_targets folder
    ### format: drugbank ID, Iteration_query-def, hit-def, Hsp_evalue (Hsp_num = 1), most_significant
    ### drugbankID: keep all drugbank ID for one queries
    ### Iteration_query_def: query proteins (one per iteration)
    ### Hit_def: zero or many Hits per iteration query proteins (sorted), many high score parts(hsp) per Hit (sorted)
    ### Hit_def: the name of Hits (one per hit) (keep all qualitfied hits for query but flag the most significant one)
    ### BlastOutput-->BlastOutput_iterations-->Iteration-->Iteration_hits-->Hit-->Hit_hsps-->Hsp-->Hsp_evalue
    ### Iteration_query_def (one->many) hit-def (one->one) Hsp_evalue (one->many) drugbank ID
    with open(old_file) as microbe_file:
        first_line = microbe_file.readline()
        if not len(first_line):
            return
        microbe = " ".join(first_line.split("[")[1].split()[:2])

    logger.debug("Processing file %s" % old_file)
    blast_file_part = ".".join(old_file.split("/")[-1].split(".")[:2])
    folder, filename = os.path.split(os.path.split(old_file)[0])
    src = os.path.join(folder, "blast_results", blast_file_part + ".pep.xml")
    if not os.path.exists(src):
        return

    parse_tree = ET.parse(src)
    root = parse_tree.getroot()
    drug_targets = pd.DataFrame(columns=['query_proteins', 'target_uniprot', 'target_proteins', 'drugs', 'e_value'])

    for iteration in root.findall('.//Iteration'):

        Iteration_query_def = iteration.find('Iteration_query-def').text

        for hit in iteration.findall('Iteration_hits/Hit'):

            hit_dict = {'query_proteins': Iteration_query_def}

            hsp_evalue = hit.find(".//*[Hsp_num='1']").find('Hsp_evalue').text
            if float(hsp_evalue) > 1e-60:
                continue
            hit_dict['e_value'] = hsp_evalue

            hit_def = hit.find('Hit_def').text.split("|")[1]
            pattern = 'DB\d+'
            db_drugs = re.findall(pattern=pattern, string=hit_def)
            hit_dict['drugs'] = ",".join(db_drugs)

            hit_uniprot, hit_proteins = hit_def.split('(DB')[0].split(maxsplit=1)
            hit_dict['target_uniprot'] = hit_uniprot
            hit_dict['target_proteins'] = hit_proteins

            i = len(drug_targets)
            drug_targets.loc[i, :] = pd.Series(hit_dict)

    microbiome_file = os.path.join(folder, "drug_targets", microbe)
    print_header = False if os.path.exists(microbiome_file) else True
    logger.debug("acquiring lock")
    cur_file = open(microbiome_file, 'a+')

    start_time = time.time()
    while True:
        try:
            if time.time() - start_time > 3600:
                break
            fcntl.flock(cur_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
            break
        except BlockingIOError:
            time.sleep(1)

    logger.debug("acquired lock")
    drug_targets.to_csv(microbiome_file, index=False, header=print_header, mode='a+')
    fcntl.flock(cur_file, fcntl.LOCK_UN)
    logger.debug("Lock is released")

if __name__ == '__main__':

    # Parse command line arguments
    argparser = argparse.ArgumentParser()
    argparser.add_argument('directory', help='The fsa and xml file in this directory')
    args = argparser.parse_args()
    input = args.directory

    old_files = []
    for root, dirs, files in os.walk(input):
        for old_file in files:

            if old_file.endswith(".fsa"):
                cur_file = os.path.join(root, old_file)
                old_files.append(cur_file)
                folder, filename = os.path.split(os.path.split(cur_file)[0])
                if not os.path.exists(os.path.join(folder, "drug_targets")):
                    os.mkdir(os.path.join(folder, "drug_targets"))

    print("get all files")
    print(str(old_files))

    logger.debug("Setting up scheduler")
    scheduler = BackgroundScheduler(job_defaults={'misfire_grace_time': 1500*60})
    scheduler.add_executor(ThreadPoolExecutor(10))
    try:

        for old_file in old_files:
            scheduler.add_job(process_file, args = [old_file])

        scheduler.start()

        while len(scheduler.get_jobs()):
            time.sleep(1)
            continue
    except:
        raise

    finally:
        scheduler.print_jobs()
        scheduler.shutdown(wait=True)
        scheduler.print_jobs()
