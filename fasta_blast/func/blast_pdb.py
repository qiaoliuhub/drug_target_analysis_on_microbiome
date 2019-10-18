#!/usr/bin/env python

import logging
import os
from datetime import datetime
import pandas as pd
from fasta_blast.config import blast_pdb_config
from fasta_blast.utils.dwload_xtract_url import DownloadExtractURL
from fasta_blast.utils.pseq_analysis import PseqAnalysis
from apscheduler.schedulers.background import BackgroundScheduler
from apscheduler.executors.pool import ThreadPoolExecutor
import time
import random

logging.basicConfig()
logger = logging.getLogger("Blast_PDB")
logger.setLevel(logging.DEBUG)

def process_file(file, analysis_result_folder, blast_result, blast_results_folder):

    # get fasta sequences from each file, blast against pdb and save to output file in xml format
    start = datetime.now()
    dirname, filename = os.path.split(file)

    # Blasted results was saved in the .xml file under blastresult_outputfile_path directory
    basename = filename[:-len(".fsa")]
    new_file = basename + ".xml"
    blastresult_outputfile_path = os.path.join(blast_results_folder, new_file)

    threshold = 0.1 ** blast_pdb_config.threshold
    PseqAnalysis.blast_against_db(dirname, filename, blastresult_outputfile_path, threshold, blast_pdb_config.db)

    # parse the xml files and analyze them
    analysis_file = basename + ".analysis"
    parsedfile_path = os.path.join(analysis_result_folder, analysis_file)

    xml_file_dir, xml_file = os.path.split(blastresult_outputfile_path)
    aligned_records, total_records, percentage = PseqAnalysis.analyze_write_aligned_total_result(xml_file_dir, xml_file,
                                                                                                 parsedfile_path)
    blast_result.append((xml_file, aligned_records, total_records, percentage))

    end = datetime.now()
    diff = end - start
    logger.debug("Blasting and analysis takes %s seconds" % diff.seconds)

if __name__ == "__main__":

    # download and extract files
    DownloadExtractURL.download(blast_pdb_config.url, blast_pdb_config.downloaded_file_path)
    DownloadExtractURL.extract_tarfile(blast_pdb_config.downloaded_file_path, blast_pdb_config.dest_file_path)

    # get all the fasta sequences from folders
    sorted_files = sorted(DownloadExtractURL.files_end_with(blast_pdb_config.dest_file_path, ".fsa"))

    # found the parts needed to be processed
    start = int(float(blast_pdb_config.part_num - 1) / blast_pdb_config.num_of_parts * len(sorted_files))
    end = int(float(blast_pdb_config.part_num) / blast_pdb_config.num_of_parts * len(sorted_files))

    # Build the folder to save results
    ts = datetime.now()
    suffix = "_" + ts.strftime("%Y%m%d%H")
    result_folder = os.path.join(blast_pdb_config.data_dir, "blast_" + str(blast_pdb_config.db) + suffix)
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
        logger.debug("Made new folder: blast_" + str(blast_pdb_config.db) + suffix)

    blast_results_folder = os.path.join(result_folder, "blast_results")
    if not os.path.exists(blast_results_folder):
        os.makedirs(blast_results_folder)
        logger.debug("Made new folder: blast_results")

    analysis_result_folder = os.path.join(result_folder, "analysis")
    if not os.path.exists(analysis_result_folder):
        os.makedirs(analysis_result_folder)
        logger.debug("Made new folder: analysis")

    blast_result = []
    logger.debug("Setting up scheduler")
    scheduler = BackgroundScheduler(job_defaults={'misfire_grace_time': 1500*2400})
    scheduler.add_executor(ThreadPoolExecutor(10))
    random.shuffle(sorted_files) 
    try:
        for file in sorted_files[start: end]:
            scheduler.add_job(process_file, args=[file, analysis_result_folder, blast_result, blast_results_folder])

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

    try:
        logger.debug("Writing summary result out")
        summary_file = os.path.join(analysis_result_folder, "summary.analysis")
        output = pd.DataFrame(blast_result, columns=["specie_file_name", "aligned", "total", "aligned/total"])
        output.to_csv(summary_file, index=False)
        logger.debug("Writen summary result out successfully")
    except Exception as e:
        logger.debug("Fail to write summary result out, caused by %s" %str(e))
