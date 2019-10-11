import logging
from fasta_blast.utils.dwload_xtract_url import DownloadExtractURL
from fasta_blast.config import count_aligned_total_percent_config
from fasta_blast.utils.pseq_analysis import PseqAnalysis
import os
from datetime import datetime
import pandas as pd

logging.basicConfig()
logger = logging.getLogger("Count Aligned total percentage")
logger.setLevel(logging.DEBUG)


if __name__ == "__main__":

    count_results = []

    evalue = 0.1 ** 65

    for file in DownloadExtractURL.files_end_with(count_aligned_total_percent_config.input_xml_folder, ".xml"):

        try:
            start = datetime.now()
            xml_file_dir, xml_file = os.path.split(file)
            logger.debug("Start to analyze " + str(xml_file))
            basename = xml_file[:-len(".xml")]
            parsedfile_path = os.path.join(count_aligned_total_percent_config.output_analysis_folder, basename + ".analysis")
            aligned_records, total_records, percentage = PseqAnalysis.analyze_write_aligned_total_result(xml_file_dir, xml_file, parsedfile_path, evalue)
            count_results.append((basename, aligned_records, total_records, percentage))
            end = datetime.now()
            diff = end - start
            logger.debug("Analyzed " + str(xml_file) + " successfully and it takes %s seconds" % str(diff.seconds))

        except Exception as e:
            logger.debug("Fail to analyze " + str(xml_file) + ", caused by %s" %str(e))

    try:
        logger.debug("Writing summary result out")
        summary_file = os.path.join(count_aligned_total_percent_config.output_analysis_folder, "summary.csv")
        output = pd.DataFrame(count_results, columns=["specie_file_name", "aligned", "total", "aligned/total"])
        output.to_csv(summary_file, index=False)
        logger.debug("Writen summary result out successfully")
    except Exception as e:
        logger.debug("Fail to write summary result out, caused by %s" %str(e))

