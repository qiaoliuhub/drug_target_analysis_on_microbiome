# percentage versus e-value
import logging
import os

from fasta_blast.config import aligned_protein_in_pdb_config
from fasta_blast.utils.dwload_xtract_url import DownloadExtractURL
from fasta_blast.utils.pseq_analysis import PseqAnalysis

logging.basicConfig()
logger = logging.getLogger("aligned_protein_in_pdb")
logger.setLevel(logging.DEBUG)

if __name__ == "__main__":

    # each element in results = [evalue, total_seq, aligned_seq, percentage]
    distri_results = [[i, 0, 0, 0] for i in range(aligned_protein_in_pdb_config.maximum_evalue_index + 1)]
    aligned_seqs = [0] * (aligned_protein_in_pdb_config.maximum_evalue_index + 1)
    total_seqs = 0

    for file in DownloadExtractURL.files_end_with(aligned_protein_in_pdb_config.input_xml_folder, ".xml"):

        # process each file and update results
        dirname, filename = os.path.split(file)
        logger.debug("Finding aligned sequence")
        total, aligned_seqs = PseqAnalysis.aligned_seq_num_distribution_of_evalue(dirname, filename, aligned_seqs)
        total_seqs += total
        logger.debug("Found aligned sequence with evalues successfully")
        logger.debug("Total_seqs: %s, Aligned_seqs: %s" % (total_seqs, aligned_seqs))


    for i in range(len(distri_results)):
        distri_results[i][1] = total_seqs
        distri_results[i][2] = aligned_seqs[i]
        distri_results[i][3] = float(distri_results[i][2]) / (distri_results[i][1] + 0.1)

    logger.debug("Total_seqs, Aligned_seqs, Percentage: %s" % str(distri_results))

    distri_file_path = os.path.join(aligned_protein_in_pdb_config.result_folder, "distribution_summary.distri")
    with open(distri_file_path, "w") as output:
        output.write(str(distri_results))
