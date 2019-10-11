import logging
import math
import os

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline

from fasta_blast.utils.name_repeat_handler import NameRepeatHandler

logging.basicConfig()
logger = logging.getLogger("PseqAnalysis")
logger.setLevel(logging.DEBUG)

class PseqAnalysis:

    @staticmethod
    def blast_against_db(source_dirname, source_filename, outputfile_path, threshold, db):

        # Prevent repetitive processing
        if os.path.exists(outputfile_path):
            logger.debug("Blasted results file %s exists", outputfile_path)
            return

        source_file = os.path.join(source_dirname, source_filename)
        NameRepeatHandler.preprocess(source_filename, outputfile_path, "_blast")
        logger.debug("Blasting against "+ str(db) +" database for file %s", source_file)

        if not os.path.exists(outputfile_path):
            # build the command
            blastp_cline = NcbiblastpCommandline(query=source_file, db=db, evalue=threshold, outfmt=5, out=outputfile_path)
            blastp_cline()

        # Prevent repetitive processing
        NameRepeatHandler.add_to_map(source_filename, outputfile_path, "_blast")
        logger.debug("Done with blasting")

    @staticmethod
    def get_total_aligned_in_blast_records(blast_records, evalue):

        # Analyze a blast records to get the total number of seqs and the number of aligned seqs

        total_records, aligned_records = 0, 0
        for blast_record in blast_records:
            total_records += 1
            aligned = False
            for alignment in blast_record.alignments:

                for hsp in alignment.hsps:
                    if hsp.expect < evalue:
                        aligned = True
                        break

                if aligned:
                    break

            if aligned:
                aligned_records += 1

        return total_records, aligned_records

    @staticmethod
    def analyze_write_aligned_total_result(source_dirname, source_filename, parsedfile_path, evalue = 0.1 ** 3):

        if os.path.exists(parsedfile_path):
            logger.debug("Parsed file %s exists", parsedfile_path)
            return

        # Prevent repetitive processing
        NameRepeatHandler.preprocess(source_filename, parsedfile_path, "_analysis")

        # get aligned and total sequences number results in the file
        xml_file = os.path.join(source_dirname, source_filename)
        logger.debug("Analyzing xml files %s", xml_file)
        try:
            result_handle = open(xml_file)
            blast_records = NCBIXML.parse(result_handle)
            total_records, aligned_records = 0,0
            total_records, aligned_records = PseqAnalysis.get_total_aligned_in_blast_records(blast_records, evalue)

        except Exception as e:
            logger.error("Fail to analyze file %s, caused by %s", xml_file, str(e))
            error_file = os.path.join(os.getcwd(), "error_file")
            with open(error_file, "w") as err:
                err.write("Fail to analyze file %s" %xml_file)

        finally:
            result_handle.close()

        # Writing the result to one file
        logger.debug("Writing results to %s", parsedfile_path)
        with open(parsedfile_path, 'w+') as output:
            if total_records:
                output.write("%s, %d, %d, %f" %(source_filename, aligned_records, total_records, float(aligned_records)/(total_records+0.1)))
        # Prevent repetitive processing
        NameRepeatHandler.add_to_map(source_filename, parsedfile_path, "_analysis")
        return aligned_records, total_records, float(aligned_records)/(total_records+0.1)

    @staticmethod
    def __generate_partial_aligned_records(blast_records, global_min_evalue, partial_aligned_records):

        # generate a list of aligned seq nums distribution of evalue partially (check
        # aligned_seq_num_distribution_of_evalue for full explanation

        total_records = 0
        for blast_record in blast_records:

            # for each sequence, there are multiple alignments
            total_records += 1
            # find the lowest evalue for this sequence
            cur_min_evalue = 1

            for alignment in blast_record.alignments:

                # logger.debug("alignment: %s" %alignment.title)
                for hsp in alignment.hsps:
                    # logger.debug("evalue: %s" %hsp.expect)
                    if hsp.expect < global_min_evalue:
                        cur_min_evalue = global_min_evalue
                        break
                    if hsp.expect < cur_min_evalue:
                        cur_min_evalue = hsp.expect
                if cur_min_evalue <= global_min_evalue:
                    break

            # logger.debug("min_evalue is %s" % cur_min_evalue)
            max_threshold = min(int(math.floor(-math.log10(cur_min_evalue))), len(partial_aligned_records) - 1)
            # logger.debug("max_threshold is %s" % max_threshold)

            # add 1 to partial_aligned_record indexed with max_threshold
            partial_aligned_records[max_threshold] += 1

        return total_records


    @staticmethod
    def aligned_seq_num_distribution_of_evalue(dirname, filename, results):

        # process a file and add total_seqs and aligned_seqs to the results

        xml_file = os.path.join(dirname, filename)
        logger.debug("Analyzing xml files %s", xml_file)
        try:
            result_handle = open(xml_file)
            blast_records = NCBIXML.parse(result_handle)
            logger.debug("Successfully parsed the xml file %s", xml_file)

            # each entry i in partial_aligned_results represents the number of sequence whose evalue is just exactly less than 0.1 ** i,
            # for example, if a best sequence alignment evalue is 0.002, then add 1 to partial_aligned_results[2], not to others
            partial_aligned_records = [0] * len(results)
            global_min_evalue = 0.1 ** len(results)
            total_records = PseqAnalysis.__generate_partial_aligned_records(blast_records, global_min_evalue, partial_aligned_records)

            # total_records denotes the total sequences in this file, aligned_results denotes the number of sequence which has
            # alignment better than that threshold
            aligned_results = PseqAnalysis.__wrap_up_partial_results(partial_aligned_records)

            # update results
            for i in range(len(results)):
                results[i] += aligned_results[i]

            logger.debug("total_records: %s, results: %s" %(total_records, results))
            return total_records, results

        except Exception as e:
            logger.error("Fail to analyze file %s, caused by %s", file, str(e))
            error_file = os.path.join(os.getcwd(), "error_file")
            with open(error_file, "w") as err:
                err.write("Fail to analyze file %s" % file)

            return 0, results

        finally:
            result_handle.close()

    @staticmethod
    def __wrap_up_partial_results(partial_results):

        results = [0] * len(partial_results)
        cur_sum = 0

        for i in range(len(partial_results) - 1, -1, -1):
            cur_sum += partial_results[i]
            results[i] = cur_sum

        return results
