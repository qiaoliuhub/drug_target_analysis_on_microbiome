import os
import tarfile
import logging
import urllib

logging.basicConfig()
logger = logging.getLogger("DownloadExtractURL")
logger.setLevel(logging.DEBUG)

class DownloadExtractURL:

    @staticmethod
    def __reporthook(count, blockSize, totalSize):

        total_blocks_count = totalSize / blockSize
        base_count = total_blocks_count / 20
        percent = count * (float(blockSize) / totalSize) * 100
        if count % base_count == 0:
            logger.debug("Downloaded %s percent", str(percent))

    @staticmethod
    def download(url, download_file_path):

        # download the zipped file to file_path
        if not os.path.exists(download_file_path):
            # when the file was not downloaded yet
            try:
                logger.debug("Downloading file to %s", download_file_path)
                urllib.urlretrieve(url, download_file_path, reporthook = DownloadExtractURL.__reporthook)
                logger.debug("Downloaded file to %s successfully", download_file_path)
            except Exception as e:
                logger.error("Fail to download files, caused by %s", str(e))
        else:
            logger.debug("File %s exists", download_file_path)

    @staticmethod
    def __fasta_files(members):
        for tarinfo in members:
            if os.path.splitext(tarinfo.name)[1] == ".fsa":
                yield tarinfo

    @staticmethod
    def extract_tarfile(source_file_path, dest_file_path):

        if not os.path.exists(dest_file_path):
            # if source file has not been extracted before
            try:
                logger.debug("Extracting files")
                tar = tarfile.open(source_file_path)
                tar.extractall(dest_file_path)
                logger.debug("Extracted files successfully")
            except Exception as e:
                logger.error("Fail to extract files, caused by %s", str(e))
            finally:
                tar.close()

        else:
            logger.debug("Extracted file %s exists", dest_file_path)

    @staticmethod
    def files_end_with(folder, pattern):
        results = []
        files = [x for x in os.walk(folder)]
        for file_list in files:
            for file in file_list[2]:
                if file.endswith(pattern):
                    candidate = os.path.join(file_list[0], file)
                    results.append(candidate)
        logger.debug("Found %d fasta files", len(results))

        return results
