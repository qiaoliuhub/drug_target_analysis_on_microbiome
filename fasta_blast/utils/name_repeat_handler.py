import shutil
import logging

logging.basicConfig()
logger = logging.getLogger("NmpPepHandler")
logger.setLevel(logging.DEBUG)

class NameRepeatHandler:

    input_output_map = {}
    # {"first few number in input_file name" => "processed file path"}

    @classmethod
    def preprocess(cls, input_file, output_file_path, tag = ""):
        # check if input_file has been processed before
        ids = input_file.split(".")[0] + tag
        if ids in cls.input_output_map:

            logger.debug("%s exists" %ids)
            processed_file_path = cls.input_output_map[ids]
            with open(processed_file_path, "r") as src_file, open(output_file_path, "w") as dest_file:
                shutil.copyfileobj(src_file, dest_file)

            logger.debug("%s was copied to %s" %(processed_file_path, output_file_path))

        else:
            logger.debug("No %s was found" %ids)

    @classmethod
    def add_to_map(cls, input_file, output_file_path, tag = ""):

        # add new relationship
        ids = input_file.split(".")[0] + tag
        if ids not in cls.input_output_map:

            cls.input_output_map[ids] = output_file_path
            logger.debug("%s was added to map" % ids)
