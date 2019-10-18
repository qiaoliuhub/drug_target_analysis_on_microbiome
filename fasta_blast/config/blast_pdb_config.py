import os

# how many parts do you want to seperate the fasta files in the folder to
num_of_parts = 1

# which part of the folder to run against, starting from 1
part_num = 1

# threshold for e_value (0 < threshold < 6)
threshold = 2

# database name to use (directory path has been configured in environment variable)
db = "pdb"

# This is the url to download, this is not need to be specified if dest_file_path is not empty
url = ""

# Directory to save output data (data_dir/filename is the unzipped downloaded data directory)
data_dir = "/data/saturn/a/qliu/blast_pdb"

# filename (data_dir/filename is the unzipped downloaded data directory)
filename = "hmp_ref"

# Current working directory
cur_directory = os.getcwd()

# The file path to save the downloaded zipped fasta file, these fasta files need to be blasted
downloaded_file_path = os.path.join("/home/qliu/blast_pdb", "temp.tar.gz")

# The place to save all unzipped fasta data, these fasta files need to be blasted
dest_file_path = os.path.join(data_dir, filename)
