import os
import random
import zipfile

util_dir = "util_files"
graphs_dir = "graphs/"

def get_all_files_in_dir_as_list(dir):
    # get all the files in the directory
    files = os.listdir(dir)
    # get all the files in the directory as a list
    files = list(files)
    return files

def extract_file_from_zip(zip_file_path, dir_to_save, file_to_extract = 'ncbi_dataset/data/genomic.fna'):
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        # zip_ref.extractall('temp')
        zip_ref.extract(file_to_extract, dir_to_save)

def take_random_sample(x, n):
    # get a random sample of size n from list x
    return random.sample(x, n)