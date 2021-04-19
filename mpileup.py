import pandas as pd
import numpy as np
import os
import re
import shutil
import subprocess
import pysam
import csv
import snptools as snp


def samtools_mpileup(path_to_files, file_extension, path_to_reference_file):

    # for files inside individual mouse folder
    # have list of folders containing prefixed files
    # Logic to do it for all folders
    # folders=['mouse_1']
    # for folder in folders:
    parent_folder = snp.get_file_basename(path_to_files)
    print(parent_folder)
    print(path_to_files)
    command = f'samtools mpileup -uf {path_to_reference_file} {path_to_files}/*{file_extension} > {parent_folder}.raw.bcf' 
    print('Executing:',command)
    subprocess.call([command],shell=True)


def main():
    samtools_mpileup('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/mouse_1','pfx.sorted.bam','/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/Combined.fasta')

if __name__ == '__main__':
    main()

