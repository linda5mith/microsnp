import pandas as pd
import numpy as np
import os
import re
import shutil
import subprocess
import pysam
import csv
import snptools as snp


#samtools mpileup /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/Combined.fasta


def samtools_mpileup(path_to_files, file_extension, path_to_reference_file):
    '''Run bcftools mpileup on the bam contents of multiple folders. Path to faidx indexed reference sequence file path_to_reference_file 
    needs to be included.'''
    for folder in os.listdir(path_to_files):
        path_to_folder=os.path.join(path_to_files,folder)
        if os.path.isdir(path_to_folder):
            output_file=f'{path_to_folder}/{folder}.raw.bcf'
            if not os.path.exists(output_file):
                try:
                    command=f'bcftools mpileup -Ou -f {path_to_reference_file}  {path_to_folder}/*.bam | bcftools call --ploidy 1 -mv -Ob -o {path_to_folder}/{folder}.raw.bcf'
                    print('Executing:',command)
                    subprocess.call([command],shell=True)
                except Exception as e:
                    print(e)


    # for files inside individual mouse folder
    # have list of folders containing prefixed files
    # Logic to do it for all folders
    # folders=['mouse_1']
    # for folder in folders:
    # parent_folder = snp.get_file_basename(path_to_files)
    # print(parent_folder)
    # print(path_to_files)
    # command = f'samtools mpileup -uf {path_to_reference_file} {path_to_files}/*{file_extension} > {parent_folder}.raw.bcf' 
    # print('Executing:',command)
    # subprocess.call([command],shell=True)


def main():
    samtools_mpileup('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','sorted.bam','/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/Combined.fasta')

if __name__ == '__main__':
    main()

