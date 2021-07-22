import snptools as snp
import sys
import os
import re
import csv
import subprocess
import pandas as pd
import argparse
import textwrap
import glob
import shutil

def read_params(args):
    FUNCTION_MAP = {'bcftools_mpileup_single' : bcftools_mpileup_single,
                    'bcftools_mpileup_multi' : bcftools_mpileup_multi,
                    'filter_vcf_qual' : filter_vcf_qual,
                    'bcftools_compress_vcf' : bcftools_compress_vcf,
                    'htsfile': htsfile,
                    'filter_bcf_by_species': filter_bcf_by_species,
                    'get_number_of_variants' : get_number_of_variants,
                    'create_species_folders': snp.create_species_folders,
                    'move_files_to_species_folder' : snp.move_files_to_species_folder,
                    'bcftools_isec' : bcftools_isec,
                    'bcf_to_vcf' : bcf_to_vcf,
                    'get_allele_freq' : get_allele_freq,
                    'get_depth' : get_depth,
                    'get_site_mean_depth' : get_site_mean_depth,
                    'get_site_quality' : get_site_quality,
                    'get_missing_prop_per_site' : get_missing_prop_per_site,
                    'get_het' : get_het,
                    'rename_file_extension' : rename_file_extension,
                    'bcftools_isec' : bcftools_isec
                    }

    p = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''    ---------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    Microbial SNP pipeline using bcftools and SNPEff. 

    --------------------- SNP calling pipeline --------------------------------------------------------------------------------------------------------------------
    1. The first step is to create a mpileup of genomes vs. the reference genome. To make a mpileup consisting of one sample file vs. reference genome use:
    python mpileup.py bcftools_mpileup_single /external_HDD4/linda/unc_mouse_trial/snp_pipeline .sorted.bam /external_HDD4/linda/unc_mouse_trial/snp_pipeline/Combined.fasta

    To create a mpileup of multiple samples genomes (located in a single folder) vs. reference use: 
    ./mpileup.py bcftools_mpileup_multi /external_HDD4/linda/unc_mouse_trial/snp_pipeline .sorted.bam /external_HDD4/linda/unc_mouse_trial/snp_pipeline/Combined.fasta

    2. Compress files with bcftools bgzip and then index 
    ./mpileup.py bcftools_compress_vcf /external_HDD4/linda/unc_mouse_trial/snp_pipeline/ .fltq.vcf
    
    ./mpileup.py bcftools_index_vcf /external_HDD4/linda/unc_mouse_trial/snp_pipeline/ .fltq.vcf.gz

    3. Test if files have been correctly bgzipped with htsfile
    ./mpileup.py htsfile /external_HDD4/linda/unc_mouse_trial/snp_pipeline/ .fltq.vcf.gz

    4. Filter individual species from quality checked bgzipped vcf files
    ./mpileup.py filter_bcf_by_species /external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1 .fltq.vcf.gz

    -------------------------------------------------------------------'''),
    epilog='Linda Smith https://github/flannsmith/snp-calling')

    p.add_argument('command', choices=FUNCTION_MAP.keys())
    # if command is 'bcftools_mpileup_single' or 'bcftools_mpileup_multi'
                    

    p.add_argument('--i', metavar='input_dir', type=dir_path, nargs='+', default=None, help='Path to files to process. If parent directory specified, files with given --file_extension are recursively found in subdirectories.')
    p.add_argument('--file_ext', metavar='file_extension', type=str, nargs='+', default=None, help='File extension of files e.g. \'.vcf\'  \'.vcf.gz\'  \'flt.vcf.gz\'')
    p.add_argument('--o', metavar='output_dir', type=dir_path, nargs='+', default=os.getcwd(), help='Path to output files. Default is current directory.')
    p.add_argument('--qual', metavar='qual', type=int, default=None, nargs='+', help='QUAL score to filter vcf file by.')
    p.add_argument('--path_to_ref', metavar='path_to_reference', type=file_path, nargs='+', default=os.getcwd(), help='Path to reference fasta file.')

    arg = p.parse_args()
    func = FUNCTION_MAP[arg.command]
    func()

def dir_path(path_to_dir):
    if os.path.isdir(path_to_dir):
        return path_to_dir
    else:
        raise NotADirectoryError(path_to_dir)

def file_path(path_to_file):
    if os.path.isfile(path_to_file):
        return path_to_file
    else:
        print("Could not open/read file:", path_to_file)
        sys.exit()

def bcftools_mpileup_single(args):
    '''Runs bcftools mpileup on a single BAM file against the reference. Path to faidx indexed reference sequence path_to_reference_file 
    must be included.'''
    # args=sys.argv
    files = snp.os_walk(args.path_to_files, args.file_extension)
    for f in files:
        # get subject_ID
        subject_path=snp.get_file_dir(f)
        subject_clean=snp.get_file_basename(subject_path)
        # get sample_ID
        sample_id=snp.get_file_basename(f)
        sample_id_clean=sample_id.split('_')[0]
        command=f'bcftools mpileup -Ou -f {args.path_to_reference_file}  {f} | bcftools call --ploidy 1 -Ou -mv |  bcftools filter -s LowQual -e \'%QUAL<20\'  > {subject_path}/{subject_clean}_{sample_id_clean}.flt.vcf'
        print(command)
        subprocess.call([command],shell=True)


def main():
    pars = read_params(sys.argv)
    print(pars)