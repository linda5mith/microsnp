import snptools as snp
import sys
import os
import re
import csv
import subprocess
import pandas as pd
import argparse
import textwrap
import time

# python argparse_test.py /external_HDD4/linda/unc_mouse_trial/argparse_test .sorted.bam /external_HDD4/linda/unc_mouse_trial/snp_pipeline/Combined.fasta

def action(args):
    print(args.path_to_files)
    print(args.file_extension)
    print(args.path_to_reference_file)

def parse_args():
    # Top-level parser
    parser = argparse.ArgumentParser(description='Variant calling program for microbial metagenomics.\n\
        To understand parameters for individual functions e.g. mpileup_single run\n: python argparse_test.py mpileup_single -h')
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    subparser = parser.add_subparsers(dest='action')

    # Parent parser for functions that demand path_to_files, file_extension, path_to_ref, path_to_output be defined
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("path_to_files", type=str, help='Path to files')
    parent_parser.add_argument("file_extension", metavar='file_ext', type=str, help='File extension e.g. \'.bam\' \'.sorted.bam\'')

    mpileup_parser = argparse.ArgumentParser(add_help=False)
    mpileup_parser.add_argument("path_to_ref", metavar='path_to_ref', type=str, help='Path to indexed reference fasta file.')
    mpileup_parser.add_argument("--path_to_output", metavar='out_dir', type=str, help='Path to output mpileups.')

    # Create the parser for the mpileup_single command
    bcftools_single = subparser.add_parser('mpileup_single',parents=[parent_parser,mpileup_parser], help='Run an mpileup on a single sample vs. a reference genome.',
    description='Example: python argparse.py mpileup_single /external_HDD4/linda/unc_mouse_trial/bam_files .sorted.bam /external_HDD4/linda/unc_mouse_trial/snp_pipeline/reference.fasta /external_HDD4/linda/unc_mouse_trial/mpileups')
    bcftools_single.set_defaults(func=mpileup_single)

    # Create the parser for the mpileup_multi command
    bcftools_multi = subparser.add_parser('mpileup_multi',parents=[parent_parser,mpileup_parser], help='run an mpileup on multiple samples vs. a reference genome.',
    description='Example: python argparse.py mpileup_multi /external_HDD4/linda/unc_mouse_trial/bam_files .sorted.bam /external_HDD4/linda/unc_mouse_trial/snp_pipeline/reference.fasta /external_HDD4/linda/unc_mouse_trial/mpileups')
    bcftools_multi.set_defaults(func=mpileup_multi)

    bcftools_compress = subparser.add_parser('compress_vcf',help='Compresses all vcf files in path ending in file_extension to BGZF-compressed variant calling data.')
    bcftools_compress.add_argument("path_to_files", type=str, help='Path to files')
    bcftools_compress.add_argument("file_extension",metavar='file_ext', type=str, help='File extension e.g. .vcf')
    bcftools_compress.required = True

    bcftools_index = subparser.add_parser('index_vcf',parents=[parent_parser],help='Indexes all files with specified file_extension in path.')
    bcftools_index.required = True
    bcftools_index.set_defaults(func=index_vcf)

    htsfile_parser = subparser.add_parser('htsfile',parents=[parent_parser],help='Output compression status all files in path_to_files with specified file_extension')
    htsfile_parser.required = True
    htsfile_parser.set_defaults(func=htsfile)

    args = parser.parse_args()
    #print(args) 
    return args.func(args)

def mpileup_single(args):
    '''Runs bcftools mpileup on a single BAM file against the reference. Path to faidx indexed reference sequence path_to_reference_file 
    must be included.'''
    files = snp.os_walk(args.path_to_files, args.file_extension)
    #print(files)
    for f in files:
        # get subject_ID
        subject_path=snp.get_file_dir(f)
        subject_clean=snp.get_file_basename(subject_path)
        # get sample_ID
        sample_id=snp.get_file_basename(f)
        sample_id_clean=sample_id.split('_')[0]
        command=f'bcftools mpileup -Ou --max-depth 8000 --redo-BAQ --min-MQ 30 --min-BQ 30 --per-sample-mF \
            --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
            -f {args.path_to_ref} {f} | bcftools call --ploidy 1 --multiallelic-caller --variants-only -Ou -mv |  bcftools filter -s LowQual -e \'%QUAL<20\'  > {args.path_to_output}/{subject_clean}_{sample_id_clean}.flt.vcf'
        print(command)
        try:
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def mpileup_multi(args): #ath_to_files, file_extension, path_to_reference_file, depth
    '''Runs bcftools mpileup on multiple bam files (all bam files present in a folder). Path to faidx indexed reference sequence file path_to_reference_file 
    needs to be included.'''
    files = snp.os_walk(args.path_to_files, args.file_extension)
    folders = []
    for f in files:
        # get subject_ID
        dir_path=snp.get_file_dir(f)
        folders.append(dir_path)
    folder_paths=set(folders)
    for folder in folder_paths:
        dirname=snp.get_file_basename(folder)
        try:
            #command=f'bcftools mpileup -Ou --max-depth 8000 --min-MQ 30 --min-BQ 30 -f {args.path_to_ref}  {f} | bcftools call --ploidy 1 -Ou -mv |  bcftools filter -s LowQual -e \'%QUAL<20\'  > {args.path_to_output}/{subject_clean}_{sample_id_clean}.flt.vcf'
            command=f'bcftools mpileup -Ou --max-depth 8000 --min-MQ 30 --min-BQ 30 -f {args.path_to_ref} {folder}/*{args.file_extension} | bcftools call --ploidy 1 -mv -Ou |  bcftools filter -s LowQual -e \'%QUAL<20\'  > {args.path_to_output}/{dirname}.flt.vcf'
            print('Executing:',command+'\n')
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def compress_vcf(path_to_files, file_extension):
    '''Compresses all vcf files in path ending in file_extension to BGZF-compressed variant calling data.'''
    files = snp.os_walk(path_to_files, file_extension)
    for f in files:
        basename=snp.get_output_name(f)
        outdir=snp.get_file_dir(f)
        vcf_out_file=basename+'.fltq.vcf.gz'
        out_file_path=os.path.join(outdir,vcf_out_file)
        command=f'bcftools view -Oz -o {out_file_path} {f}'
        try:
            print('Executing:',command)
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def index_vcf(path_to_files, file_extension):
    '''Indexes all files with specified file_extension in path.'''
    files = snp.os_walk(path_to_files, file_extension)
    for f in files:
        command=f'bcftools index {f}'
        try:
            print('Executing:',command)
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def htsfile(path_to_files, file_extension):
    '''Output compression status all files in path_to_files with specified file_extension.'''
    files = snp.os_walk(path_to_files, file_extension)
    for f in files:
        command=f'htsfile {f}'
        try:
            print(command)
            subprocess.call([command],shell=True)
            print('\n')
        except Exception as e:
            print(e)



def main():
    parse_args()

if __name__ == '__main__':
    #t0 = time.time()
    main()
    #sys.stderr.write('Elapsed time to run SNP_tools: {} s\n'.format( (time.time()-t0) ) )

