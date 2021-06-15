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
    parser = argparse.ArgumentParser(description='Program for parsing mpileup arguments')
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    subparser = parser.add_subparsers(dest='subparser_name')

    # Create the parser for the bcftools_mpileup_single command
    bcftools_parser = subparser.add_parser('bcftools_mpileup_single', help='run an mpileup')
    bcftools_parser.add_argument("path_to_files", type=str, help='Path to files')
    bcftools_parser.add_argument("file_extension", type=str, help='File extension e.g. .vcf')
    bcftools_parser.add_argument("path_to_reference_file", type=str, help='path_to_reference_file')
    bcftools_parser.required = True
    bcftools_parser.set_defaults(func=bcftools_mpileup_single)

    args = parser.parse_args()
    #print(args) 
    return args.func(args)

def bcftools_mpileup_single(args):
    '''Runs bcftools mpileup on a single BAM file against the reference. Path to faidx indexed reference sequence path_to_reference_file 
    must be included.'''
    # args=sys.argv
    files = snp.os_walk(args.path_to_files, args.file_extension)
    print(files)
    # print(files)
    # for f in files:
    #     # get subject_ID
    #     subject_path=snp.get_file_dir(f)
    #     subject_clean=snp.get_file_basename(subject_path)
    #     # get sample_ID
    #     sample_id=snp.get_file_basename(f)
    #     sample_id_clean=sample_id.split('_')[0]
    #     command=f'bcftools mpileup -Ou -f {args.path_to_reference_file}  {f} | bcftools call --ploidy 1 -Ou -mv |  bcftools filter -s LowQual -e \'%QUAL<20\'  > {subject_path}/{subject_clean}_{sample_id_clean}.flt.vcf'
    #     print(command)
    #     subprocess.call([command],shell=True)
        
    
def main():
    parse_args()

if __name__ == '__main__':
    #t0 = time.time()
    main()
    #sys.stderr.write('Elapsed time to run SNP_tools: {} s\n'.format( (time.time()-t0) ) )

