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
    bcftools_parser.add_argument("path_to_files", help='Path to files')
    bcftools_parser.add_argument("file_extension", help='File extension e.g. .vcf')
    bcftools_parser.add_argument("path_to_reference_file",help='path_to_reference_file')
    bcftools_parser.required = True
    bcftools_parser.set_defaults(func=action)

    args = parser.parse_args()
    print(args)
     

def bcftools_mpileup_single(args):
    '''Runs bcftools mpileup on a single BAM file against the reference. Path to faidx indexed reference sequence path_to_reference_file 
    must be included.'''
    # args=sys.argv
    files = snp.os_walk(args.path_to_files, args.file_extension)
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
#     # pars = read_params(sys.argv)
#     # print(pars)

#     # # create the top-level parser
#     # parser = argparse.ArgumentParser()
#     # subparsers = parser.add_subparsers()

#     # # create the parser for the "bcftools_mpileup_single" command
#     # parser_mpileup_s = subparsers.add_parser('bcftools_mpileup_single')
#     # parser_mpileup_s.add_argument('--path_to_files', type=str)
#     # parser_mpileup_s.add_argument('--file_extension', type=str)
#     # parser_mpileup_s.add_argument('--path_to_reference_file', type=str)
#     # parser_mpileup_s.set_defaults(func=bcftools_mpileup_single)

#     # args = parser.parse_args('bcftools_mpileup_single --path_to_files 1 --file_extension 2 --path_to_reference_file 3'.split())
#     # args.func(args)


if __name__ == '__main__':
    #t0 = time.time()
    main()
    #sys.stderr.write('Elapsed time to run SNP_tools: {} s\n'.format( (time.time()-t0) ) )

