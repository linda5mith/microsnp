import snptools as snp
import mpileup as mp
import sys
import os
import re
import csv
import subprocess
import pandas as pd
import argparse
import textwrap
import time
from collections import defaultdict


def parse_args():
    # Top-level parser
    parser = argparse.ArgumentParser(description='Variant calling program for microbial metagenomics.\n\
        Python wrapper around bcftools functions for ease of handling multiple samples and species. To understand parameters for individual functions e.g. mpileup_single run\n: python microsnp.py mpileup_single -h')
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    subparser = parser.add_subparsers(dest='action')

    # Parent parser for functions that demand path_to_files, file_extension, path_to_ref, path_to_output be defined
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("path_to_files", type=str, help='Path to files')
    parent_parser.add_argument("file_extension", metavar='file_ext', type=str, help='File extension e.g. \'.bam\' \'.sorted.bam\' \'.vcf\' \'fltq.vcf.gz\'')

    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument("path_to_output", metavar='out_dir', type=str, help='Path to direct output to.')

    mpileup_parser = argparse.ArgumentParser(add_help=False)
    mpileup_parser.add_argument("path_to_ref", metavar='path_to_ref', type=str, help='Path to indexed reference fasta file.')
    #mpileup_parser.add_argument("--path_to_output", metavar='out_dir', type=str, help='Path to output mpileups.')

    dp_parser = argparse.ArgumentParser(add_help=False)
    dp_parser.add_argument("DP", metavar='DP', type=int, help='Minimum depth to filter FORMAT/DP in vcf files.')
    dp_parser.add_argument("expression", metavar='expr', type=str, help='MIN to return only samples that match the minimum DP entered. AVG to return samples that have a depth that equates to the average DP across all samples.')

    gq_parser = argparse.ArgumentParser(add_help=False)
    gq_parser.add_argument("GQ", metavar='GQ', type=int, help='Minimum genotype quality to filter FORMAT/GQ in vcf files.')

    # Create the parser for the mpileup_single command
    bcftools_single = subparser.add_parser('mpileup_single',parents=[parent_parser,mpileup_parser,output_parser], help='Run an mpileup on a single sample vs. a reference genome.',
    description='Example: python argparse.py mpileup_single /external_HDD4/linda/unc_mouse_trial/bam_files .sorted.bam /external_HDD4/linda/unc_mouse_trial/snp_pipeline/reference.fasta /external_HDD4/linda/unc_mouse_trial/mpileups')
    bcftools_single.set_defaults(func=mpileup_single)

    # Create the parser for the mpileup_multi command
    bcftools_multi = subparser.add_parser('mpileup_multi',parents=[parent_parser,mpileup_parser,output_parser], help='run an mpileup on multiple samples vs. a reference genome.',
    description='Example: python argparse.py mpileup_multi /external_HDD4/linda/unc_mouse_trial/bam_files .sorted.bam /external_HDD4/linda/unc_mouse_trial/snp_pipeline/reference.fasta /external_HDD4/linda/unc_mouse_trial/mpileups')
    bcftools_multi.set_defaults(func=mpileup_multi)

    # Need to add description variable with example of usage
    bcftools_compress = subparser.add_parser('compress_vcf',help='Compresses all vcf files in path ending in file_extension to BGZF-compressed variant calling data.')
    bcftools_compress.add_argument("path_to_files", type=str, help='Path to files')
    bcftools_compress.add_argument("file_extension",metavar='file_ext', type=str, help='File extension e.g. .vcf')
    bcftools_compress.required = True

    bcftools_index = subparser.add_parser('index_vcf',parents=[parent_parser],help='Indexes all files with specified file_extension in path.')
    bcftools_index.required = True
    bcftools_index.set_defaults(func=index_vcf)

    htsfile_parser = subparser.add_parser('htsfile',parents=[parent_parser],help='Output compression status for all files in path_to_files with specified file_extension')
    htsfile_parser.required = True
    htsfile_parser.set_defaults(func=htsfile)

    filterpass =  subparser.add_parser('filter_pass',parents=[parent_parser],description='Filters variants from vcf which passed filter to create a new file.')
    filterpass.required = True
    filterpass.set_defaults(func=filter_pass)

    filterdpgq =  subparser.add_parser('filter_dp_gq',parents=[parent_parser,dp_parser,gq_parser],description='Filters variants that have DP > min_dp and GQ > min_gq.')
    filterdpgq.required = True
    filterdpgq.set_defaults(func=filter_dp_gq)

    filterdp =  subparser.add_parser('filter_dp',parents=[parent_parser,dp_parser],description='Filters variants across all values and samples that have DP > DP specified for given expression e.g. EXPR(FORMAT/DP)>DP\
    MIN(FORMAT/DP)>10, AVG(FORMAT/DP)>10')
    filterdp.required = True
    filterdp.set_defaults(func=filter_dp)

    filterspecies = subparser.add_parser('filter_species',parents=[parent_parser],description='Finds all unique species/chromosomes from #CHROM column in a vcf file and outputs the filtered variants to a new vcf file with the naming convention <species>.fltq.vcf')
    filterspecies.required = True
    filterspecies.set_defaults(func=filter_species)

    bcftoolsnorm = subparser.add_parser('bcftools_norm',parents=[parent_parser,mpileup_parser],description='Normalise (split multi-allelic calls + left-align indels) each VCF and set a unique identifier for ID field')
    bcftoolsnorm.required = True
    bcftoolsnorm.set_defaults(func=bcftools_norm)

    bcftoolsmergennorms = subparser.add_parser('bcftools_merge_norms',parents=[parent_parser,output_parser],description='Merges all normalised species files in path_to_files to create normalized, merged vcf for all subjects and samples.')
    bcftoolsmergennorms.required = True
    bcftoolsmergennorms.set_defaults(func=bcftools_merge_norms)

    vcftoolsdepth = subparser.add_parser('get_depth',parents=[parent_parser],help='Output vcftools depth report of file.')
    vcftoolsdepth.required = True
    vcftoolsdepth.set_defaults(func=mp.get_depth)

    avgdepth = subparser.add_parser('get_avg_depth',parents=[parent_parser],help='Get average depth for all samples from vcftools idepth summary file.')
    avgdepth.required = True
    avgdepth.set_defaults(func=mp.get_avg_depth)

    args = parser.parse_args()
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
            -f {args.path_to_ref} {f} | bcftools call --ploidy 1 --multiallelic-caller --format-fields GQ,GP --variants-only -Ou -mv |  bcftools filter -s LowQual -e \'%QUAL<20\'  > {args.path_to_output}/{subject_clean}_{sample_id_clean}.flt.vcf'
        print(command)
        try:
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def mpileup_multi(args): #path_to_files, file_extension, path_to_reference_file
    '''Runs bcftools mpileup on multiple bam files (all bam files present in a folder). Path to faidx indexed reference sequence file path_to_reference_file 
    needs to be included.'''
    files = snp.os_walk(args.path_to_files, args.file_extension)
    folders = []
    for f in files:
        # get subject_ID
        dir_path=snp.get_file_dir(f)
        folders.append(dir_path)
    folder_paths=set(folders)
    #print(len(folder_paths))
    for folder in folder_paths:
        dirname=snp.get_file_basename(folder)
        try:
            command=f'bcftools mpileup -Ou --max-depth 8000 --redo-BAQ --min-MQ 30 --min-BQ 30 --per-sample-mF\
            --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR\
            -f {args.path_to_ref} {folder}/*{args.file_extension} | bcftools call --ploidy 1 --multiallelic-caller --format-fields GQ,GP --variants-only -Ou -mv |  bcftools filter -s LowQual -e \'%QUAL<20\'  > {args.path_to_output}/{dirname}.flt.vcf'
            print(command+'\n')
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def filter_pass(args):
    '''Filters only variants that have passed the applied filter'''
    files = snp.os_walk(args.path_to_files, args.file_extension)
    for f in files:
        basename=snp.get_output_name(f)
        out_dir=snp.get_file_dir(f)
        command = f'bcftools view -Oz -f PASS {f} > {out_dir}/{basename}.fltq.vcf.gz'
        try:
            print(command)
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def filter_dp_gq(args):
    '''Filters variants that have DP > 10 and GQ > 15'''
    files = snp.os_walk(args.path_to_files, args.file_extension)
    for f in files:
        basename=snp.get_output_name(f)
        out_dir=snp.get_file_dir(f)
        command = f'bcftools view -Oz -i  \'MIN(FMT/DP)>{args.DP} & MIN(FMT/GQ)>{args.GQ}\' {f} > {out_dir}/{basename}.merged.norm.fltdpgq.vcf.gz'
        try:
            print('Executing:',command)
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def filter_dp(args):
    '''Filters variants that have DP > min_DP'''
    files = snp.os_walk(args.path_to_files, args.file_extension)
    expr = args.expression
    expr = expr.upper()
    valid_expr = ['MIN','AVG']
    if expr not in valid_expr:
        print(f'Expression: {expr} not valid.')
    else:
        for f in files:
            basename=snp.get_output_name(f)
            out_dir=snp.get_file_dir(f)
            command = f'bcftools view -Oz -i  \'{expr}(FORMAT/DP)>{args.DP}\' {f} > {out_dir}/{basename}.merged.norm.fltdp.vcf.gz'
            try:
                print('Executing:',command)
                subprocess.call([command],shell=True)
            except Exception as e:
                print(e)

def compress_vcf(args):
    '''Compresses all vcf files in path ending in file_extension to BGZF-compressed variant calling data.'''
    files = snp.os_walk(args.path_to_files, args.file_extension)
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

def index_vcf(args):
    '''Indexes all files with specified file_extension in path.'''
    files = snp.os_walk(args.path_to_files, args.file_extension)
    for f in files:
        command=f'bcftools index {f}'
        try:
            print(command)
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def htsfile(args):
    '''Output compression status all files in path_to_files with specified file_extension.'''
    files = snp.os_walk(args.path_to_files, args.file_extension)
    for f in files:
        command=f'htsfile {f}'
        try:
            print(command)
            subprocess.call([command],shell=True)
            print('\n')
        except Exception as e:
            print(e)

def filter_species(args):
    '''Finds all unique species/chromosomes from #CHROM column in a bcf/vcf file and outputs the filtered variants to a new bcf file
    with the naming convention <species>.fltq.bcf'''
    files = snp.os_walk(args.path_to_files, args.file_extension)
    for f in files:
    # Extract all unique chromosomes/species from bcf file
        command=f'bcftools view {f} |  grep -v "#" | cut -f 1 | uniq | sort'
        unique_chrom=snp.save_process_output(command)
        for chrom in unique_chrom:
            #print(chrom)
            species=chrom.split('_')
            species_prefix = '_'.join([species[0],species[1],species[2]]) 
            outfile = snp.get_output_name(f)
            out_file_name = f'{outfile}_{species_prefix}.fltqs.vcf.gz'
            out_dir = snp.get_file_dir(f)
            full_outfile_path=os.path.join(out_dir,out_file_name)
            filter_command=f'bcftools view -Oz {f} --regions {chrom} -o {full_outfile_path}'
            try:    
                print(f'{filter_command}\n')
                subprocess.call([filter_command],shell=True)
            except Exception as e:
                print(e)

def bcftools_norm(args):
    '''Normalise (split multi-allelic calls + left-align indels) each VCF and set a unique identifier for ID field'''
    files = snp.os_walk(args.path_to_files,args.file_extension)
    #print(files)
    for f in files:
        file_base = snp.get_output_name(f)
        out_file = f'{file_base}.norm.fltqs.vcf.gz'
        #print(out_file)
        file_dir = snp.get_file_dir(f)
        out_file_path = os.path.join(file_dir,out_file)
        try:
            command = f'bcftools norm -m-any --check-ref w -f {args.path_to_ref} {f} | bcftools annotate -x ID --set-id +\'%CHROM\_%POS\_%REF\_%FIRST_ALT\' -Oz > {out_file_path}'
            print(f'{command}')
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)
        try:
            command2 = f'tabix -p vcf {out_file_path}'
            print(f'{command2}\n')
            subprocess.call([command2],shell=True)
        except Exception as e:
            print(e)

def bcftools_merge_norms(args):
    '''Merges all normalised files in path for a particular species.'''
    files = snp.os_walk(args.path_to_files, args.file_extension)
    # Find all file names by species and group them
    groups = defaultdict(list)
    species=[]
    for f in files:  
        basename = snp.get_output_name(f)
        out_dir = snp.get_file_dir(f)
        try:
            filename_list = basename.split('_')
            #print(filename_list)
            species_ID = '_'.join([filename_list[2],filename_list[3],filename_list[4]])
            species.append(species_ID)
        except Exception as e:
            pass
        finally:
            species_ID = '_'.join([filename_list[2],filename_list[3]])
            #print(species_ID)
            species.append(species_ID)
    unique_species=set(species)
    print(unique_species)
    if args.path_to_output:
        for s in unique_species:
            try:
                command=f'bcftools merge {args.path_to_output}/*{s}.norm.fltqs.vcf.gz -Oz > {args.path_to_output}/{s}.merged.norm.fltqs.vcf.gz'
                print(command+'\n')
                subprocess.call([command],shell=True)
            except Exception as e:
                print(e)
                command=f'bcftools merge {args.path_to_output}/*{s}*.norm.fltqs.vcf.gz -Oz {args.path_to_output}/{s}.merged.norm.fltqs.vcf.gz'
                print(command+'\n')
                subprocess.call([command],shell=True)
                
def main():
    parse_args()

if __name__ == '__main__':
    main()

