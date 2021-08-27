#!/usr/bin/env python
__author__ = ('Linda Smith (linda.smith@ucc.ie)')
__version__ = '1.0.0'
__date__ = '10 June 2021'


import microsnp.snptools as snp
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
from collections import defaultdict


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

    2. Filter quality of sequences
    ./mpileup.py filter_vcf_qual /external_HDD4/linda/unc_mouse_trial/snp_pipeline .flt.vcf 20

    3. Compress files with bcftools bgzip and then index 
    ./mpileup.py bcftools_compress_vcf /external_HDD4/linda/unc_mouse_trial/snp_pipeline/ .fltq.vcf
    
    ./mpileup.py bcftools_index_vcf /external_HDD4/linda/unc_mouse_trial/snp_pipeline/ .fltq.vcf.gz

    4. Test if files have been correctly bgzipped with htsfile
    ./mpileup.py htsfile /external_HDD4/linda/unc_mouse_trial/snp_pipeline/ .fltq.vcf.gz

    5. Filter individual species from quality checked bgzipped vcf files
    ./mpileup.py filter_bcf_by_species /external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1 .fltq.vcf.gz

    6 Create species folders to house extracted sequences for each species from vcf files
    ./mpileup.py create_species_folders /external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1 /external_HDD4/linda/unc_mouse_trial/genomes/species_sequences.txt
    #snp.create_species_folders('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1', '/external_HDD4/linda/unc_mouse_trial/genomes/species_sequences.txt')

    7. Move species files to corresponding parent folder
    #snp.move_files_to_species_folder('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1','.fltqs.vcf','/external_HDD4/linda/unc_mouse_trial/genomes/species_sequences.txt')

    # -------------------- 03_06_21 Find common variants between multiple files -----------------------------------------------------------------------------
    
    # Index filtered species files
    #bcftools_compress_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltqs.vcf')
    #bcftools_index_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltq.vcf.gz')

    #bcftools_isec('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltq.vcf.gz','isec_out','+2')

    #bcftools_isec('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1', file_extension, isec_outdir, nfiles_output_pos)

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
    #return vars(p.parse_args())


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

def filter_vcf_qual(path_to_files, file_extension, qual=20):
    '''Filters out snp calls from vcf file with quality below specified qual. Default threshold is 20. Outputs file with fltq.vcf extension.'''
    files = snp.os_walk(path_to_files, file_extension)
    for f in files:
        # get subject_ID
        subject_path=snp.get_file_dir(f)
        subject_clean=snp.get_file_basename(subject_path)
        # get sample_ID
        sample_id=snp.get_output_name(f)
        sample_id_clean=sample_id.split('_')[2]
        command=f'bcftools view -e \'QUAL<{qual}\' {f} > {subject_path}/{subject_clean}_{sample_id_clean}.fltq.vcf'
        try:
            print(command)
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def bcftools_mpileup_multi(path_to_files, file_extension, path_to_reference_file, depth):
    '''Runs bcftools mpileup on multiple bam files (all bam files present in a folder). Path to faidx indexed reference sequence file path_to_reference_file 
    needs to be included.'''
    for folder in os.listdir(path_to_files):
        path_to_folder=os.path.join(path_to_files,folder)
        if os.path.isdir(path_to_folder):
            output_file=f'{path_to_folder}/{folder}.raw.bcf'
            try:
                command=f'bcftools mpileup -Ou -f {path_to_reference_file}  {path_to_folder}/*.bam | bcftools call --ploidy 1 -mv -Ob -o {path_to_folder}/{folder}.raw.bcf'
                print('Executing:',command)
                subprocess.call([command],shell=True)
            except Exception as e:
                print(e)

def index_bcf(path_to_files, file_extension):
    '''Indexes bcf files'''
    files=snp.os_walk(path_to_files, file_extension)
    for f in files:
        command=f'bcftools index {f}'
        try:
            print(command)
            subprocess.call([command],shell=True)
            print('\n')
        except Exception as e:
            print(e)

def filter_bcf_by_species(path_to_files, file_extension):
    '''Finds all unique species/chromosomes from #CHROM column in a bcf/vcf file and outputs the filtered variants to a new bcf file
    with the naming convention <species>.fltq.bcf'''
    files = snp.os_walk(path_to_files, file_extension)
    for f in files:
    # Extract all unique chromosomes/species from bcf file
        command=f'bcftools view {f} |  grep -v "#" | cut -f 1 | uniq | sort'
        unique_chrom=snp.save_process_output(command)
        for chrom in unique_chrom:
            species=chrom.split('_')
            species_prefix = '_'.join([species[0],species[1],species[2]]) 
            outfile = snp.get_output_name(f)
            out_file_name = f'{outfile}_{species_prefix}.fltqs.vcf'
            out_dir = snp.get_file_dir(f)
            full_outfile_path=os.path.join(out_dir,out_file_name)
            filter_command=f'bcftools view {f} --regions {chrom} > {full_outfile_path}'
            print(filter_command)
            print('\n')
            #subprocess.call([filter_command],shell=True)

def get_number_of_variants(path_to_files, file_extension, path_to_output_file=os.getcwd()):
    '''Counts the number of variants for each chromosome in all bcf files or files with specified file_extension
     in path_to_files and outputs result to csv file.'''
    d = {}
    files=snp.os_walk(path_to_files, file_extension)
    for f in files:
        sample_species=snp.get_file_basename(full_path_to_file).split('.')[0]
        command=f'bcftools view {full_path_to_file} | wc -l'
        num_variants=snp.save_process_output(command)
        d[sample_species]=num_variants[0]
    with open(f'{path_to_output_file}/num_variants.csv', 'w') as csv_file:  
        writer = csv.writer(csv_file)
        for key, value in d.items():
            writer.writerow([key, value])
    df=pd.read_csv(f'{path_to_output_file}/num_variants.csv')
    df.columns=['sample_species','variants']
    df=df.sort_values(by='sample_species')
    print(df)
    return df

def bcf_to_vcf(path_to_files, file_extension):
    '''Searches for all bcf files in path with given file_extension and converts them to vcf.'''
    files=snp.os_walk(path_to_files, file_extension)
    for f in files:
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command =f'bcftools view -Oz -o {path_to_output_file}/{output_file_name}.vcf.gz {f}'
        try:
            print('Executing:',command)
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def get_allele_freq(path_to_files, file_extension):
    '''Uses vcftools to output the frequency of alleles for a chromosome in a vcf file.'''
    #vcftools --gzvcf $SUBSET_VCF --freq2 --out $OUT --max-alleles 2
    files=snp.os_walk(path_to_files, file_extension)
    for f in files:
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command=f'vcftools --gzvcf {f} --freq2 --out {path_to_output_file}/{output_file_name} --max-alleles 2'
        try:
            print('Executing:',command,'\n')
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def get_depth(args):
    ''''Calculates the mean depth of coverage per individual.'''
    #vcftools --gzvcf $SUBSET_VCF --depth --out $OUT
    files=snp.os_walk(args.path_to_files, args.file_extension)
    for f in files:
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command=f'vcftools --gzvcf {f} --depth --out {path_to_output_file}/{output_file_name}_depth_profile'
        try:
            print('Executing:',command,'\n')
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def get_avg_depth(args):
    ''''Calculates the average depth of coverage for all samples in vcftools idepth summary.'''
    files=snp.os_walk(args.path_to_files, args.file_extension)
    d = defaultdict(list)
    for f in files:
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command = 'awk \'{ total += $3 } END { print total/NR }\'' + ' ' +f
        avg_depth=snp.save_process_output(command)
        d[output_file_name].append(avg_depth)
    print(d)
    df = pd.DataFrame.from_dict(d)
    print(df.T)

def get_site_mean_depth(path_to_files, file_extension):
    '''Estimates the mean depth of coverage for each site.'''
    #vcftools --gzvcf $SUBSET_VCF --site-mean-depth --out $OUT
    files=snp.os_walk(path_to_files, file_extension)
    for f in files:
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command=f'vcftools --gzvcf {f} --site-mean-depth --out {path_to_output_file}/{output_file_name}_site_mean_depth'
        try:
            print('Executing:',command,'\n')
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def get_site_quality(path_to_files, file_extension):
    '''Extracts the site quality score for each site.'''
    #vcftools --gzvcf $SUBSET_VCF --site-quality --out $OUT
    files=snp.os_walk(path_to_files, file_extension)
    for f in files:
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command=f'vcftools --gzvcf {f} --site-quality --out {path_to_output_file}/{output_file_name}_site_quality'
        try:
            print('Executing:',command,'\n')
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def get_missing_prop_per_site(path_to_files, file_extension):
    ''''Another individual level statistic - calculates the proportion of missing data per sample.'''
    #vcftools --gzvcf $SUBSET_VCF --missing-indv --out $OUT
    files=snp.os_walk(path_to_files, file_extension)
    for f in files:
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command=f'vcftools --gzvcf {f} --missing-indv --out {path_to_output_file}/{output_file_name}_missing_indv'
        try:
            print('Executing:',command,'\n')
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def get_missing_prop_per_indiv(path_to_files, file_extension):
    ''''Gets proportion of missing data per site rather than per individual'''
    #vcftools --gzvcf $SUBSET_VCF --missing-site --out $OUT
    files=snp.os_walk(path_to_files, file_extension)
    for f in files:
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command=f'vcftools --gzvcf {f} --missing-site --out {path_to_output_file}/{output_file_name}_missing_site'
        try:
            print('Executing:',command,'\n')
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def get_het(path_to_files, file_extension):
    '''Calculate heterozygosity and inbreeding coefficient per individual'''
    #vcftools --gzvcf $SUBSET_VCF --het --out $OUT
    files=snp.os_walk(path_to_files, file_extension)
    for f in files:
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command=f'vcftools --gzvcf {f} --het --out {path_to_output_file}/{output_file_name}_het'
        try:
            print('Executing:',command,'\n')
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def prokka_annotate(path_to_files, file_extension):
    '''Annotates files in path file specified file_extension. Note: Sequence ID in FASTA file must be less than 24 characters.'''
    # conda activate prokka
    files=snp.os_walk(path_to_files,file_extension)
    for f in files:
        prefix=snp.get_output_name(f)
        outdir=os.path.join(snp.get_file_dir(f),prefix)
        command = f'prokka {f} --outdir {outdir}/ --prefix {prefix}'
        try:
            print(command)
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def clean_prokka(path_to_references, file_extension):
    '''Deletes lines 2-10 in prokka generated gbk file as these cause formatting issues.'''
    files=snp.os_walk(path_to_references, file_extension)
    for f in files:
        print('Cleaning:',f)
        with open(f,'r+') as handle:
            gbk=handle.readlines()
            handle.seek(0)
            # Delete lines 2:10
            for i in gbk:
                if i not in list(range(2,11)):
                    handle.write(i)
            handle.truncate()
        
def build_snpeff_db_gbk(path_to_snpeff_installation, file_extension):
    '''Note 0_0 SnpEff path to data folder is hardcoded so you have to build reference database inside snpEff program folder.
    Builds a SnpEff reference database for each gbk file with specified file_extension in path_to_references. 
    Note: must first add genomes to snpEff.config and gbk files to /data/programs/snpEff/data.
    '''
    files = snp.os_walk(path_to_snpeff_installation,file_extension)
    for f in files:
        dirname = snp.get_file_dir(f)
        path_to_data_dir = os.path.join(path_to_snpeff_installation,'data',dirname)
        command = f'java -jar /data/programs/snpEff/snpEff.jar build -genbank -v {path_to_data_dir}'
        try:
            print('Executing:',command)
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)
    
def rename_file_extension(path_to_files, file_extension, new_file_extension):
    '''Renames all files in path with specified file_extension to new_file_extension.'''
    files = snp.os_walk(path_to_files, file_extension)
    for f in files:
        basename=snp.get_output_name(f)
        outdir=snp.get_file_dir(f)
        vcf_out_file=basename+new_file_extension
        out_file_path=os.path.join(outdir,vcf_out_file)
        command=f'mv {f} {out_file_path}'
        try:
            print('Executing:',command)
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def bcftools_compress_vcf(path_to_files, file_extension):
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

def bcftools_index_vcf(path_to_files, file_extension):
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
            #print('Executing:',command)
            subprocess.call([command],shell=True)
            print('\n')
        except Exception as e:
            print(e)

def filter_vcf_by_col(path_to_files, file_extension):
    '''DOESNT WORK -- EMPTY FILES?? Uses bcftools on gzipped vcf files e.g. file_extension is vcf.gz '''
    files = snp.os_walk(path_to_files, file_extension)
    for f in files:
        # extract all column names from each vcf file 
        command=f'bcftools query -l {f}' 
        vcf_columns=snp.save_process_output(command)
        print('\n')
        path_to_output_file=snp.get_file_dir(f)
        mouse_path=snp.get_file_dir(path_to_output_file)
        mouse=snp.get_file_basename(mouse_path)
        species = snp.get_file_basename(path_to_output_file)
        species_short = species.split('_')
        spec = '_'.join([species_short[0],species_short[1],species_short[2]])
        for col in vcf_columns:
            column=snp.get_file_basename(col)
            sample_name=snp.get_output_name(column)
            try:
                sample_name=sample_name.split('_')[0]
            except:
                sample_name=sample_name
            output_filename=mouse+'_'+sample_name + '_' + spec+'.vcf.gz'
            #grab mouse ID to put in file_name
            output_file_path=os.path.join(path_to_output_file,output_filename)
            #create a new vcf file containing only that column
            command2=f'bcftools view {f} --regions {col} -Oz -o {os.path.join(path_to_output_file,output_file_path)}'
            try:
                print('Executing:',command2)
                subprocess.call([command2],shell=True)
                print('\n')
            except Exception as e:
                print(e)

def bcftools_isec(path_to_files, file_extension, isec_outdir, nfiles_output_pos):
    '''Finds all files in path with specified file_extension that are inside a folder and applies a minimal bcftools isec on them
    to create intersection and complements of the sets saving the output to isec_outdir/ The set confguration is specified by nfiles_output_pos
    e.g. +2 means find intersection between 2 files or more. Other options: output positions present in this many (=), this many or more (+), this many or fewer (-), or the exact same (~) files.
    Output file name includes subject and species name.'''
    for folder in os.listdir(path_to_files):
        if 'mouse' in folder:
            path_to_mouse_folder=os.path.join(path_to_files,folder)
            for dir_ in os.listdir(path_to_mouse_folder):
                path_to_species_dir=os.path.join(path_to_mouse_folder,dir_) 
                if os.path.isdir(path_to_species_dir):
                    mouse = snp.get_file_basename(path_to_mouse_folder)
                    species = snp.get_output_name(path_to_species_dir)
                    try:
                        command = f'bcftools isec {path_to_species_dir}/*{file_extension} -p {os.path.join(path_to_species_dir,isec_outdir)} -n {nfiles_output_pos} -o {mouse}_{species}.vcf'
                        print(command, '\n')
                        subprocess.call([command],shell=True)
                    except Exception as e:
                        print(e)

def find_max_isec(path_to_files, file_extension, isec_outdir):
    files = snp.os_walk(path_to_files, file_extension)
    path_to_isec_files = []
    for f in files:
        if isec_outdir in f:
            vcf_path = snp.get_file_dir(f)
            path_to_isec_files.append(vcf_path)
    vcf_isecs = []
    for g in path_to_isec_files:
        vcf_isecs.append(max(glob.glob(f'{g}/????.vcf')))
    uniq_vcf_isec = set(vcf_isecs)
    return uniq_vcf_isec

def find_isec_files(path_to_files, file_extension, isec_outdir):
    '''Returns the last vcf file in isec output e.g. 0004.vcf which would contain the intersection of 0001.vcf, 0002.vcf, 0003.vcf'''
    uniq_vcf_isec = find_max_isec(path_to_files, file_extension, isec_outdir)
    for vcf in uniq_vcf_isec:
        isec_dir=snp.get_file_dir(vcf)
        vcf_name=snp.get_file_basename(vcf)
        species=snp.get_file_dir(isec_dir)
        species_base=snp.get_file_basename(species)
        species_l=species_base.split('_')
        try:
            species_short='_'.join([species_l[0],species_l[1],species_l[2]])
        except:
            species_short='_'.join([species_l[0],species_l[1]])
        species_dir=snp.get_file_dir(isec_dir)
        subject_dir=snp.get_file_dir(species_dir)
        subject=snp.get_file_basename(subject_dir)
        # new filename path
        new_filename=f'{isec_dir}/{subject}_{species_short}_{vcf_name}'
        # rename vcf file to include subject_species name
        os.rename(vcf,new_filename)

def find_renamed_isec(path_to_files, file_extension, isec_outdir, isec_cp_dest):
    renamed_isec = []
    files = snp.os_walk(path_to_files, file_extension)
    for f in files:
        if isec_outdir in f:
            if 'mouse' in snp.get_file_basename(f):
                renamed_isec.append(f)
    # copy all renamed isec files to a single location for downloading
    for fle in renamed_isec:
        shutil.copy(fle, isec_cp_dest)

# Since bcftools isec doesn't tell you from which sample the variants are common for
# So need to do a whole other series of steps 
# https://www.biostars.org/p/298361/#298464
# https://www.biostars.org/p/9477018/
    
def bcftools_norm(path_to_files, file_extension, path_to_reference):
    '''Normalise (split multi-allelic calls + left-align indels) each VCF and set a unique identifier for ID field'''
    files = snp.os_walk(path_to_files,file_extension)
    #print(files)
    for f in files:
        file_base = snp.get_output_name(f)
        out_file = f'{file_base}.norm.fltq.vcf.gz'
        #print(out_file)
        file_dir = snp.get_file_dir(f)
        out_file_path = os.path.join(file_dir,out_file)
        try:
            command = f'bcftools norm -m-any --check-ref w -f {path_to_reference} {f} | bcftools annotate -x ID --set-id +\'%CHROM\_%POS\_%REF\_%FIRST_ALT\' -Oz > {out_file_path}'
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

def bcftools_merge_norms(path_to_files, file_extension, path_to_output_files=''):
    '''Merges all normalised files in path for a particular subject.'''
    files = snp.os_walk(path_to_files, file_extension)
    # for folder in os.listdir(path_to_files):
    #     path_to_folder=os.path.join(path_to_files,folder)
    #     #print('HERE')
    #     if os.path.isdir(path_to_folder):
    #         files = snp.os_walk(path_to_folder,file_extension)
    #print(files)
    f = files[0]
    # prepare file output name and path
    species_path=snp.get_file_dir(f)
    species=snp.get_file_basename(species_path)
    subject_path=snp.get_file_dir(species_path)
    subject=snp.get_file_basename(subject_path)
    species_l = species.split('_')
    try:
        species_short='_'.join([species_l[0],species_l[1],species_l[2]])
    except:
        species_short='_'.join([species_l[0],species_l[1]])
    out_name = f'{subject}_{species_short}_merged.norm.fltq.vcf.gz'
    out_path=(os.path.join(species_path,out_name))
    print(out_path)
    if path_to_output_files != '':
        command = f'bcftools merge {path_to_files}/*{file_extension} -Oz > {os.path.join(out_name,out_path)}'
    else:
        command = f'bcftools merge {path_to_files}/*{file_extension} -Oz > {out_path}'
    try:
        print(f'{command}\n') 
        subprocess.call([command],shell=True) 
    except Exception as e:
        print(e)
    command2 = f'tabix -p vcf {out_path}'
    try:
        print(f'{command2}\n') 
        subprocess.call([command2],shell=True) 
    except Exception as e:
        print(e)

def validate_num_sample_files(path_to_files, file_extension, path_to_txt):
    '''Finds the difference in contents of folder compared to list of samples present in txt file'''
    files = snp.os_walk(path_to_files,file_extension)
    handle = open(path_to_txt,'r')
    samples = [line.split('\n') for line in handle.readlines()]
    print(samples)
    # Read through all files, extract sample name ID in file, then see is it in sample list
    for f in files:
        basename = snp.get_file_basename(f)
        file_split=basename.split('_')
        samplename=file_split[2]
        if samplename not in samples:
            print(f)
            # move/delete file  

    
def main():
    #bcftools_isec('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltq.vcf.gz','isec_out','+2')
    #find_isec_files('/external_HDD4/linda/unc_mouse_trial/snp_pipeline','.vcf','isec_out')
    #find_max_isec('/external_HDD4/linda/unc_mouse_trial/snp_pipeline','.vcf','isec_out')
    #find_renamed_isec('/external_HDD4/linda/unc_mouse_trial/snp_pipeline','.vcf','isec_out','/external_HDD4/linda/unc_mouse_trial/snp_pipeline/isec_intersections')

    # Problem with Imtechella compressed vcf files - run again compression and indexing again

    #bcftools_compress_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_15/Imtechella_halotolerans_length_3113269', '.fltqs.vcf')
    #bcftools_index_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_15/Imtechella_halotolerans_length_3113269', '.fltq.vcf.gz')

    # move all files ending in fltq.vcf.gz to all_samples_merged

    #find . -name "*Imtechella*.fltq.vcf.gz"  -exec mv {} /external_HDD4/linda/unc_mouse_trial/snp_pipeline/all_samples_merged/Imtechella_halotolerans \;
    #find . -name "*Imtechella*.fltq.vcf.gz.csi" -exec mv {} /external_HDD4/linda/unc_mouse_trial/snp_pipeline/all_samples_merged/Imtechella_halotolerans \;

    #bcftools_norm('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/all_samples_merged/Imtechella_halotolerans','.fltq.vcf.gz','/external_HDD4/linda/unc_mouse_trial/snp_pipeline/Combined.fasta')
    bcftools_merge_norms('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/all_samples_merged/Imtechella_halotolerans','.norm.fltq.vcf.gz')

    #validate_num_sample_files('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/all_samples_merged/NC_017316.1_Enterococcus_faecalis_OG1RF','.vcf.gz','/external_HDD4/linda/unc_mouse_trial/genomes/dev/all_samples.txt')





    # pars = read_params(sys.argv)
    # print(pars)

    #print(params)
    # --------------------- Pipeline 2.0 14_05_21 ----------------------------------------------------------------------------------------------------------
    # --------------------- SNP calling --------------------------------------------------------------------------------------------------------------------
    #bcftools_mpileup_single('/external_HDD4/linda/unc_mouse_trial/snp_pipeline', '.sorted.bam', '/external_HDD4/linda/unc_mouse_trial/snp_pipeline/Combined.fasta')
    #filter_vcf_qual('/external_HDD4/linda/unc_mouse_trial/snp_pipeline', '.flt.vcf', 20)

    # Compress files with bgzip and then index before filtering can be applied
    #bcftools_compress_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltq.vcf')
    #bcftools_index_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltq.vcf.gz')

    #htsfile('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/all_samples_merged/Imtechella_halotolerans', 'norm.fltq.vcf.gz')

    # Filter fltq.vcf.gz by species 
    #filter_bcf_by_species('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1', '.fltq.vcf.gz')

    # create species folders 
    #snp.create_species_folders('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1', '/external_HDD4/linda/unc_mouse_trial/genomes/species_sequences.txt')

    # move species files to corresponding folder
    #snp.move_files_to_species_folder('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1','.fltqs.vcf','/external_HDD4/linda/unc_mouse_trial/genomes/species_sequences.txt')

    # -------------------- 03_06_21 Find common variants between multiple files -----------------------------------------------------------------------------
    
    # Index filtered species files
    #bcftools_compress_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltqs.vcf')
    #bcftools_index_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltq.vcf.gz')

    #bcftools_isec('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltq.vcf.gz','isec_out','+2')

    #bcftools_isec('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1', file_extension, isec_outdir, nfiles_output_pos)

    # -------------------------------------------------------------------------------------------------------------------------------------------------------

    #bcftools_mpileup_multi('/external_HDD4/linda/unc_mouse_trial/snp_pipeline','sorted.bam','/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/Combined.fasta')
    #index_bcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline','.bcf')
    # --------------------- Filter VCF quality -------------------------------------------------------------------------------------------------------------
    #filter_qual_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline','.bcf',20)
    #index_bcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline','.vcf.gz')
    #filter_bcf_by_species('/external_HDD4/linda/unc_mouse_trial/snp_pipeline','flt.vcf.gz')
    
    #get_number_of_variants('/external_HDD4/linda/unc_mouse_trial/snp_pipeline','flt.bcf','/external_HDD4/linda/unc_mouse_trial/snp_pipeline')

    # --------------------- Create folders for each species inside each mouse folder -----------------------------------------------------------------------
    #snp.create_species_folders('/external_HDD4/linda/unc_mouse_trial/snp_pipeline','/external_HDD4/linda/unc_mouse_trial/genomes/species_sequences.txt') 
    #snp.move_files_to_species_folder('/external_HDD4/linda/unc_mouse_trial/snp_pipeline','.flt.bcf','/external_HDD4/linda/unc_mouse_trial/genomes/species_sequences.txt')

    # --------------------- Convert bcf to vcf  ------------------------------------------------------------------------------------------------------------
    #bcf_to_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline','.flt.bcf') 

    # --------------------- Vcftools statistics  -----------------------------------------------------------------------------------------------------------
    # get_allele_freq('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')
    # get_depth('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')
    # get_site_mean_depth('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')
    # get_site_quality('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')
    # get_missing_prop_per_site('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')
    # get_missing_prop_per_indiv('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')
    # get_het('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')

    # --------------------- Annotate using Prokka ----------------------------------------------------------------------------------------------------------
    #prokka_annotate('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/reference_genomes_db','.fasta')
    #clean_prokka('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/reference_genomes_db','.gbk')

    # --------------------- SnpEff pipeline ----------------------------------------------------------------------------------------------------------------
    #build_snpeff_db_gbk('/data/programs/snpEff','.gbk')

    #snp.create_sample_folders('/external_HDD4/linda/unc_mouse_trial/genomes/mouse_samples.csv','/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2')

    #rename_file_extension('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/', 'fltq.vcf.gz','.vcf')
    #bcftools_compress_vcf('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/', '_fltq.vcf')
    #bcftools_index_vcf('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/', 'fltq.vcf.gz')
    #htsfile('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/', 'fltq.vcf.gz')

    #filter_vcf_by_col('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/','fltq.vcf.gz')


if __name__ == '__main__':
    main()

