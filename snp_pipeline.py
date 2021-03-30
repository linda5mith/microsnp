import pandas as pd
import numpy as np
import os
import re
import shutil
import subprocess
import pysam
import csv


def create_sample_folders(path_to_csv, path_to_output_folder=os.getcwd()):
    '''Function taking csv as input and creating a folder for each specimen (column header/name) and corresponding subfolders (column values)
     which represent the samples/sequences/timepoints for that specimen.
    Column names become the parent folder with the column values becoming the subfolders'''
    df=pd.read_csv(path_to_csv)
    sample_dict={}
    for item in df.columns:
        sample_dict[item]=df[item].values.tolist()
    keys=sample_dict.keys()
    try:
        for item in keys:
            parent_folder = path_to_output_folder+'/'+item
            if os.path.exists(parent_folder):
                print(parent_folder, ' directory already exists.')
            else:
                os.mkdir(parent_folder)
                print('\n')
                print(f'Creating folder: {parent_folder}')
                subfolders=sample_dict.get(item)
                for f in subfolders:
                    try:  
                        subfolder_path=os.path.join(parent_folder,f)
                        if os.path.exists(subfolder_path):
                            print('Subfolder: ',subfolder_path, ' already exists.')
                        else:
                            os.mkdir(subfolder_path)
                            print('Creating subfolder: ',subfolder_path)
                    except Exception as e:  
                        print('Error: ',e,' occurred trying to create subfolder: ',f)
    except Exception as e:
        print(e)    
    return sample_dict 

def access_folder_contents(path_to_folder,file_extension):
    '''Returns files with specified extension from inside one folder.'''
    files=[]
    for f in os.listdir(path_to_folder):
        if re.search(r'{}$'.format(file_extension), f):
            path_to_file=os.path.join(path_to_folder,f)
            files.append(path_to_file)
    return files  

def access_subfolder_contents(path_to_parent_folder,file_extension):
    '''Specify path to parent folder containing sub directories and the file extension to extract. Returns path to all files
    with desired extension. '''
    files=[]
    for folder in os.listdir(path_to_parent_folder):
        # Need to specify full path to subfolder!
        path_to_subfolders=os.path.join(path_to_parent_folder,folder)
        if os.path.isdir(path_to_subfolders):
            contents=os.listdir(path_to_subfolders)
            for fle in contents:
                if re.search(r'{}$'.format(file_extension),fle):
                    match=os.path.join(path_to_subfolders,fle)
                    files.append(match)
    return files  
    
def add_file_prefix_to_chrom(path_to_parent_folder, file_extension, path_to_species_file, path_to_output_files):
    handle=open(path_to_species_file,'r')
    species=handle.readlines()
    files=access_folder_contents(path_to_parent_folder,file_extension)
    prefixes=[]
    for f in files:
        fle=f.split('/')[-1]
        prefix=fle.split('_')[0]
        # !!!!!! Note for future production 0_0 add regex to extract text after last occurrence of '/'
        output_f=f.split('.')[2]
        output_file=output_f.split('/')[-1]
        for species_read in species:
            s=species_read.strip()
            print(f'Prefixing file: {f}')
            print(f'Adding prefix {prefix} to chromosome: {s}')
            command='samtools view '+ f+' | sed -e '+'s/'+s+'/{}'.format(prefix)+'_'+s+'/g | samtools reheader - '+f+' > '+ path_to_output_files+'/'+output_file+'_pfx.bam'
            print('\n')


def gather_bam_files(path_to_csv, path_to_bam_files, output_dir=os.getcwd()):
    '''Uses column names/headers from csv file and searches inside specified folder for bam files ending
     with a given extension and moves them to their output directory '''

    df=pd.read_csv(path_to_csv)
    bam_files=os.listdir(path_to_bam_files)
    print(type(bam_files))
    # for specimen in df.columns:
    #     samples=(df[specimen].values)
    #     for s in samples:

def add_file_prefix_to_chrom(path_to_parent_folder, file_extension, path_to_species_file, path_to_output_files):
    '''Adds prefix of file name to chromosome column in bam file using awk and outputs to sspecified directory. '''
    handle=open(path_to_species_file,'r')
    species_chr=handle.readlines()
    files=access_folder_contents(path_to_parent_folder,file_extension)
    # regex to extract file name (matches everything before last occurrence of '/')
    pattern=r'.*\/'
    for f in files:
        fle=re.sub(pattern,'',f)
        prefix=fle.split('_')[0]
        print(f'Adding prefix {prefix} to chromosome column of {fle}')
        output_file=fle.split('.')[0]
        output_file_path=os.path.join(path_to_output_files,output_file+'_pfx.sorted.bam')
        print(f'Outputting file to:{output_file_path}')
        if not os.path.exists(output_file_path):
            command= 'samtools view '+ f+' | awk -F\'\\t\' -vOFS=\'\\t\' \'{ $3 = \"'+prefix+'\" $3 }1\' > '+output_file_path
            print(f'Executing command {command}\n')
            subprocess.call(command,shell=True)
            input_bam=pysam.AlignmentFile(f)
            # Adding new header to awk modified bam file
            new_head = input_bam.header.to_dict()
            for seq in new_head['SQ']:
                seq['SN'] = prefix + '_' + seq['SN']
            with pysam.AlignmentFile(output_file_path, "w", header=new_head) as outf:
                for read in input_bam.fetch():
                    outf.write(read)


def prefix_chr_pysam(path_to_parent_folder, file_extension, path_to_species_file, path_to_output_files):
    '''Deprecated for inefficiency: use add_file_prefix_to_chrom
    Uses pysam module to edit a region of the chromosome column. path_to_parent_folder is the path to folder containing the bam/sam files,
    file_extension is the extension of the files you want to edit i.e. 'sorted.bam', path_to_species_file is the path to the file which contains
    all unique species hits or chromosome reference ids that appear in the bam file.'''
    handle=open(path_to_species_file,'r')
    species_chr=handle.readlines()
    files=access_folder_contents(path_to_parent_folder,file_extension)
    # regex to extract file name (matches everything before last occurrence of '/')
    pattern=r'.*\/'
    for s in species_chr:
        species=s.strip()
        for f in files:
            fle=re.sub(pattern,'',f)
            prefix=fle.split('_')[0]
            output_file=fle.split('.')[0]
            full_output_path=path_to_output_files+'/'+output_file+'_prf.sorted.bam'
            print(f'Filtering species lines from:{fle} located at {f}')
            print(f'Appending prefix: {prefix} to chromosome: {species}')
            print(f'Outputting file to: {full_output_path}\n')
            input_bam=pysam.AlignmentFile(f)
            new_head = input_bam.header.to_dict()
            for seq in new_head['SQ']:
                seq['SN'] = prefix + '_' + seq['SN']
            #create output BAM with newly defined header
            with pysam.AlignmentFile(full_output_path, "w", header=new_head) as outf:
                for read in input_bam.fetch():
                    prefixed_chrom = prefix + '_' + read.reference_name
                    a = pysam.AlignedSegment(outf.header)
                    a.query_name = read.query_name
                    a.query_sequence = read.query_sequence
                    a.reference_name = prefixed_chrom
                    a.flag = read.flag
                    a.reference_start = read.reference_start
                    a.mapping_quality = read.mapping_quality
                    a.cigar = read.cigar
                    a.next_reference_id = read.next_reference_id
                    a.next_reference_start = read.next_reference_start
                    a.template_length = read.template_length
                    a.query_qualities = read.query_qualities
                    a.tags = read.tags
                    outf.write(a)


def gather_files_by_name(path_to_bam_files,file_extension,path_to_csv,path_to_output_dir):
    '''Using a csv file for reference, finds all files corresponding to sample name and moves them to their corresponding parent folder.
    If a folder for the subject does not exist it will be created in specified output directory.'''
    files=access_folder_contents(path_to_bam_files,file_extension)
    pattern=r'.*\/'
    df=pd.read_csv(path_to_csv)
    subjects=df.columns
    data_dict={}
    for subject in subjects:
        key = subject
        data_dict.setdefault(key, [])
        for value in df[subject].values:
            if pd.notna(value): 
                data_dict[key].append(value)
    print(f'CSV sample data converted to: {data_dict}\n')
    keys=data_dict.keys()
    for fle in os.listdir(path_to_bam_files):
        f=re.sub(pattern,'',fle)
        for key in keys:
            samples=data_dict.get(key)
            for sample in samples:
                if sample in fle:
                    print(f'SAMPLE:::{sample} found in {fle}\n')
                    output_folder=os.path.join(path_to_output_dir,key)
                    if not os.path.exists(output_folder):
                        os.mkdir(output_folder)
                    else:
                        full_path_to_file=os.path.join(path_to_bam_files,fle)
                        full_path_to_output_file=os.path.join(path_to_bam_files,key,fle)
                        #print(full_path_to_output_file)
                        if not os.path.exists(full_path_to_output_file):
                            print(f'Moving file {f} to {full_path_to_output_file}\n')
                            shutil.move(full_path_to_file, os.path.join(path_to_output_dir,key))
                        else:
                            print(f'{full_path_to_output_file} already exists in {os.path.join(path_to_output_dir,key)}\n')
                                                                   

def samtools_merge(path_to_bam_files, file_extension, path_to_csv, output_dir):
    '''Reads csv containing sample names (values) by subject (column). Looks for sam/bam files with the sample name in the file name and 
    merges these bam files into one large bam file.'''
    files=access_folder_contents(path_to_bam_files,file_extension)
    # base filename
    pattern=r'.*\/'
    df=pd.read_csv(path_to_csv)
    subjects=df.columns
    d=df.to_dict('series')
    for subject in subjects:
        values=d.get(subject)
        for v in values:
            if pd.notna(v):
                for f in files:
                    fle=re.sub(pattern,'',f)
                    print('PRINGINT:::',fle)
                    # If column value is a substring of the filename
                    if v in fle:
                        print(f'{v} IS IN {fle}')
                        # if subject folder does not exist in output path make one
                        print(f'OUTPUT DIR: {os.path.join(output_dir,subject)}')
                        # else move the files there 
                        #output_file_path=os.path.join(output_dir,subject,fle)
                        #print(output_file_path)
                    
                    # if v in f:
                    #     print('\n')
                        #if not os.path.exists(output_file_path):
                            #print('Moving file '+ fle +' to '+ output_file_path)
                            #shutil.move(f, os.path.join(output_dir,subject,fle))
    # Then merge all files with file extension
    # for folder in output folder


def main():
    gather_files_by_name('/external_HDD4/linda/unc_mouse_trial/genomes/prefixed_bam','sorted.bam','/external_HDD4/linda/unc_mouse_trial/genomes/mouse_samples.csv','/external_HDD4/linda/unc_mouse_trial/genomes')
    #samtools_merge('/external_HDD4/linda/unc_mouse_trial/genomes/prefixed_bam','sorted.bam','/external_HDD4/linda/unc_mouse_trial/genomes/mouse_samples.csv','/external_HDD4/linda/unc_mouse_trial/genomes')
    
    #add_file_prefix_to_chrom('/external_HDD4/Tom/S.A.3_MouseTrial/Genomes/Round_2','.sorted.bam','../species_sequences.txt','/external_HDD4/linda/unc_mouse_trial/genomes/prefixed_bam')
    #create_sample_folders('test.csv','/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline')
    #gather_bam_files('test.csv','/external_HDD4/Tom/S.A.3_MouseTrial/Genomes/Round_2')


if __name__ == '__main__':
    main()
    

