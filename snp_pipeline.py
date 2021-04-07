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
                print(f"Creating folder {parent_folder}")
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


def add_file_prefix_to_chrom(path_to_parent_folder, file_extension, path_to_output_files):
    '''Adds prefix of file name to chromosome column in bam file using awk and outputs to sspecified directory. 
    Future work should take sample csv file as input and use samples as prefixes instead of splitting on file name.'''
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


def get_file_basename(path_to_file):
    '''Returns filename from a full path e.g. /external_HDD4/linda/unc_mouse_trial/genomes/mouse_1/mouse_1_merged.bam returns mouse_1_merged.bam'''
    pattern=r'.*\/'
    f=re.sub(pattern,'',path_to_file)
    return f

def get_file_dir(path_to_file):
    '''Returns directory from a full path e.g. /external_HDD4/linda/unc_mouse_trial/genomes/mouse_1/mouse_1_merged.bam returns /external_HDD4/linda/unc_mouse_trial/genomes/mouse_1/'''
    f=os.path.dirname(os.path.realpath(path_to_file))
    return f

def gather_files_by_name(path_to_bam_files, file_extension, path_to_csv, path_to_output_dir):
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
            full_path_input_file=os.path.join(path_to_bam_files,fle)
            full_path_output_file=os.path.join(path_to_output_dir,key,fle)
            # If the file is present in the directory
            if os.path.exists(full_path_input_file):
                try:
                    for sample in samples:
                        if sample in fle:
                            print(f'Sample name: {sample} found in {fle}')
                            output_folder=os.path.join(path_to_output_dir,key)
                            if not os.path.exists(output_folder) and not os.path.exists(full_path_output_file):
                                os.mkdir(output_folder)
                                print(f'Moving file {f} to {full_path_output_file}\n')
                                shutil.move(full_path_input_file, output_folder)
                            elif not os.path.exists(full_path_output_file):
                                print(f'Moving file {f} to {full_path_output_file}\n')
                                shutil.move(full_path_input_file, output_folder)
                            else:
                                print(f'File {f} already exists in {full_path_output_file}\n')
                                os.remove(full_path_input_file)
                except Exception as e:
                    print(e)
                                                                   

def samtools_merge(path_to_bam_files, file_extension, output_dir, template_header_file):
    '''Looks for sam/bam files with specified file_extention and merges them into one large bam file.
    Need to specify header file or else headers from all files get merged into composite header'''
    # Grab header of template bam file, output to temp file to pass to samtools merge
    # Header file needs to be in SAM format?
    #command = f'samtools view -H {template_header_file} > header.sam' 
    template_bam_header_path=template_header_file
    #subprocess.call([command],shell=True)
    files = os.listdir(path_to_bam_files)
    files = sorted(files)
    for folder in os.listdir(path_to_bam_files):
        if os.path.isdir(os.path.join(path_to_bam_files,folder)):
            full_folder_path=os.path.join(path_to_bam_files,folder)
            print(f'Accessing files from {full_folder_path}\n')
            # Merge all files inside folder with file_extension
            command = f'samtools merge -h {template_bam_header_path} {full_folder_path}/{folder}_merged.bam {full_folder_path}/*{file_extension}'
            try:
                print('Executing',command+'\n')
                subprocess.call([command],shell=True)
            except Exception as e:
                print(f'{e}\n')
       

def samtools_idx(path_to_parent_folder,file_extension):
    '''Indexes all files with specified extension within subfolders of parent folder. 
    To open a file using pysam the file must first be indexed.'''
    files=access_subfolder_contents(path_to_parent_folder,file_extension)
    files=sorted(files)
    for f in files:
        command = f'samtools index {f}' # no need to specify output directory, samtools creates index file in same directory as input file
        print(f'Executing: {command}\n')
        try:
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)


def reheader_file(path_to_file):
    '''Iterates through bam input file and finds unique sequence names to place in new header.'''
    input_bam=pysam.AlignmentFile(path_to_file)
    header = input_bam.header.to_dict()

    command= f'samtools view {path_to_file} | awk \'{{print $3}}\' | sort | uniq'
    print(f'Executing {command}')
    output= subprocess.check_output([command],
            stderr=subprocess.STDOUT,
            shell=True)
    print(output.split('\n'))

    #command = samtools view mouse_2_merged.bam | awk '{print $3}' | sort | uniq
    # save output of subprocess call to list

    # for seq in new_head['SQ']:
    #     seq['SN'] = prefix + '_' + seq['SN']
    # #create output BAM with newly defined header
    # with pysam.AlignmentFile(full_output_path, "w", header=new_head) as outf:
   

def filter_bam_by_species(path_to_parent_folder, file_extension, file_with_species):
    '''Filters or extracts lines in a bam file pertaining to a given species from an input txt file and outputs them to a new bam file. '''
    files=access_subfolder_contents(path_to_parent_folder,file_extension)
    files=sorted(files)
    handle=open(file_with_species,'r')
    species=handle.readlines()
    for f in files:
        try:
            input_bam=pysam.AlignmentFile(f)
            path_to_species_folder=get_file_dir(f)
            filename=get_file_basename(f)
            if not os.path.exists(path_to_species_folder):
                os.mkdir(path_to_species_folder)
            elif os.path.exists(path_to_species_folder):
                for s in species:
                    s=s.strip()
                    species_list=s.split('_')
                    if species_list:
                        # Shortening folder name by taking first 3 field identifiers e.g. Enterococcus_phage_A2_length_149431 -> Enterococcus_phage_A2
                        species_folder = '_'.join([species_list[0],species_list[1],species_list[2]])
                        # Naming convention of output file will have parent folder prefix e.g. mouse_1_Enterococcus_phage_A2_merged.bam
                        parent_folder=get_file_basename(path_to_species_folder)
                        output_file=f'{parent_folder}_{species_folder}_merged.bam'
                        output_folder_path=os.path.join(path_to_species_folder,species_folder)
                        output_file_path=os.path.join(output_folder_path,output_file)
                        if os.path.exists(output_file_path):
                            pass
                        else: 
                            # Create new header for file
                            output_bam = pysam.AlignmentFile(output_file_path, "wb", template=input_bam)
                            print(f'Filtering reads containing {s} from {filename} and outputting to {output_file_path}\n')
                            for read in input_bam.fetch():
                                if s in read.reference_name: 
                                    output_bam.write(read)
        except Exception as e: 
            print(e)

def main():
    #add_file_prefix_to_chrom('/external_HDD4/Tom/S.A.3_MouseTrial/Genomes/Round_2','.sorted.bam','/external_HDD4/linda/unc_mouse_trial/genomes/prefixed_bam')
    #gather_files_by_name('/external_HDD4/linda/unc_mouse_trial/genomes/prefixed_bam','pfx.sorted.bam','/external_HDD4/linda/unc_mouse_trial/genomes/mouse_samples.csv','/external_HDD4/linda/unc_mouse_trial/genomes')
    #samtools_merge('/external_HDD4/linda/unc_mouse_trial/genomes/','_pfx.sorted.bam','/external_HDD4/linda/unc_mouse_trial/genomes','/external_HDD4/linda/unc_mouse_trial/genomes/dev/header.sam')
    #samtools_idx('/external_HDD4/linda/unc_mouse_trial/genomes/','merged.bam')

    reheader_file('/external_HDD4/linda/unc_mouse_trial/genomes/mouse_1/mouse_1_merged.bam')

    # Next filter species from merged file and place into corresponding species folder
    #filter_bam_by_species('/external_HDD4/linda/unc_mouse_trial/genomes/','merged.bam','/external_HDD4/linda/unc_mouse_trial/genomes/species_sequences.txt')

    # SNP pipeline
    #reheader_file('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/UNC2FT4158_vs_combined.sorted.bam.bai')


if __name__ == '__main__':
    main()
    

