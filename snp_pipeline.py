import pandas as pd
import os
import re
import shutil
import subprocess
import pysam


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
                print('Creating folder: ',parent_folder)
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
            print('Prefixing file: ',f)
            print('Adding prefix ',prefix,' to chromosome: ',s)
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

    #        # for f in bam_files:
    #             # Regex with multiple criteria e.g. name match + extension
    #             if(f.find(str(s)==0)) and re.search('\.sorted.bam$', f):
    #                 print(f)
                # if f contains the sample ID and end in file extension copy file to folder it belongs to
                # find its corresponding folder and move it there

def prefix_chr_pysam(path_to_parent_folder, file_extension, path_to_species_file, path_to_output_files):
    '''Uses pysam module to edit a region of the chromosome column.'''
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
            print("Filtering species lines from:",fle,'located at',f)
            print('Appending prefix: ',prefix,'to chromosome:',species)
            print('Outputting file to:',full_output_path+'\n')
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


def main():
    prefix_chr_pysam('/external_HDD4/Tom/S.A.3_MouseTrial/Genomes/Round_2','.sorted.bam','../species_sequences.txt','/external_HDD4/linda/unc_mouse_trial/genomes/prefixed_bam')
    #create_sample_folders('test.csv','/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline')
    #gather_bam_files('test.csv','/external_HDD4/Tom/S.A.3_MouseTrial/Genomes/Round_2')


if __name__ == '__main__':
    main()
    

# Those bam files which were aligned against the combined sequences file just contain an alignment score and other stuff so if I filter the reads from each 