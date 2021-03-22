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
                








def main():
    create_sample_folders('test.csv','/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline')
    gather_bam_files('test.csv','/external_HDD4/Tom/S.A.3_MouseTrial/Genomes/Round_2')


if __name__ == '__main__':
    main()
    

