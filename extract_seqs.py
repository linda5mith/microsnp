import os
import re
import shutil
import subprocess
import pysam

path_to_bam_files = '/external_HDD4/Tom/S.A.3_MouseTrial/Genomes/Round_2'
entries = os.listdir(path_to_bam_files)

# Read mouse sample ID's from file / make mouse dictionary
mouse_IDS = {
    'mouse_1':'UNC2FT197,UNC2FT2113,UNC2FT4153,UNC2FT55198,UNC2FT7533A,UNC2lu2,UNC2lu3',
    'mouse_2':'UNC2FT198,UNC2FT2114,UNC2FT4154,UNC2FT55199,UNC2FT7534A,UNC2lu13',
    'mouse_3':'UNC2FT199,UNC2FT2115,UNC2FT4155,UNC2FT55200,UNC2FT7535A,UNC2lu17,UNC2lu18',
    'mouse_4':'UNC2FT2116,UNC2FT4156,UNC2FT55201',
    'mouse_5':'UNC2FT1101,UNC2FT2117,UNC2FT4157,UNC2FT55202,UNC2FT7537A,UNC2lu29,UNC2lu30',
    'mouse_6':'UNC2FT1102,UNC2FT2118,UNC2FT4158,UNC2FT55203,UNC2lu40,UNC2u37',
    'mouse_7':'UNC2FT1103,UNC2FT2119,UNC2FT4159,UNC2FT7539A,UNC2lu49',
    'mouse_8':'UNC2FT125,UNC2FT29,UNC2FT4145,UNC2FT55190,UNC2FT7540A,UNC2lu22',
    'mouse_9':'UNC2FT126,UNC2FT210,UNC2FT4146,UNC2FT55191,UNC2FT7541A,UNC2lu34',
    'mouse_10':'UNC2FT127,UNC2FT211,UNC2FT4147,UNC2FT55192,UNC2FT7542A,UNC2lu42,UNC2lu43,UNC2lu44',
    'mouse_11':'UNC2FT128,UNC2FT212,UNC2FT4148,UNC2FT55193,UNC2FT7543A,UNC2lu52,UNC2lu53,UNC2lu54',
    'mouse_12':'UNC2FT129,UNC2FT213,UNC2FT4149,UNC2FT55194,UNC2FT7544A',
    'mouse_13':'UNC2FT130,UNC2FT214,UNC2FT4150,UNC2FT55195,UNC2FT7545A,UNC2lu57',
    'mouse_14':'UNC2FT131,UNC2FT4151,UNC2FT55196,UNC2FT7546A',
    'mouse_15':'UNC2FT132,UNC2FT216,UNC2FT4152,UNC2FT55197,UNC2FT7547A,UNC2lu75'
}

# Create folder for individual mouse bam files
def create_folders(mouse_IDS):
    dict_keys=mouse_IDS.keys()
    path = os.getcwd()
    try:
        for item in dict_keys:
            mouse_path = path+'/'+item
            if os.path.exists(mouse_path):
                print(mouse_path, ' directory already exists')
            else:
                os.mkdir(mouse_path)       
    except OSError:
        print ("Creation of the directory failed")


# For each key in dict extract sample IDs, find files and insert into corresponding directory
def move_bam_files(dict_keys):
    for item in dict_keys:
        mouse_IDS_list=mouse_IDS.get(item).split(',')
        for bam in entries:
            if (bam.split('_')[0]) in mouse_IDS_list:
                # Extract only the files ending in 'sorted.bam' 
                match=(re.search("\.sorted.bam$", bam))
                if match:
                    print(bam)
                    # copy files into item folder
                    shutil.copy(path_to_bam_files+'/'+bam,os.path.join(path, item))
    

# merge bam files inside each folder to create combined bam file

# samtools merge /external_HDD4/linda/unc_mouse_trial/bacterial_genomes/mouse_11/mouse_1/mouse_1_merged.bam input1.bam input2.bam
# samtools merge mouse_1_merged.bam *.bam

def samtools_merge(dict_keys):
    mouse_path = os.getcwd()
    for item in dict_keys:
        path=os.path.join(mouse_path,item)
        command='samtools merge '
        c2=command+path+'_merged.bam *.bam'
        # files=os.listdir(path)
        # for f in files:
        #     f=f+' '
        #     #f=mouse_path+'/'+f+' '
        #     c2+=f
        #     #print(c2)
        print('Executing ',c2)
        subprocess.call([command],shell=True)
        print('Finished ', c2)
    return 'Finished'
    

# Next step extract specific species from merged file to make a new file for it
# Feed into SNP 

def read_combined_seqs(path_to_file):
    combined_seqs = open(path_to_file,'r')
    f=combined_seqs.read()
        
    # Does bam file contain the '>' character - No
    # Strip leading '>' from text 
    clean_seqs=[]
    for line in f:
        line=line.strip('>')
        clean_seqs.append(line)

    with open("species_sequences.txt","w") as f:
        for line in clean_seqs:
            f.writelines(line)
    handle=open("species_sequences.txt","r")
    f=handle.readlines()
    return f


handle=open("species_sequences.txt","r")
species_file=handle.readlines()
#print(species_file)

def create_species_folders():
    mouse_path = os.getcwd()
    for folder in os.listdir(mouse_path):
        if os.path.isdir(folder):
            # new path is mouse_path/folder
            for line in species_file:
                #print(line)
                #Keep only first 3 fields e.g. Enterococcus_phage_A2
                species_list=line.split('_')
                species_str = '_'.join([species_list[0],species_list[1],species_list[2]])
                mouse_folder=os.path.join(mouse_path, folder)
                species_folder=os.path.join(mouse_folder,species_str)
                os.mkdir(species_folder)


# Merged file needs to be indexed before it can be read + filtered with pysam
# test_file = '/external_HDD4/linda/unc_mouse_trial/genomes/mouse_1/mouse_1_merged.bam'
# pysam.index(test_file)


def index_bam_files(path_to_bam_files):
    '''Indexes bam files ending in sorted.bam. Need to specify path to folder containing files.'''
    if os.path.isdir(path_to_bam_files):
        files=os.listdir(path_to_bam_files)
        print(files)
        for f in files:
            match=(re.search("\.sorted.bam$", f))
            if match:
                try:
                    print('Found: ',f)
                    full_file_path=os.path.join(path_to_bam_files,f)
                    # To add: Check if indexed file already exists
                    print('Indexing: ',full_file_path)
                    pysam.index(full_file_path)
                except Exception as e:
                    print(e)


def index_bams(path_to_folders):
    '''Recursively indexes bam files ending in sorted.bam. Need to specify path to parent directory to search within subfolders for files.'''
    for folder in os.listdir(path_to_folders):
        if os.path.isdir(folder):
            print('\n')
            print("Entering folder: ",folder)
            files = os.listdir(folder)
            print(files)
            folder_path=os.path.join(path_to_folders,folder)
            for f in files:
                match=(re.search("\.sorted.bam$", f))
                if match:
                    try:
                        print('Found sorted file: ',f)
                        full_file_path=os.path.join(folder_path,f)
                        print('Indexing file: ',full_file_path)
                        pysam.index(full_file_path)
                    except Exception as e:
                        print(e)


def extract_species_bam(path_to_folders, file_with_species):
    '''Filters or extracts lines pertaining to a given species input txt file and outputs them to a new bam file. 
    Future work is to add parameters to specifiy path_to_output files to'''

    contents=os.listdir(path_to_folders)
    # Input file containing species to extract individual reads for
    handle=open(file_with_species,'r')
    f=handle.readlines()
    for species in f:
        s=species.strip()
        species_list=s.split('_')
        if s:
            species_folder = '_'.join([species_list[0],species_list[1],species_list[2]])
            for folder in contents:
                if os.path.isdir(folder):
                    path=os.path.join(path_to_folders,folder)
                    subfolder_content=os.listdir(folder)
                    species_folder_path=os.path.join(path,species_folder)
                    for fle in subfolder_content:
                        match=(re.search("\_merged.bam$", fle))
                        if match:
                            # Grab full path to bam file
                            path_to_bam=os.path.join(path,fle)
                            # Specify naming convention for filtered output bam file
                            output_name=folder+'_'+species_folder+'.bam'
                            # Full path for outputting file to 
                            output_file_path=os.path.join(species_folder_path,output_name)
                            # Search through bam file and fetch lines containing species s
                            input_bam=pysam.AlignmentFile(path_to_bam,'rb')
                            output_bam = pysam.AlignmentFile(output_file_path, "wb", template=input_bam)
                            print('Filtering species ',s, 'from file: ',fle)
                            print('Outputting reads to ',output_file_path)
                            print('\n')
                            for read in input_bam.fetch(s):
                                output_bam.write(read)
                            input_bam.close()
                            output_bam.close()

def call_snp(path_to_folder):
    '''Finds files with specified name/extension and commences samtools snp calling pipeline on them'''

    # SNP pipeline steps:
    # 1. Samtools sort
    # 2. Samtools index
    # 3. Samtools mpileup
    # 4. bcftools call
    # 5. Generating vcf file bcftools view my-var.bcf | vcfutils.pl varFilter - > my.var-final.vcf
    
# Need to make samtools helper functions

def access_subfolder_contents(path_to_parent_folder,file_extension):
    for folder in os.listdir(path_to_parent_folder):
        if os.path.isdir(folder):
            contents=os.listdir(folder)
            # Need to specify full path to subfolder!
            path=os.path.join(path_to_parent_folder,folder)
            for subfolder in contents:
                full_path=os.path.join(path,subfolder)
                if os.path.isdir(full_path):
                    subcontents=os.listdir(full_path)
                    contents=[]
                    for f in subcontents:
                        if re.search(r'\.{}$'.format(file_extension),f):
                            full_path_to_content=os.path.join(full_path,f)
                            contents.append(full_path_to_content)
    return contents
                                                    
                      
def samtools_sort(path_to_parent_folder,file_extension):
    access_subfolder_contents(path_to_parent_folder,file_extension)
    #print(contents)




def main():
    #clean_species=read_combined_seqs('species_sequences.txt')
    #print(clean_species)
    #index_bam_files('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/mouse_1')
    #index_bams('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/')
    
    #extract_species_bam('/external_HDD4/linda/unc_mouse_trial/genomes/','species_sequences.txt')
    #samtools_sort('/external_HDD4/linda/unc_mouse_trial/genomes/','.bam')
    access_subfolder_contents('/external_HDD4/linda/unc_mouse_trial/genomes/','bam')

if __name__ == '__main__':
    main()


# Check size of Imtechella_halotolerans_length_3113269 files as run time was longer to extract sequences