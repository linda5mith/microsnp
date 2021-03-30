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

def samtools_merge(dict_keys,path_to_bam_files):
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


def create_species_folders(path_to_input_file):
    handle=open(path_to_input_file,"r")
    species_file=handle.readlines()
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

# fetch called on bamfile without index
# rb is for bam
# r is for sam

def prefix_chr_pysam(path_to_parent_folder, file_extension, path_to_species_file):
    '''Uses pysam module to edit a region of the chromosome.'''
    handle=open(path_to_species_file,'r')
    species_chr=handle.readlines()
    files=access_folder_contents(path_to_parent_folder,file_extension)
    # regex to extract file name (matches everything before last occurrence of '/')
    pattern=r'.*\/'
    for s in species_chr:
        species=s.split()
        for f in files:
            fle=re.sub(pattern,'',f)
            prefix=fle.split('_')[0]
            print('\n')
            print('FILEEE::: ',fle)
            print('PREFIX::: ',prefix)
            print('\n')
            output_file=fle.split('.')[0]
            input_bam=pysam.AlignmentFile(f,'rb')
            for read in input_bam.fetch(reference=species):
                print(read.reference_id)
                print(read.query_name)
                prefixed_chrom=prefix + '_' +input_bam.get_reference_name(read.reference_id)
                print(prefixed_chrom)

    
# Need to retain folder name / prefix folder name onto read name in bam file
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

def call_snp(path_to_parent_folder,file_extension):
    '''Finds files with specified name/extension and commences samtools snp calling pipeline on them'''
    files=access_subfolder_contents(path_to_parent_folder,file_extension)
    # Add prefix to files

    # SNP pipeline steps:
    # 1. Samtools sort
    samtools_sort(files)
    # 2. Samtools index
    # 3. Samtools mpileup
    # 4. bcftools call
    # 5. Generating vcf file bcftools view my-var.bcf | vcfutils.pl varFilter - > my.var-final.vcf
    
# Need to make samtools helper functions
                   
def samtools_sort(list_of_bam_files):
    for fle in list_of_bam_files:
        #output_directory=os.path.dirname(os.path.abspath(fle))
        output_file=''.join(fle.split('.')[0]+'_sorted.bam')
        print("Sorting:",fle)
        pysam.sort("-o", output_file, fle)

def samtools_idx(list_of_bam_files):
    for fle in list_of_bam_files:
        #output_directory=os.path.dirname(os.path.abspath(fle))
        # Need to find files ending in 'sorted'
        output_file=''.join(fle.split('.')[0]+'_idx.bam')
        print("Sorting:",fle)  
        pysam.index("-o", output_file, fle)

def access_subfolder_contents(path_to_parent_folder,file_extension):
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

def access_folder_contents(path_to_folder,file_extension):
    '''Returns files with specified extension from inside one folder.'''
    files=[]
    for f in os.listdir(path_to_folder):
        if re.search(r'{}$'.format(file_extension), f):
            path_to_file=os.path.join(path_to_folder,f)
            files.append(path_to_file)
    return files                                                      


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
        print('Adding prefix',prefix,'to chromosome column of',fle)
        output_file=fle.split('.')[0]
        output_file_path=os.path.join(path_to_output_files,output_file+'_pfx.sorted.bam')
        print('Outputting file to:',output_file_path)
        if not os.path.exists(output_file_path):
            command= 'samtools view '+ f+' | awk -F\'\\t\' -vOFS=\'\\t\' \'{ $3 = \"'+prefix+'\" $3 }1\' > '+output_file_path
            print('Executing command ',command+'\n')
            subprocess.call(command,shell=True)
            input_bam=pysam.AlignmentFile(f)
            # Adding new header to awk modified bam file
            new_head = input_bam.header.to_dict()
            for seq in new_head['SQ']:
                seq['SN'] = prefix + '_' + seq['SN']
            with pysam.AlignmentFile(output_file_path, "w", header=new_head) as outf:
                for read in input_bam.fetch():
                    outf.write(read)
            input_bam.close()
            outf.close()
                

def move_files_to_folder(path_to_files, path_to_output_files):
    pass

def main():
    add_file_prefix_to_chrom('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline','sorted.bam','../species_sequences.txt','/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/awk_test_prefix')



    #add_file_prefix_to_chrom('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline','.sorted.bam','../species_sequences.txt','/external_HDD4/linda/unc_mouse_trial/genomes/prefixed_bam')
    

    # #samtools view -h mouse_1_Allobacillus_halotolerans_length.bam | sed -e s/Allobacillus_halotolerans_length_2700297/test_Allobacillus_halotolerans_length_2700297/g | head -10



    #call_snp('/external_HDD4/linda/unc_mouse_trial/genomes/','merged.bam')
    #prefix_bam_reads('/external_HDD4/Tom/S.A.3_MouseTrial/Genomes/Round_2/UNC2FT29_vs_combined.sam.bam.sorted.bam','Allobacillus_halotolerans_length_2700297')


    # samtools view /external_HDD4/Tom/S.A.3_MouseTrial/Genomes/Round_2/UNC2FT29_vs_combined.sam.bam.sorted.bam | sed -e 's/Enterococcus_phage_Sw5_length_143759/Test_Enterococcus_phage_Sw5_length_143759/g'

    #call_snp('/external_HDD4/linda/unc_mouse_trial/genomes/','merged.bam')
 

    #clean_species=read_combined_seqs('species_sequences.txt')
    #print(clean_species)
    #index_bam_files('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/mouse_1')
    #index_bams('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/')
    #extract_species_bam('/external_HDD4/linda/unc_mouse_trial/genomes/','species_sequences.txt')
    #samtools_sort('/external_HDD4/linda/unc_mouse_trial/genomes/','.bam')
    # files=access_subfolder_contents('/external_HDD4/linda/unc_mouse_trial/genomes/','merged.bam')
    # print(files)
    

if __name__ == '__main__':
    main()


# Check size of Imtechella_halotolerans_length_3113269 files as run time was longer to extract sequences


