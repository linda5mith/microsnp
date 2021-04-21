import snptools as snp
import os
import subprocess
import re

#  bcftools view var.raw.bcf.gz | bcftools filter -e'"Enterococcus"' 

# Count the number of variants by species and output to file 
# bcftools view var.raw.bcf.gz | bcftools filter -e'"Enterococcus"' | wc -l 

# MOUSE 10 too many files ---> have a look

def bcftools_mpileup(path_to_files, file_extension, path_to_reference_file):
    '''Run bcftools mpileup on the bam contents of multiple folders. Path to faidx indexed reference sequence file path_to_reference_file 
    needs to be included.'''
    for folder in os.listdir(path_to_files):
        path_to_folder=os.path.join(path_to_files,folder)
        if os.path.isdir(path_to_folder):
            output_file=f'{path_to_folder}/{folder}.raw.bcf'
            if not os.path.exists(output_file):
                try:
                    command=f'bcftools mpileup -Ou -f {path_to_reference_file}  {path_to_folder}/*.bam | bcftools call --ploidy 1 -mv -Ob -o {path_to_folder}/{folder}.raw.bcf'
                    print('Executing:',command)
                    subprocess.call([command],shell=True)
                except Exception as e:
                    print(e)


def index_bcf(path_to_files, file_extension):
    '''Indexes bcf files'''
    for folder in os.listdir(path_to_files):
        path_to_folder=os.path.join(path_to_files,folder)
        if os.path.isdir(path_to_folder):
            files_in_folder=os.listdir(path_to_folder)
            for f in files_in_folder:
                #if file ends in file_extension 
                if re.search(r'{}$'.format(file_extension),f):
                    full_path_to_file = os.path.join(path_to_folder,f)
                    command=f'bcftools index {full_path_to_file}'
                    print('Executing:',command)
                    subprocess.call([command],shell=True)


def filter_bcf_by_species(path_to_files, file_extension):
    '''Finds the unique species/chromosomes from #CHROM column in a bcf file and outputs the filtered variants for that species to a new bcf file
    with the naming convention <species>.flt.bcf'''
    for folder in os.listdir(path_to_files):
        path_to_folder=os.path.join(path_to_files,folder)
        if os.path.isdir(path_to_folder):
            files_in_folder=os.listdir(path_to_folder)
            for f in files_in_folder:
                #if file ends in bcf (or specified file_extension)
                if re.search(r'{}$'.format(file_extension),f):
                    print('\n')
                    full_path_to_file = os.path.join(path_to_folder,f)
                    # Extract all unique chromosomes/species from file
                    command=f'bcftools view {full_path_to_file} |  grep -v "#" | cut -f 1 | uniq | sort'
                    unique_chrom=snp.save_process_output(command)
                    for chrom in unique_chrom:
                        species=chrom.split('_')
                        species_prefix = '_'.join([species[0],species[1],species[2]])
                        output_file_name = f'{folder}_{species_prefix}.flt.bcf'
                        full_path_to_output_file=os.path.join(path_to_folder,output_file_name)
                        filter_command=f'bcftools view {full_path_to_file} --regions {chrom} > {full_path_to_output_file}'
                        print('Executing:',filter_command)
                        subprocess.call([filter_command],shell=True)
                    
def get_number_of_variants(path_to_files,file_extension):
    # count the number of variants for each chromosome in bcf file and output to text file
    pass

def get_allele_frequency(path_to_files, file_extension, path_to_species_file):
    '''Uses vcftools to output the frequency of alleles for a chromosome in a vcf file.'''
    species_handle=open(path_to_species_file,'r')
    for species in path_to_species_file:
        for folder in os.listdir(path_to_files):
            path_to_folder=os.path.join(path_to_files,folder)
            if os.path.isdir(path_to_folder):
                output_file=f'{path_to_folder}/{folder}.raw.bcf'
                #if not os.path.exists(output_file):


    # Filter vcf file to extract species/create vcf file for each species



    # On each species file calculate allele frequency
    # Calculate mean depth of coverage per mouse sample 



    #vcftools --bcf mouse_12.raw.bcf --freq --chr Enterococcus_phage_A2_length_149431
    # vcftools --gzvcf input_file.vcf.gz --freq --chr 1 --out chr1_analysis

    #vcftools --gzvcf $SUBSET_VCF --depth --out $OUT

    # vcftools --gzvcf $SUBSET_VCF --depth --out $OUT


def main():
    #bcftools_mpileup('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','sorted.bam','/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/Combined.fasta')
    #index_bcf('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.bcf')
    filter_bcf_by_species('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.bcf')

if __name__ == '__main__':
    main()

