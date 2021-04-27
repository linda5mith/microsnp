import snptools as snp
import os
import subprocess
import re
import pandas as pd
import csv	

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

                    
def get_number_of_variants(path_to_files, file_extension, path_to_output_file=os.getcwd()):
    '''Counts the number of variants for each chromosome in all bcf files or with specified file_extension
     in specified path and outputs to csv file'''
    d = {}
    for folder in os.listdir(path_to_files):
        path_to_folder=os.path.join(path_to_files,folder)
        if os.path.isdir(path_to_folder):
            files_in_folder=os.listdir(path_to_folder)
            for f in files_in_folder:
                #if file ends in bcf (or specified file_extension)
                if re.search(r'{}$'.format(file_extension),f):
                    full_path_to_file = os.path.join(path_to_folder,f)
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
    '''Searches for all bcf files in path with given file_extension and converts them to vcf'''
    files=snp.access_subfolder_contents(path_to_files, file_extension)
    for f in files:
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command =f'bcftools view -Oz -o {path_to_output_file}/{output_file_name}.vcf.gz {f}'
        try:
            subprocess.call([command],shell=True)
            print('Executing:',command)
        except Exception as e:
            print(e)

def get_allele_freq(path_to_files, file_extension):
    '''Uses vcftools to output the frequency of alleles for a chromosome in a vcf file.'''
    #vcftools --gzvcf $SUBSET_VCF --freq2 --out $OUT --max-alleles 2
    files=snp.access_subfolder_contents(path_to_files, file_extension)
    for f in files:
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command=f'vcftools --gzvcf {f} --freq2 --out {path_to_output_file}/{output_file_name} --max-alleles 2'
        try:
            print('Executing:',command,'\n')
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def get_depth(path_to_files, file_extension):
    ''''Calculates the mean depth of coverage per individual.'''
    #vcftools --gzvcf $SUBSET_VCF --depth --out $OUT
    files=snp.access_subfolder_contents(path_to_files, file_extension)
    for f in files:
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command=f'vcftools --gzvcf {f} --depth --out {path_to_output_file}/{output_file_name}_depth_indiv'
        try:
            print('Executing:',command,'\n')
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

def get_site_mean_depth(path_to_files, file_extension):
    '''Estimates the mean depth of coverage for each site.'''
    #vcftools --gzvcf $SUBSET_VCF --site-mean-depth --out $OUT
    files=snp.access_subfolder_contents(path_to_files, file_extension)
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
    files=snp.access_subfolder_contents(path_to_files, file_extension)
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
    files=snp.access_subfolder_contents(path_to_files, file_extension)
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
    files=snp.access_subfolder_contents(path_to_files, file_extension)
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
    files=snp.access_subfolder_contents(path_to_files, file_extension)
    for f in files:
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command=f'vcftools --gzvcf {f} --het --out {path_to_output_file}/{output_file_name}_het'
        try:
            print('Executing:',command,'\n')
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)

# Filter variants from bcf files which have QUAL <=30 (1 in 1000 chance of being incorrect)


# Filter variants according to parameters present in https://speciationgenomics.github.io/filtering_vcfs/ 

# Generate gbk files using prokka
def prokka_annotate(path_to_references,file_extension):
    '''Annotates files in path file specified file_extension'''
    # source ~/anaconda3/etc/profile.d/conda.sh
    # conda activate prokka
    files=snp.access_folder_contents(path_to_references,file_extension)
    for f in files:
        prefix=snp.get_output_name(f)
        outdir=snp.get_file_dir(f)
        command = f'prokka {f} --outdir {outdir} --prefix {prefix}'
        print(command)


# Build databases for references for SnpEff
def build_snpeff_db(path_to_references, file_extension, path_to_snpeff_installation):
    '''Builds a SnpEff reference database for each file with specified file_extension in path_to_references'''
    files=snp.access_folder_contents(path_to_references,file_extension)
    for f in files:
        try:
            command = f'java -jar {path_to_snpeff_installation}/snpEff.jar build {f} -v'
            #java -Xmx8g -jar snpEff.jar download -c path/to/snpEff/snpEff.config -v GRCh37.75
            print('Executing:',command)
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)


def main():
    # --------------------- SNP calling --------------------------------------------------------------------------------------------------------------------
    #bcftools_mpileup('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','sorted.bam','/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/Combined.fasta')
    #index_bcf('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.bcf')
    #filter_bcf_by_species('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.bcf')
    #get_number_of_variants('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','flt.bcf','/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2')

    # --------------------- Create folders for each species inside each mouse folder -----------------------------------------------------------------------
    #snp.create_species_folders('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','/external_HDD4/linda/unc_mouse_trial/genomes/species_sequences.txt')
    #snp.move_files_to_species_folder('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.flt.bcf','/external_HDD4/linda/unc_mouse_trial/genomes/species_sequences.txt')

    # --------------------- Convert bcf to vcf  ------------------------------------------------------------------------------------------------------------
    #bcf_to_vcf('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.flt.bcf') 

    # --------------------- Vcftools statistics  -----------------------------------------------------------------------------------------------------------
    # get_allele_freq('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')
    # get_depth('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')
    # get_site_mean_depth('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')
    # get_site_quality('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')
    # get_missing_prop_per_site('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')
    # get_missing_prop_per_indiv('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')
    # get_het('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2','.vcf.gz')

    # --------------------- Annotate using Prokka ----------------------------------------------------------------------------------------------------------
    prokka_annotate('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/reference_genomes_db','.fasta')


    # --------------------- SnpEff pipeline ----------------------------------------------------------------------------------------------------------------
    #build_snpeff_db('/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/reference_genomes_db','.fasta','/data/programs/snpEff')


if __name__ == '__main__':
    main()

