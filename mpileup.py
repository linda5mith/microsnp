import snptools as snp
import os
import subprocess
import re
import pandas as pd
import csv	

#  bcftools view var.raw.bcf.gz | bcftools filter -e'"Enterococcus"' 

# Count the number of variants by species and output to file 
# bcftools view var.raw.bcf.gz | bcftools filter -e'"Enterococcus"' | wc -l 


def bcftools_mpileup_single(path_to_files, file_extension, path_to_reference_file):
    '''Runs bcftools mpileup on a single bam file with the reference. Path to faidx indexed reference sequence file path_to_reference_file 
    needs to be included.'''
    files = snp.os_walk(path_to_files, file_extension)
    for f in files:
        # get mouseID
        subject_path=snp.get_file_dir(f)
        subject_clean=snp.get_file_basename(subject_path)
        # get sample_ID
        sample_id=snp.get_file_basename(f)
        sample_id_clean=sample_id.split('_')[0]
        command=f'bcftools mpileup -Ou -f {path_to_reference_file}  {f} | bcftools call --ploidy 1 -Ou -mv |  bcftools filter -s LowQual -e \'%QUAL<20\'  > {subject_path}/{subject_clean}_{sample_id_clean}.flt.vcf'
        print(command)
        subprocess.call([command],shell=True)

# bcftools view -e 'QUAL<20' mouse_1_UNC2FT2113.flt.vcf

def filter_vcf_qual(path_to_files, file_extension, qual=20):
    '''Filters out snp calls from vcf file with quality below specified qual. Default threshold is 20.'''
    files = snp.os_walk(path_to_files, file_extension)
    for f in files:
        # get mouse_ID
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



def bcftools_mpileup_multi(path_to_files, file_extension, path_to_reference_file):
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
    '''Finds the unique species/chromosomes from #CHROM column in a bcf/vcf file and outputs the filtered variants for that species to a new bcf file
    with the naming convention <species>.fltq.bcf'''
    files = snp.os_walk(path_to_files, file_extension)
    print(files)
    for f in files:
    # # Extract all unique chromosomes/species from bcf file
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
        print(f)
        output_file_name=snp.get_output_name(f)
        path_to_output_file=snp.get_file_dir(f)
        command =f'bcftools view -Oz -o {path_to_output_file}/{output_file_name}.vcf.gz {f}'
        try:
            print('Executing:',command)
            #subprocess.call([command],shell=True)
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
def filter_qual_vcf(path_to_files, file_extension, min_qual):
    files=snp.os_walk(path_to_files,file_extension)
    for f in files:
        filename = snp.get_output_name(f)
        path_to_output = snp.get_file_dir(f)
        full_filename=filename+'_flt'+'.vcf.gz'
        full_output_path = os.path.join(path_to_output,full_filename)
        try:
            # htsfile file.vcf.gz # prints file type
            # z for compressed vcf, b for compressed bcf
            command = f'bcftools view -i \'%QUAL>={min_qual}\' {f} -O z -o {full_output_path}'
            #command = f'vcftools --gzvcf {f} --minQ {min_qual} --recode --stdout | gzip -c > {full_output_path}'
            print('Executing:',command)
            subprocess.call([command],shell=True)
        except Exception as e:
            print(e)


# mv file.vcf.gz plain.vcf
# bcftools view -Oz -o compressed.vcf.gz plain.vcf
# htsfile compressed.vcf.gz
# bcftools index compressed.vcf.gz


# Generate gbk files using prokka
def prokka_annotate(path_to_references,file_extension):
    '''Annotates files in path file specified file_extension. Note: Sequence ID in FASTA file must be less than 24 characters.'''
    # source ~/anaconda3/etc/profile.d/conda.sh
    # conda activate prokka
    files=snp.access_folder_contents(path_to_references,file_extension)
    for f in files:
        prefix=snp.get_output_name(f)
        outdir=os.path.join(snp.get_file_dir(f),prefix)
        command = f'prokka {f} --outdir {outdir}/ --prefix {prefix}'
        try:
            subprocess.call([command],shell=True)
            print('Executing:',command)
        except Exception as e:
            print(e)

def clean_prokka(path_to_references,file_extension):
    '''Deletes lines 2-10 in prokka generated gbk file as these cause formatting issues.'''
    files=snp.access_subfolder_contents(path_to_references, file_extension)
    print(files)
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
        
# Build databases for references for SnpEff
def build_snpeff_db_gbk(path_to_snpeff_installation, file_extension):
    '''Note 0_0 SnpEff path to data folder is hardcoded so you have to build database inside snpEff program folder.
    Builds a SnpEff reference database for each gbk file with specified file_extension in path_to_references. 
    Note: must first add genomes to snpEff.config and gbk files to /data/programs/snpEff/data.
    '''
    files = snp.os_walk(path_to_snpeff_installation,file_extension)
    print(files)
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

def recompress_bgzip_vcfs(path_to_files, file_extension):
    '''This function is to counteract error you get if your vcf files are not properly compressed. E.g. overcome error:
    Failed to open file.vcf.gz: not compressed with bgzip. Converts files with file_extension vcf.gz back to vcf, recompresses and indexes them.'''
    rename_file_extension(path_to_files, file_extension)
    bcftools_compress_vcf(path_to_files, file_extension)
    bcftools_index_vcf(path_to_files, file_extension)
    htsfile(path_to_files, file_extension)

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
    e.g. +2 means find intersection between 2 files or more. Other options: output positions present in this many (=), this many or more (+), this many or fewer (-), or the exact same (~) files'''
    for folder in os.listdir(path_to_files):
        path_to_folder=os.path.join(path_to_files,folder)
        if os.path.isdir(path_to_folder):
            try:
                command = f'bcftools isec {path_to_folder}/*{file_extension} -p {os.path.join(path_to_folder,isec_outdir)} -n {nfiles_output_pos}'
                print(command, '\n')
                subprocess.call([command],shell=True)
            except Exception as e:
                print(e)


def main():
    # --------------------- Pipeline 2.0 14_05_21 ----------------------------------------------------------------------------------------------------------
    # --------------------- SNP calling --------------------------------------------------------------------------------------------------------------------
    #bcftools_mpileup_single('/external_HDD4/linda/unc_mouse_trial/snp_pipeline', '.sorted.bam', '/external_HDD4/linda/unc_mouse_trial/snp_pipeline/Combined.fasta')
    #filter_vcf_qual('/external_HDD4/linda/unc_mouse_trial/snp_pipeline', '.flt.vcf', 20)

    # Compress files with bgzip and then index before filtering can be applied
    #bcftools_compress_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltq.vcf')
    #bcftools_index_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltq.vcf.gz')
    #htsfile('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltq.vcf.gz')

    # Filter fltq.vcf.gz by species 
    # filter_bcf_by_species('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1', '.fltq.vcf.gz')

    # create species folders 
    # snp.create_species_folders('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1', '/external_HDD4/linda/unc_mouse_trial/genomes/species_sequences.txt')

    # move species files to corresponding folder
    # snp.move_files_to_species_folder('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/','.fltqs.vcf','/external_HDD4/linda/unc_mouse_trial/genomes/species_sequences.txt')


    # -------------------- 03_06_21 Find common variants between multiple files -----------------------------------------------------------------------------
    
    # Index filtered species files
    #bcftools_compress_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltqs.vcf')
    #bcftools_index_vcf('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/', '.fltq.vcf.gz')


    bcftools_isec('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1', '.fltq.vcf.gz','isec_out','+2')

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

