import snptools as snp
import re
import pandas as pd
from bs4 import BeautifulSoup

import snptools as snp


def scrape_snpeff_html_to_txt(path_to_files, file_extension, species, path_to_output_files):
    files = snp.os_walk(path_to_files,file_extension)
    for f in files:
        if species in f:
            print('Scraping:',f)
            output_filename=snp.get_output_name(f)
            soup = BeautifulSoup(open(f), "html.parser")
            results = soup.find_all('table', attrs={'border': '0'})
            output_file = open(f'{path_to_output_files}/{output_filename}.txt', "w")
            for line in results:
                print('Outputting scrapes to:',output_file)
                output_file.write(line.text)
            output_file.close()
            results2 = soup.find_all('table', attrs={'border': '1'})
            output_file = open(f'{path_to_output_files}/{output_filename}.txt', "a") # a for append mode
            for line in results2:
                print('Outputting scrapes to:',output_file)
                output_file.write(line.text)
            print('\n')  


def scrape_snpeff_txt_to_csv_2(path_to_files,file_extension):
    data_dict={}
    files=snp.access_folder_contents(path_to_files, file_extension)
    for fle in files:
        with open(fle,'r') as f:
            for line in f:
                print(line)
                line=line.strip()
                if re.match(r"^Command line arguments$",line):
                    command=next(f)
                    command=command.strip()
                    mouse=command.split('/')[6]
                    data_dict['mouse']=mouse
                elif re.match(r"^Genome$",line):
                    genome=next(f)
                    genome=genome.strip()
                    data_dict['Genome']=genome
                elif re.match(r"^Number of variants",line):
                    numVariants=next(f)
                    numVariants=numVariants.strip()
                    data_dict['numVariants']=numVariants
                elif re.match(r"^Number of known variants",line):
                    skip=next(f)
                    numknownVariants=next(f)
                    numknownVariants=numknownVariants.strip()
                    data_dict['numKnownVariants']=numknownVariants
                elif re.match(r"^SNP$",line):
                    print(line)
                    #skip=next(f)
                    snps=next(f)
                    snps=snps.strip()
                    data_dict['SNP']=snps
    print(data_dict)


# Make csv file from data
def scrape_snpeff_txt_to_csv(path_to_files,file_extension):
    data_dict={}
    files=snp.access_folder_contents(path_to_files, file_extension)
    for fle in files:
        with open(fle,'r') as f:
            for line in f:
                line=line.strip()
                if re.match(r"^Command line arguments$",line):
                    command=next(f)
                    command=command.strip()
                    mouse=command.split('/')[6]
                    data_dict['mouse']=mouse
                elif re.match(r"^Genome$",line):
                    genome=next(f)
                    genome=genome.strip()
                    data_dict['Genome']=genome
                elif re.match(r"^Number of variants",line):
                    numVariants=next(f)
                    numVariants=numVariants.strip()
                    data_dict['numVariants']=numVariants
                elif re.match(r"^Number of known variants",line):
                    skip=next(f)
                    numknownVariants=next(f)
                    numknownVariants=numknownVariants.strip()
                    data_dict['numKnownVariants']=numknownVariants
                elif re.match(r"^SNP",line):
                    snps=next(f)
                    snps=snps.strip()
                    data_dict['SNP']=snps
                elif re.match(r"^MNP",line):
                    mnp=next(f)
                    mnp=mnp.strip()
                    data_dict['MNP']=mnp
                elif re.match(r"^INS",line):
                    ins=next(f)
                    ins=ins.strip()
                    data_dict['INS']=ins
                elif re.match(r"^DEL",line):
                    delete=next(f)
                    delete=delete.strip()
                    data_dict['DEL']=delete
                elif re.match(r"^MIXED",line):
                    mixed=next(f)
                    mixed=mixed.strip()
                    data_dict['MIXED']=mixed
                elif re.match(r"^INV",line):
                    inv=next(f)
                    inv=inv.strip()
                    data_dict['INV']=inv
                elif re.match(r"^DUP",line):
                    dup=next(f)
                    dup=dup.strip()
                    data_dict['DUP']=dup
                elif re.match(r"^BND",line):
                    bnd=next(f)
                    bnd=bnd.strip()
                    data_dict['BND']=bnd
                elif re.match(r"^Number of effects$",line):
                    numEffects=next(f)
                    numEffects=numEffects.strip()
                    #print(numEffects)
                    data_dict['numEffects']=numEffects
                elif re.match(r"^Genome total length$",line):
                    genome_length=next(f)
                    genome_length=genome_length.strip()
                    #print(genome_length)
                    data_dict['genome_length']=genome_length
                elif re.match(r"^Variant rate$",line):
                    variant_rate=next(f)
                    variant_rate=variant_rate.strip()
                    #print(variant_rate)
                    data_dict['variant_rate']=variant_rate
                elif re.match(r"^HIGH$",line):
                    skip=next(f)
                    high_count=next(f)
                    high_count=high_count.strip()
                    high_perc=next(f)
                    high_perc=high_perc.strip()
                    data_dict['high_count']=high_count
                    data_dict['high_%']=high_perc
                elif re.match(r"^LOW$",line):
                    skip=next(f)
                    low_count=next(f)
                    low_count=low_count.strip()
                    low_perc=next(f)
                    low_perc=low_perc.strip()
                    data_dict['low_count']=low_count
                    data_dict['low_%']=low_perc
                elif re.match(r"^MODERATE$",line):
                    skip=next(f)
                    mod_count=next(f)
                    mod_count=mod_count.strip()
                    mod_perc=next(f)
                    mod_perc=mod_perc.strip()
                    data_dict['mod_count']=mod_count
                    data_dict['mod_%']=mod_perc
                elif re.match(r"^MODIFIER$",line):
                    skip=next(f)
                    modif_count=next(f)
                    modif_count=modif_count.strip()
                    modif_perc=next(f)
                    modif_perc=modif_perc.strip()
                    data_dict['modif_count']=modif_count
                    data_dict['modif_%']=modif_perc
                elif re.match(r"^MISSENSE$",line):
                    skip=next(f)
                    missense_count=next(f)
                    missense_count=missense_count.strip()
                    missense_perc=next(f)
                    missense_perc=missense_perc.strip()
                    data_dict['missense_count']=missense_count
                    data_dict['missense_%']=missense_perc
                elif re.match(r"^NONSENSE$",line):
                    skip=next(f)
                    nonsense_count=next(f)
                    nonsense_count=nonsense_count.strip()
                    nonsense_perc=next(f)
                    nonsense_perc=nonsense_perc.strip()
                    data_dict['nonsense_count']=nonsense_count
                    data_dict['nonsense_%']=nonsense_perc
                elif re.match(r"^SILENT$",line):
                    skip=next(f)
                    silent_count=next(f)
                    silent_count=silent_count.strip()
                    silent_perc=next(f)
                    silent_perc=silent_perc.strip()
                    data_dict['silent_count']=silent_count
                    data_dict['silent_%']=silent_perc
                elif re.match(r"^conservative_inframe_insertion$",line):
                    skip=next(f)
                    conservative_inframe_insertion_count=next(f)
                    conservative_inframe_insertion_count=conservative_inframe_insertion_count.strip()
                    conservative_inframe_insertion_perc=next(f)
                    conservative_inframe_insertion_perc=silent_perc.strip()
                    data_dict['conservative_inframe_insertion_count']=conservative_inframe_insertion_count
                    data_dict['conservative_inframe_insertion_%']=conservative_inframe_insertion_perc
                elif re.match(r"^downstream_gene_variant$",line):
                    skip=next(f)
                    downstream_gene_variant_count=next(f)
                    downstream_gene_variant_count=downstream_gene_variant_count.strip()
                    downstream_gene_variant_perc=next(f)
                    downstream_gene_variant_perc=downstream_gene_variant_perc.strip()
                    data_dict['downstream_gene_variant_count']=downstream_gene_variant_count
                    data_dict['downstream_gene_variant_%']=downstream_gene_variant_perc
                elif re.match(r"^frameshift_variant$",line):
                    skip=next(f)
                    frameshift_variant_count=next(f)
                    frameshift_variant_count=frameshift_variant_count.strip()
                    frameshift_variant_perc=next(f)
                    frameshift_variant_perc=frameshift_variant_perc.strip()
                    data_dict['frameshift_variant_count']=frameshift_variant_count
                    data_dict['frameshift_variant_%']=frameshift_variant_perc
                elif re.match(r"^intergenic_region$",line):
                    skip=next(f)
                    intergenic_region_count=next(f)
                    intergenic_region_count=intergenic_region_count.strip()
                    intergenic_region_perc=next(f)
                    intergenic_region_perc=intergenic_region_perc.strip()
                    data_dict['intergenic_region_count']=intergenic_region_count
                    data_dict['intergenic_region_%']=intergenic_region_perc
                elif re.match(r"^missense_variant$",line):
                    skip=next(f)
                    missense_variant_count=next(f)
                    missense_variant_count=missense_variant_count.strip()
                    missense_variant_perc=next(f)
                    missense_variant_perc=missense_variant_perc.strip()
                    data_dict['missense_variant_count']=missense_variant_count
                    data_dict['missense_variant_%']=missense_variant_perc
                elif re.match(r"^splice_region_variant$",line):
                    skip=next(f)
                    splice_region_variant_count=next(f)
                    splice_region_variant_count=splice_region_variant_count.strip()
                    splice_region_variant_perc=next(f)
                    splice_region_variant_perc=splice_region_variant_perc.strip()
                    data_dict['splice_region_variant_count']=splice_region_variant_count
                    data_dict['splice_region_variant_%']=splice_region_variant_perc
                elif re.match(r"^start_lost$",line):
                    skip=next(f)
                    start_lost_count=next(f)
                    start_lost_count=start_lost_count.strip()
                    start_lost_perc=next(f)
                    start_lost_perc=start_lost_perc.strip()
                    data_dict['start_lost_count']=start_lost_count
                    data_dict['start_lost_%']=start_lost_perc
                elif re.match(r"^stop_gained$",line):
                    skip=next(f)
                    stop_gained_count=next(f)
                    stop_gained_count=stop_gained_count.strip()
                    stop_gained_perc=next(f)
                    stop_gained_perc=stop_gained_perc.strip()
                    data_dict['stop_gained_count']=stop_gained_count
                    data_dict['stop_gained_%']=stop_gained_perc
                elif re.match(r"^stop_lost$",line):
                    skip=next(f)
                    stop_lost_count=next(f)
                    stop_lost_count=stop_lost_count.strip()
                    stop_lost_perc=next(f)
                    stop_lost_perc=stop_lost_perc.strip()
                    data_dict['stop_lost_count']=stop_lost_count
                    data_dict['stop_lost_%']=stop_lost_perc
                elif re.match(r"^synonymous_variant$",line):
                    skip=next(f)
                    synonymous_variant_count=next(f)
                    synonymous_variant_count=synonymous_variant_count.strip()
                    synonymous_variant_perc=next(f)
                    synonymous_variant_perc=synonymous_variant_perc.strip()
                    data_dict['synonymous_variant_count']=synonymous_variant_count
                    data_dict['synonymous_variant_%']=synonymous_variant_perc
                elif re.match(r"^upstream_gene_variant$",line):
                    skip=next(f)
                    upstream_gene_variant_count=next(f)
                    upstream_gene_variant_count=upstream_gene_variant_count.strip()
                    upstream_gene_variant_perc=next(f)
                    upstream_gene_variant_perc=upstream_gene_variant_perc.strip()
                    data_dict['upstream_gene_variant_count']=upstream_gene_variant_count
                    data_dict['upstream_gene_variant_%']=upstream_gene_variant_perc
                elif re.match(r"^DOWNSTREAM$",line):
                    skip=next(f)
                    DOWNSTREAM_count=next(f)
                    DOWNSTREAM_count=DOWNSTREAM_count.strip()
                    DOWNSTREAM_perc=next(f)
                    DOWNSTREAM_perc=DOWNSTREAM_perc.strip()
                    data_dict['DOWNSTREAM_count']=DOWNSTREAM_count
                    data_dict['DOWNSTREAM_%']=DOWNSTREAM_perc
                elif re.match(r"^EXON$",line):
                    skip=next(f)
                    EXON_count=next(f)
                    EXON_count=EXON_count.strip()
                    EXON_perc=next(f)
                    EXON_perc=EXON_perc.strip()
                    data_dict['EXON_count']=EXON_count
                    data_dict['EXON_%']=EXON_perc
                elif re.match(r"^INTERGENIC$",line):
                    skip=next(f)
                    INTERGENIC_count=next(f)
                    INTERGENIC_count=INTERGENIC_count.strip()
                    INTERGENIC_perc=next(f)
                    INTERGENIC_perc=INTERGENIC_perc.strip()
                    data_dict['INTERGENIC_count']=INTERGENIC_count
                    data_dict['INTERGENIC_%']=INTERGENIC_perc
                elif re.match(r"^UPSTREAM$",line):
                    skip=next(f)
                    UPSTREAM_count=next(f)
                    UPSTREAM_count=UPSTREAM_count.strip()
                    UPSTREAM_perc=next(f)
                    UPSTREAM_perc=UPSTREAM_perc.strip()
                    data_dict['UPSTREAM_count']=UPSTREAM_count
                    data_dict['UPSTREAM_%']=UPSTREAM_perc
                
            print(data_dict)        

            # Make csv of each file
            output_name=snp.get_output_name(fle)
            output_dir=snp.get_file_dir(fle)
            # print(output_name)
            # print(output_dir)
            print('\n')
            df = pd.DataFrame.from_dict(data_dict, orient='index')
            df = df.transpose()
            df.to_csv(f'{output_dir}/{output_name}.csv', index=False)


# Append/concat multiple csv file contents together
def append_data(path_to_dfs,file_extension):
    files=snp.access_folder_contents(path_to_dfs, file_extension)
    mydf = pd.DataFrame()
    for f in files:
        df=pd.read_csv(f)
        mydf=mydf.append(df)
    mydf.to_csv('NC_017316.1_Enterococcus_faecalis_OG1RF_all_snps.csv')
    return mydf


def main():
    #scrape_snpeff_html_to_txt('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/','.html','NC_017316','/external_HDD4/linda/unc_mouse_trial/snp_pipeline/html_reports')
    #scrape_snpeff_txt_to_csv('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/html_reports','.txt')
    append_data('/external_HDD4/linda/unc_mouse_trial/snp_pipeline/html_reports','csv')


if __name__ == '__main__':
    main()
   

# NC_017316.1_Enterococcus_faecalis_OG1RF_complete_genome_length_2739625
# NC_011993.1_Escherichia_coli_LF82_complete_genome_length_4773108
