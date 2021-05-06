import snptools as snp
import re
from bs4 import BeautifulSoup


path_to_output='/external_HDD4/linda/unc_mouse_trial/genomes/snpeff_reports'

# files = snp.access_folder_contents('/external_HDD4/linda/unc_mouse_trial/html_snpeff','.html')
# for f in files:
#     if 'Escherichia_coli' in f:
#         print(f)
#         output_filename=snp.get_output_name(f)
#         print(output_filename)
#         soup = BeautifulSoup(open(f), "html.parser")
#         results = soup.find_all('table', attrs={'border': '0'})
#         output_file = open(f'{path_to_output}/{output_filename}.txt', "w")
#         for t in test:
#             output_file.write(t.text)
#         output_file.close()


# Make csv file from data
data_dict={}
with open('/external_HDD4/linda/unc_mouse_trial/genomes/snpeff_reports/mouse_1_NC_011993.txt','r') as f:
    for line in f:
        line=line.strip()
        if re.match(r"^Genome$",line):
            genome=next(f)
            genome=genome.strip()
            print(genome)
            data_dict['Genome']=genome
        elif re.match(r"^Number of lines (input file)$",line):
            numVariants=next(f)
            numVariants=numVariants.strip()
            print(numVariants)
            data_dict['numVariants']=numVariants
        elif re.match(r"^Number of effects$",line):
            numEffects=next(f)
            numEffects=numEffects.strip()
            print(numEffects)
            data_dict['numEffects']=numEffects
        elif re.match(r"^Genome total length$",line):
            genome_length=next(f)
            genome_length=genome_length.strip()
            print(genome_length)
            data_dict['genome_length']=genome_length
        elif re.match(r"^Variant rate$",line):
            variant_rate=next(f)
            variant_rate=variant_rate.strip()
            print(variant_rate)
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
        
print(data_dict)        







