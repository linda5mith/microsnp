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
        







