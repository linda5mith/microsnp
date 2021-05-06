import snptools as snp
from bs4 import BeautifulSoup


path_to_output='/external_HDD4/linda/unc_mouse_trial/genomes/snpeff_reports'

files = snp.access_folder_contents('/external_HDD4/linda/unc_mouse_trial/html_snpeff','.html')
for f in files:
    if 'Escherichia_coli' in f:
        print(f)
        output_filename=snp.get_output_name(f)
        print(output_filename)
        soup = BeautifulSoup(open(f), "html.parser")
        results = soup.find_all('table', attrs={'border': '0'})
        output_file = open(f'{path_to_output}/{output_filename}.txt', "w")
        for t in test:
            output_file.write(t.text)
        output_file.close()


# Make csv file from data





