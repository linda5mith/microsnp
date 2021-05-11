import argparse
import sys
import os
import snptools as snp


def read_params(args):
    p = ap.ArgumentParser( description=
            "DESCRIPTION\n"
            " SNP caller: \n"
            " Tools for automating a snp pipeline using samtools and snpEff.\n\n"
            " Python v3.6 is required. Samtools, bcftools, vcftools, SNPEff\n\n"
            "\n------------------------------------------------------------------- \n\n"
             "",
            formatter_class=ap.RawTextHelpFormatter,
            add_help=False )
    arg = p.add_argument















parser = argparse.ArgumentParser()
parser.add_argument('--path_to_html_files', type=str, help='Path to html files for scraping.')
parser.add_argument('--file_extension', type=str, help='File extension of files e.g. \'html\'')
parser.add_argument('--path_to_output_files', type=str, default=os.getcwd(), help='Path to output scraped files.')

parser.add_argument('--path_to_txt_files', type=str, help='Path to html files for scraping.')

args = parser.parse_args()
#sys.stdout.write(str(scrape_snpeff_html_to_txt(args)))


def scrape_snpeff_html_to_txt(path_to_files, file_extension, path_to_output_files):
    files = snp.access_folder_contents(path_to_files,file_extension)
    for f in files:
        if 'Enterococcus_faecalis' in f: # Escherichia_coli
            print(f)
            # output_filename=snp.get_output_name(f)
            # print(output_filename)
            # soup = BeautifulSoup(open(f), "html.parser")
            # results = soup.find_all('table', attrs={'border': '0'})
            # output_file = open(f'{path_to_output_files}/{output_filename}.txt', "w")
            # for line in results:
            #     print('Outputting scrapes to:',output_file)
            #     output_file.write(line.text)
            # output_file.close()
            # results2 = soup.find_all('table', attrs={'border': '1'})
            # output_file = open(f'{path_to_output_files}/{output_filename}.txt', "a") # a for append mode
            # for line in results2:
            #     print('Outputting scrapes to:',output_file)
            #     output_file.write(line.text)

def main():
    scrape_snpeff_html_to_txt(args.path_to_html_files, args.file_extension, args.path_to_output_files)
    #scrape_snpeff_html_to_txt(/external_HDD4/linda/unc_mouse_trial/genomes/dev,'.html','/external_HDD4/linda/unc_mouse_trial/genomes/snpeff_reports')
    #scrape_snpeff_txt_to_csv('/external_HDD4/linda/unc_mouse_trial/genomes/snpeff_reports','.txt')
    #append_data('/external_HDD4/linda/unc_mouse_trial/genomes/snpeff_reports','csv')


if __name__ == '__main__':
    main()

