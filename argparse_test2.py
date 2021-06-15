import argparse

# Same main parser as usual
parser = argparse.ArgumentParser()

# Usual arguments which are applicable for the whole script / top-level args
parser.add_argument('--verbose', help='Common top-level parameter',
                    action='store_true', required=False)

# Same subparsers as usual
subparsers = parser.add_subparsers(help='Desired action to perform', dest='action')

# Usual subparsers not using common options
parser_other = subparsers.add_parser("extra-action", help='Do something without db')

# Create parent subparser. Note `add_help=False` and creation via `argparse.`
parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument('-p', help='add db parameter', required=True)

# Subparsers based on parent

parser_create = subparsers.add_parser("create", parents=[parent_parser],
                                      help='Create something')
# Add some arguments exclusively for parser_create

parser_update = subparsers.add_parser("update", parents=[parent_parser],
                                      help='Update something')
# Add some arguments exclusively for parser_update 


def read_params(args):
    p = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''    ---------------------------------------------------------------------------------------------------------------------------------------------------------------

            1. The first step is to create a mpileup of genomes vs. the reference genome. To make a mpileup consisting of one sample file vs. reference genome use:
            python mpileup.py bcftools_mpileup_single /external_HDD4/linda/unc_mouse_trial/snp_pipeline .sorted.bam /external_HDD4/linda/unc_mouse_trial/snp_pipeline/Combined.fasta

            -------------------------------------------------------------------'''),
            epilog='Linda Smith https://github/flannsmith/snp-calling')

    f = p.add_argument_group('Function options')
    arg = f.add_argument
    FUNCTION_MAP = {'bcftools_mpileup_single' : bcftools_mpileup_single}

    arg( '--command', choices=FUNCTION_MAP.keys(), help = "Input which type of file to process.")
    func = FUNCTION_MAP[command]
    func()
    
    arg = p.add_argument
    arg('inp', metavar='INPUT_DIR', type=str, nargs='?', default=None, help='Path to files to process. If a parent directory specified, specify the \
    file extenstion by --input_type are recursively retrieved from all subdirectories in given path.')

    arg('output', metavar='OUTPUT_DIR', type=str, nargs='?', default=os.getcwd(),
         help= "Path to output files. Default is current directory.")

    g = p.add_argument_group('Required arguments')
    arg = g.add_argument
    input_type_choices = ['vcf','vcf.gz','flt.vcf.gz','fltq.vcf.gz']
    arg( '--input_type', choices=input_type_choices, help = "Input which type of file to process.")

    return vars(p.parse_args())
