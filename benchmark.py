""" 
command line wrapper to run snakefile

Author: Jimmy Capecci
"""

# import libraries
import argparse
import subprocess

def main():

    # set up parser with arguements
    parser = argparse.ArgumentParser(description='Runs benchmarking pipeline')
    parser.add_argument('-c', '--coverage', help='Sequencing Coverage', required=False,  default='30')
    parser.add_argument('-s', '--seed', help='random seed to simulate reeds', required=False, default='42')
    parser.add_argument('-r', '--chromosome', help='Chromosome to download', required=False, default='chr20')

    # parse arguements
    args = parser.parse_args()

    # set output path
    output = f'coverage{args.coverage}_seed{args.seed}'

    # run command
    subprocess.run([
        'snakemake',
        '--cores', '4',
        '--config',
        f'coverage={args.coverage}',
        f'seed={args.seed}',
        f'chromosome={args.chromosome}',
        f'folder={output}',
    ], check=True)


if __name__ == '__main__':
    main()
    