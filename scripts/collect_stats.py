
"""Append results from pipeline to results.csv

Author: Jimmy Capecci
"""

# import libraries
from pathlib import Path
import argparse

def main():

    # set up  parser for command line support
    parser = argparse.ArgumentParser(description='Adds bench marking data to csv')

    # add arguments for command 
    parser.add_argument('-i', '--input', required=True, help='Input directory where benchmarking stats is located')
    parser.add_argument('-r', '--chromosome', required=True, help='Chromosome being simulated')
    parser.add_argument('-c', '--coverage', required=True, help='Coverage used to simulate reads')
    parser.add_argument('-s', '--seed', required=True, help='Seed used to simulate reads')

    # grab arguments from command line
    args = parser.parse_args()

    # set variables to command line arguments
    data_path = args.input
    chromosome = args.chromosome
    coverage = args.coverage
    seed = args.seed

    # set output csv path and dictionary for collection
    output_csv = 'results.csv'
    data_collector = {}

    # if results.csv doesnt exist intialize it with headers
    if not Path(output_csv).exists():
        headers = ['Chromosome', 'Coverage', 'Seed', 'Tool', 'Time', 
               'Memory', 'True-pos-baseline', 'True-pos-call',
               'False-pos', 'False-neg', 'Precision', 'Sensitivity',
               'F-measure']
        with open(output_csv, 'w') as o:
            o.write(','.join(headers) + '\n')
            

    # grab txt files that was output from /usr/bin/time
    path = Path(data_path)
    time_paths = list(path.glob('*.txt'))

    # string to match time and ram
    time_indicator = 'Elapsed (wall clock) time (h:mm:ss or m:ss):'
    ram_indicator = 'Maximum resident set size (kbytes):'

    # for each text file 
    for p in time_paths:
        p_name = str(p).split('.')[0].split('/')[-1] # grab tool used 
        with open(p, 'r') as f:
            for line in f:
                if time_indicator in line:
                    t = line.split(' ')[-1].strip() # grab time

                elif ram_indicator in line:
                    r = line.split(' ')[-1].strip() # grab ram

        data_collector[p_name] = {'time': t, 'ram': r} # assign tool to time and ram

    # grab text files associated with accuracy
    data_path += '/accuracy'
    path = Path(data_path)
    acc_paths = list(path.glob('*.txt'))

    # for each text file
    for a in acc_paths:
        p_name = str(a).split('/')[-1].split('_')[0] # grab tool used
        with open(a, 'r') as f:
            f.readline()
            f.readline()
            line = f.readline().split() # grab statistics 

        # assign tool to each statistic
        data_collector[p_name]['True-pos-baseline'] = line[1]
        data_collector[p_name]['True-pos-call'] = line[2]
        data_collector[p_name]['False-pos'] = line[3]
        data_collector[p_name]['False-neg'] = line[4]
        data_collector[p_name]['Precision'] = line[5]
        data_collector[p_name]['Sensitivity'] = line[6]
        data_collector[p_name]['F-measure'] = line[-1]

    # write results to csv 
    with open('results.csv', 'a') as o:
        for key, stats_dict in data_collector.items():
            row = [chromosome, coverage, seed, key, *stats_dict.values()]
            o.write(','.join(row) + '\n')

    
if __name__ == '__main__':
    main()
