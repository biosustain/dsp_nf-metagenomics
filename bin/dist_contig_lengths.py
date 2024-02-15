#!/usr/bin/env python3
# This 
# This script needs two arguments: assembly files and output file name
# Example of running it: python dist_contig_lengths.py $assembly_files contigs_summary.txt

# Importing the assembly directory and the output file name
import argparse
import re

def main():
   
    # 2 parameters --assembly_files, --output
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--assembly_files", help="Assembly files", nargs='+', required=True)
    parser.add_argument("-o", "--output", help="Output file name", required=True)
    args = parser.parse_args()

    assemblies = args.assembly_files
    print(assemblies)
    print(type(assemblies))
    assemblies_sorted = sorted(assemblies, reverse=True)
    #assemblies_sorted = assemblies.sort()
    output_file = args.output
    

    with open(output_file, 'w') as f:
        f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % ('Sample', 'over_20000_length', 'max', 'Total_contig_length_over_20000', 'Total_contig_length', 'Difference_Total_contig_length-Total_contig_length_over_20000'))
        for sample in assemblies:
            #Function that opens file and splits data into lines. 
            data_file = open_file(sample)
            data_processed = analyse_data(data_file)
            # we add sample name to the contigs data
            in_one_line = [sample]
            in_one_line += data_processed
            as_strings = [str(x) for x in in_one_line]
            final_string = '\t'.join(as_strings)
            final_string += '\n'
            f.write(final_string)


def open_file(sample_name) -> list[str]:
    #Prints sample name
    print(sample_name)
    #Defines file location depending on sample name and opens file
    #file_location = input_args[0]+"/"+sample_name+"/final.contigs.fa"
    #file = open(file_location)
    #splits the data into lines
    with open(sample_name, 'r') as infile:
        data = infile.read().split('\n')
    return data

def analyse_data(data) -> list:
    lengths_graph = []
    lengths_real = []
    lengths_20000 = []
    k=0
    sum = 0
    for line in data:
        # we add this in case that there are empty lines, otherwise line[0] does not work
        if line == '':
                continue
        k=k+1
        if line[0] == ">":
            line_split = line.split()
            
            lengths_real.append(int(line_split[3][4:]))
            if int(line_split[3][4:]) > 20000:
                lengths_graph.append(20000)
                lengths_20000.append(int(line_split[3][4:]))
            else:
                lengths_graph.append(int(line_split[3][4:]))
                
    sumk, sumi = 0,0
    for i in lengths_20000:
        sumk = sumk + i
        for k in lengths_real: 
            sumi = sumi + k
    data_processed = [str(len(lengths_20000)), str(max(lengths_20000)/20000), str(sumk), str(sumi), str(sumi-sumk)]
    return(data_processed)
	                
main ()