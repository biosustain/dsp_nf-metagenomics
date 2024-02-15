#!/usr/bin/env python3
# This script needs two arguments: files to summarize list and output file name
# Example of running it: python files_list $assembly_files assembly_stats.txt

# Importing the assembly directory and the output file name
#import sys
import re # reg exp libraries
#input_args = sys.argv[1:]
import argparse

def main():
    # 2 parameters --assembly_files, --output
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--assembly_files", help="Assembly files", nargs='+', required=True)
    parser.add_argument("-o", "--output", help="Output file name", required=True)
    args = parser.parse_args()

    assemblies = args.assembly_files
    assemblies_sorted = sorted(assemblies)
    #assemblies_sorted = assemblies.sort()
    print(assemblies)
    print(type(assemblies))
    output_file = args.output
    
    with open(output_file, 'w') as f:
        f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('Sample', 'num_contigs', 'Total_contig_bp', 'min', 'max', 'Avg', 'N50'))
        for sample in assemblies:
            #Function that opens file and splits data into lines. 
            data_file = open_file(sample)
            print(data_file)
            summary = [sample]
            for element in re.split(",", data_file):
                print(element)
                summary += re.findall("([\d]+) contigs", element)
                summary += re.findall("total ([\d]+) bp", element)
                summary += re.findall("min ([\d]+) bp", element)
                summary += re.findall("max ([\d]+) bp", element)
                summary += re.findall("avg ([\d]+) bp", element)
                summary += re.findall("N50 ([\d]+) bp", element)
            print(summary)
            as_strings = [str(x) for x in summary]
            final_string = '\t'.join(as_strings)
            final_string += '\n'
            f.write(final_string)


def open_file(sample_name) -> list[str]:
    #Prints sample name
    print(sample_name)
    #Defines file location depending on sample name and opens file
    #file_location = input_args[0]+"/"+sample_name+"/log"
    #file = open(file_location)
    #splits the data into lines
    with open(sample_name, 'r') as infile:
        data = infile.read().split('\n')
        #data = infile.readline(-3)
    return data[-3]
	              
main ()