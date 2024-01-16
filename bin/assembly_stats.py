#!/usr/bin/env python3
# This script needs two arguments: assembly results folder and output file name
# Example of running it: python assembly_stats.py assembly assembly_stats.txt

# Importing the assembly directory and the output file name
import sys
import re # reg exp libraries
input_args = sys.argv[1:]

def main():
    from os import listdir
    samples=listdir(input_args[0]) #"/work3/apca/orange_peel/02_assembly/assembly")
    sorted_samples = sorted(samples)
                    
    with open(input_args[1], 'w') as f:
        f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('Sample', 'num_contigs', 'Total_contig_bp', 'min', 'max', 'Avg', 'N50'))
        for sample in sorted_samples:
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
    file_location = input_args[0]+"/"+sample_name+"/log"
    file = open(file_location)
    #splits the data into lines
    data = file.read().split("\n")
    return data[-3]
	                
main ()