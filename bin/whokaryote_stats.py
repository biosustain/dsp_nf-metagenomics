#!/usr/bin/env python3
# This is a script to calculate some Whokaryote statistics
# This script needs two arguments: whokaryote files and output file name
# Example of running it: python whokaryote_stats.py $predictions whokaryote_stats.txt

# Importing python libraries
import argparse
import re # reg exp libraries

def main():

# 2 parameters --predictions, --output
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--predictions", help="Whokaryote prediction files", nargs='+', required=True)
    parser.add_argument("-o", "--output", help="Output file name", required=True)
    args = parser.parse_args()

    predictions = args.predictions
    print(predictions)
    print(type(predictions))
    predictions_sorted = sorted(predictions)
    output_file = args.output
	
	#Loops through all samples 
    with open(output_file, 'w') as f:
        f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('Sample', 'Total_length_bp_Euk', 'num_contigs_Euk', 'Total_length_bp_Prok', 'num_contigs_Prok', 'Total_length_bp_Arc', 'num_contigs_Arc'))
        for file in predictions_sorted:
            #Function that opens file and splits data into lines. 
            data = open_file(file)
            #Function that loops through data and counts number of contigs and length for eukaryotes , prokaryotes and archeas.
            counts = analyse_data(data)
            print(counts)
            # We get the sample name from the file name
            sample = re.findall("([\w]+)\.featuretable_predictions_T.tsv", file)[0]
            #file.split(".")
            # we add sample name to the data
            in_one_line = [str(sample)]
			# We add all the values to one line
            for key, value in counts.items():
                in_one_line += value
			# Preparing the data as string to write in a file
            as_strings = [str(x) for x in in_one_line]
            final_string = '\t'.join(as_strings)
            final_string += '\n'
            f.write(final_string)


def open_file(sample_name) -> list:
    #Prints sample name
    print(sample_name)
    #Opens the file and reads the data lines
    #splits the data into lines
    with open(sample_name, 'r') as infile:
        data = infile.read().split('\n')
        #data = infile.readlines() # must be better but does not work for now
    return data
    
def analyse_data(data) -> dict:
    #Creates and empty dictionary
    kingdom_dict = {
        "eukaryote" : [0,0], 
        "prokaryote": [0,0], 
        "archea" : [0,0]
    }
    #Loops through each line of the data 
    for line in data:
        #Splits the line by tab, as the data is from a .tsv file 
        line_split = line.split("\t")
        #Checks kingdom of each line, and counts number of contigs for
        # each as well as adding up the total length of the contigs. 
        if line_split[-1] == "eukaryote":
            kingdom_dict["eukaryote"][0] += int(line_split[1])
            kingdom_dict["eukaryote"][1] += 1     
        elif line_split[-1] == "prokaryote":
            kingdom_dict["prokaryote"][0] += int(line_split[1])
            kingdom_dict["prokaryote"][1] += 1
        elif line_split[-1] == "archea":
            kingdom_dict["archea"][0] += int(line_split[1]) 
            kingdom_dict["archea"][1] += 1
            
    return kingdom_dict
         
main()