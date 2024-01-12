#!/usr/bin/env python3

# Importing the whokaryote directory and the output file name
import sys
input_args = sys.argv[1:]

def main():
    
	# Need to take whokaryote dir fro a params variable and also output file name!
    from os import listdir
    samples=listdir(input_args[0])
    sorted_samples = sorted(samples)
    print(sorted_samples)
	
	#Loops through all samples 
    with open(input_args[1], 'w') as f:
        f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('Sample', 'Total_length_bp_Euk', 'num_contigs_Euk', 'Total_length_bp_Prok', 'num_contigs_Prok', 'Total_length_bp_Arc', 'num_contigs_Arc'))
        for sample in sorted_samples:
            #Function that opens file and splits data into lines. 
            data_file = open_file(sample)
            #Function that loops through data and counts number of contigs and length for eukaryotes , prokaryotes and archeas.
            counts = analyse_data(data_file)
            print(counts)
            # we add sample name to the contigs data
            in_one_line = [sample]
			# We add all the values to one line
            for key, value in counts.items():
                in_one_line += value
			# Preparing the data as string to write in a file
            as_strings = [str(x) for x in in_one_line]
            final_string = '\t'.join(as_strings)
            final_string += '\n'
            f.write(final_string)
        
def open_file(sample_name) -> list[str]:
    #Defines file location depending on sample name and opens file
    file_location = "whokaryote/"+sample_name+"/featuretable_predictions_T.tsv"
    file = open(file_location)
    #splits the data into lines
    data = file.read().split("\n") 
    #Prints sample name
    print(sample_name)
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