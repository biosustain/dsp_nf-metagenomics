def main():
    #list of the names of the four samples. 
    samples = ["O1","O2","O3","O4"]
    #Loops through all samples 
    for sample in samples:
        #Function that opens file and splits data into lines. 
        data_file = open_file(sample)
        #Function that loops through data and counts number of contigs and length for eukaryotes , prokaryotes and archeas.
        counts = analyse_data(data_file)
        print(counts)
        with open("/work3/apca/orange_peel/05_whokaryote/whokaryote/"+sample+"_whokaryote_stats.txt", 'w') as f:
            f.write('%s\t%s\t%s\n' % ('Kingdom', 'Total length (bp)', '# contigs'))  
            for key, value in counts.items():
                f.write('%s\t%s\t%s\n' % (key, value[0], value[1]))
        
def open_file(sample_name) -> list[str]:
    #Defines file location depending on sample name and opens file
    file_location = "/work3/apca/orange_peel/05_whokaryote/whokaryote/"+sample_name+"/featuretable_predictions_T.tsv"
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