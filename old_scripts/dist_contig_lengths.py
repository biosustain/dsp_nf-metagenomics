import pandas as pd
import numpy

def main():
    file = open("/work3/apca/orange_peel/02_assembly/coassembly/final.contigs.fa")
    #splits the data into lines
    data = file.read().split("\n")
    #Prints sample name
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
                    
    with open("/work3/apca/orange_peel/02_assembly/distribution_contigs_coassembly.txt", 'w') as f:  
        f.write("over 20000, length "+str(len(lengths_20000))+"\n")
        f.write("max "+str(max(lengths_20000)/20000)+"\n")
        f.write("total contig length over 20000: "+str(sumk)+"\n")
        f.write("Total contig length: "+str(sumi)+"\n")
        f.write("Total contig length - total contig length over 20000: "+str(sumi-sumk)+"\n")
    f.close()
                
main ()