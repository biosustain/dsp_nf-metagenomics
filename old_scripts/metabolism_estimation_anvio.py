import numpy as np

def main():
#should change line to whatever is better name i.e read or so
#should also change "reads" to whatever it actually is

	#Open file
	file = open("04_anvio/kegg-metabolism_modules_bins_filtered_075.txt")

	#splits the data into lines
	data = file.read().split("\n")
	Leuconostoc_citreum_bin1 = []
	Eukaryotes_bin2 = []
	Lactiplantibacillus_bin3 = []
	Levilactobacillus_brevis_bin4 = []
	Unknown_bin5 = []
	Unknown_bin6 = []

	for line in data:
		line_split = line.split("\t")
		
		if "Bin_1" in line_split[1:2]:
			Leuconostoc_citreum_bin1.append(line_split[2:3])
		if "Bin_2" in line_split[1:2]:
			Eukaryotes_bin2.append(line_split[2:3])
		if "Bin_3" in line_split[1:2]:
			Lactiplantibacillus_bin3.append(line_split[2:3])
		if "Bin_4" in line_split[1:2]:
			Levilactobacillus_brevis_bin4.append(line_split[2:3])
		if "Bin_5" in line_split[1:2]:
			Unknown_bin5.append(line_split[2:3])
		if "Bin_6" in line_split[1:2]:
			Unknown_bin6.append(line_split[2:3])

	print("Leuconostoc_citreum_bin1:", Leuconostoc_citreum_bin1,"\n\n")
	print("Eukaryotes_bin2:", Eukaryotes_bin2,"\n\n")
	print("Lactiplantibacillus_bin3:", Lactiplantibacillus_bin3,"\n\n")
	print("Levilactobacillus_brevis_bin4:", Levilactobacillus_brevis_bin4,"\n\n")
	print("Unknown_bin5:", Unknown_bin5,"\n\n")
	print("Unknown_bin6:", Unknown_bin6,"\n\n")

	arrays = [Leuconostoc_citreum_bin1, Eukaryotes_bin2, Lactiplantibacillus_bin3, Levilactobacillus_brevis_bin4, Unknown_bin5, Unknown_bin6]

	matrix = [[0 for _ in range(6)] for _ in range(6)]

	for i in range(6):
		for j in range(6):
			if i == j:
				matrix[i][j] = len(arrays[i])
			else:
				count = 0
				for string in arrays[i]:
					if string in arrays[j]:
						count += 1
				matrix[i][j] = count
	
	# reversing column order in matrix
	arr = np.array([['Bin 1', 'Bin 2', 'Bin 3', 'Bin 4', 'Bin 5', 'Bin 6'],
		matrix[5],
		matrix[4], 
		matrix[3], 
		matrix[2],
		matrix[1],
		matrix[0]])
	
	# reversing column order in matrix	
	flipped_arr = np.fliplr(arr)

	# Printing on the screen
	print('\nArray after changing column order as in  Tiffanys tables:\n', flipped_arr)
	
	with open("/work3/apca/orange_peel/04_anvio/Kegg_pathways_on_bins.txt", 'w') as f:
		f.write("Array after changing column order as in  Tiffanys tables:\n"+str(flipped_arr)+"\n")
	f.close()

main()
