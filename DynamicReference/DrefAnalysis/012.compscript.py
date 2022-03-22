
with open("ref_comp.txt","r") as fileRead:
	readFile = fileRead.readlines()
    
	for line in readFile:
		splitline=line.split()
		if (float(splitline[0]) > float(splitline[1])):
			print(splitline[0] + " " + splitline[1] )
