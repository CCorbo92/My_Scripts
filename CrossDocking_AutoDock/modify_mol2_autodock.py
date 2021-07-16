# This script will modify a mol2 to match a format that satisfies mgltools

import sys

file = sys.argv[1]
		
file1 = open(file,'r')
lines  =  file1.readlines()
file1.close()



i = 0  
for line in lines:
    linesplit = line.split() #split on white space
    if (len(linesplit) == 1):
        if(linesplit[0] == "@<TRIPOS>ATOM"):
        	Var1 = i
    i = i + 1

i = 0  
for line in lines:
    linesplit = line.split() #split on white space
    if (len(linesplit) == 1):
        if(linesplit[0] == "@<TRIPOS>BOND"):
        	Var2 = i
    i = i + 1

for i in range (0,Var1 ,1):
	print(lines[i].rstrip())
print("@<TRIPOS>ATOM")

for i in range (Var1 + 1,Var2 ,1):
	splitline = lines[i].split()
	orig = splitline[7]
	append = int(splitline[6])
	done = (orig + str(append ))
	splitline[7] = done
	print(splitline[0] + " " + splitline[1] + " " + splitline[2] + " " + splitline[3] + " " + splitline[4] + " " + splitline[5] + " " + splitline[6] + " " + splitline[7] + " " + splitline[8] )


i = 0  
for line in lines:
    linesplit = line.split() #split on white space
    if (len(linesplit) == 1):
        if(linesplit[0] == "@<TRIPOS>SUBSTRUCTURE"):
        	Var3 = i
    i = i + 1

i = 0  
Var4 = len(lines)

for i in range (Var2 ,Var3 ,1):
	print(lines[i]),
print("@<TRIPOS>SUBSTRUCTURE")

for i in range (Var3 + 1,Var4 - 2,1):
	splitline = lines[i].split()
	orig = splitline[1]
	append = int(splitline[0]) 
	done = (orig + str(append))
	splitline[1] = done
	print(splitline[0] + " " + splitline[1] + " " + splitline[2] + " " + splitline[3] + " " + splitline[4] + " " + splitline[5] + " " + splitline[6] + " " + splitline[7])




