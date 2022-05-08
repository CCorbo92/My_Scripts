# This script calculates the percent identity from aligned FASTA  
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: 05/2022
# Last Edit by: Christopher Corbo
import re
import sys

#Read in family name
family=sys.argv[1]
filestring=""

#Read in previously aligned fasta
with open(family+"_aligned.fasta","r") as readfile:
	for line in readfile:
		tmpfile = filestring + line.rstrip()
		filestring = tmpfile
IndFasta=[]
#THIS LINE ASSUMES the FASTA heading only has one set of parentheses and needs to be fixed. Probably count the length of the split and do the last index
#This section is editing the FASTA file to remove header data and separate each aligned FASTA seq into its own index in a list
new_list=filestring.split(")")
for i in range(1,len(new_list)-1):
	tmpList=new_list[i].split(">")
	IndFasta.append(tmpList[0])
IndFasta.append(new_list[len(new_list)-1])

# Compare each index in the pairwise comparison one at a time and increment on matches including matched gaps if they happen to exist.
for i in range(0,len(IndFasta)):
	for j in range(i+1,len(IndFasta)):
		similarity = 0
		for k in range(len(IndFasta[0])):
			if (IndFasta[i][k] == IndFasta[j][k]):
				similarity = similarity + 1
		percsim = similarity / len(IndFasta[0])
		print(str(i) + " " + str(j) + " " + str(percsim))
		

