# This script will best first cluster molecules.
# Needs a command line arg of file to read in

import time
from joblib import Parallel, delayed
import multiprocessing
import sys

startTime = time.time()

#num_cores can change at last round but need a copy of its original value
num_cores = multiprocessing.cpu_count()
step_size = num_cores

file=sys.argv[1]

names = []
fingerprints = []

#### Reading in file ####
with open(file) as file_object:
        lines = file_object.readlines()
        for line in lines:
                maccs_split = []
                stripline = line.rstrip()
                splitline = stripline.split(",")
                names.append(splitline[0])
                maccs_raw = splitline[1]

                # Script requires that MACCS must be split into a list of bits
                for bit in maccs_raw:
                        maccs_split.append(int(bit))

                fingerprints.append(maccs_split)
#### Done processing file and formatting for subsequent use ####

#### Initialzing Values ####

m = len(names)
#Length of MACCS fp
n = 167
# This is a numerical list which will be trimmed parallel to the names and fingerprints lists, so a record of their original index is maintained
index_const = []
for i in range(m):
	index_const.append(i)

cluster_match = [0] * m
removed = 0
clusters_created = 0
# Assessing similarity and clustering
clusters = [[] * m] * m
# clustering threshold
threshold = 0.95
cluster_num = 1
# required for staggered for loop
round = 0
num_cycles = int(m / num_cores)
all_tanimoto = []

#### Done initializing values ####

#### Tanimoto Calculation is submitted in parallel and takes a single molecule and calculates its Tanimoto to the rest of the remaining list
def tanimoto_calc(k):
	# Looping through one cluster creation at a time
	# m is a dynamic variable which changes as the list gets shorter and needs to change to avoid out of range for i loop
	m = len(names)
	# Temporary place holder for similarity between two strings (Size of Positive Intersection)  
	similarity = [0] * m
	# Temporary place holder for features of two strings (Size of Positive Union)
	num_features = [0] * m
	# Temporary place holder for Tanimoto between two strings (Positive Intersection / Positive Union)
	tanimoto = [0] * m

	# Looping through remaining strings in list
	for i in range(m):
		# Looping through each bit in two comparison strings
		for j in range(n):
			if(fingerprints[k][j] == 1 or fingerprints[i][j] ==1):
				num_features[i] = num_features[i] + 1
				if (fingerprints[k][j] == fingerprints[i][j]):
					similarity[i] = similarity[i] + 1
	
		tanimoto[i] = similarity[i] / num_features[i]

	return tanimoto;
	
#### End of Tanimoto Calculation

#### Beginning Clustering ####

for x in range(num_cycles):
	finish = (round + 1) * num_cores 
	finish = finish - removed
	round = round + 1
	if (x == num_cycles - 1):
		finish = len(names)

	# The following 2 if statements are to avoid out of range problems 
	if (len(names) == 0):
		break
		print("BREAK")
	if (len(names) < num_cores):
		num_cores = len(names)

	# Where the tanimoto calculation takes places in parallel
	out = Parallel(n_jobs=num_cores)(delayed(tanimoto_calc)(k) for k in range(0,finish))
	
	for y in range(len(out)):
		all_tanimoto.append(out[y])

	# do some clustering on this round of all_tanimoto
	for k in range(len(all_tanimoto)):
		m = len(names)
		cluster = [] * m
		index_remove = [] * m
		for j in range(m):
			if (all_tanimoto[k][j] > threshold):
				cluster.append(names[j])
				index_remove.append(j)
				removed = removed + 1
				cluster_match[index_const[j]] = k + (step_size * x)

		#This will create a list of clusters which are lists of molecules
		clusters[cluster_num]=cluster
		cluster_num = 1 + cluster_num

		# Need to remove in reverse indexing order to not affect actual order 
		index_remove.reverse()

		# Remove items which have been clustered
		for index in index_remove:
			del fingerprints[index]
			del names[index]
			del index_const[index]
			for i in range(len(all_tanimoto)):
				del all_tanimoto[i][index]

	#This resets all tanimoto after each round	
	all_tanimoto = []

#### Done clustering ####

for x in range(len(clusters)):
   print( "Cluster #: " + str(x))
   print(clusters[x])
   print(len(clusters[x]))

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))
