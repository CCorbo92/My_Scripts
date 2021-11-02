# This script will best first cluster molecules.
# Needs a command line arg of file to read in

import sys


file=sys.argv[1]

name1 = []
name2 = []
namesused = []
#### Reading in file ####
with open(file) as file_object:
        lines = file_object.readlines()
        for line in lines:
           
                stripline = line.rstrip()
                splitline = stripline.split()
                name1.append(splitline[0])
                name2.append(splitline[1])

#### Done processing file and formatting for subsequent use ####

#### Initialzing Values ####
m = len(name1)
cluster_num=0
cluster_match = [0] * m
# Assessing similarity and clustering
clusters = [[] * m] * m 
for x in range(m):
   cluster= [] * m
   print(str(cluster_num))
   if not(name1[x]==name2[x] ):
      if (name2[x] not in namesused):
         cluster_num = cluster_num + 1
         clusterhead=name2[x]
         namesused.append(name2[x])
         cluster.append(name2[x])
         if (name1[x] not in namesused):
            namesused.append(name1[x])
            cluster.append(name1[x])

      if (name1[x] not in namesused and name2[x]==clusterhead):
         namesused.append(name1[x])
         cluster.append(name1[x])
   for n in range(len(cluster)):
      print(cluster[n])
#for x in range(len(clusters)):
   # these following 5 lines of code are designed to make unique removing the lower ranked duplicates
#   print( "Cluster #: " + str(x))
#   print(clusters[x]) 

