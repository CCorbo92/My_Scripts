import csv
import operator
import sys
#file will be the csvto read in
file = sys.argv[1]
with open(file) as csvfile:
    reader = csv.reader(csvfile)
    #Descriptor Score
    sortedlist = sorted(reader, key=lambda i: float(i[7]))
    print("Descriptor_Score:")
    shortlist = (sortedlist[0:1000])
    for line in shortlist:
        print(', '.join(line))

with open(file) as csvfile:
    reader = csv.reader(csvfile)
    sortedlist = sorted(reader, key=lambda i: float(i[8]))
    print("Continuous_Score:")
    shortlist = (sortedlist[0:1000])
    for line in shortlist:
        print(', '.join(line))

with open(file) as csvfile:
    reader = csv.reader(csvfile)
    sortedlist = sorted(reader, key=lambda i: float(i[11]))
    print("Footprint_Similarity_Score:")
    shortlist = (sortedlist[0:1000])
    for line in shortlist:
        print(', '.join(line))


with open(file) as csvfile:
    reader = csv.reader(csvfile)
    sortedlist = sorted(reader, key=lambda i: float(i[16]))
    print("desc_FPS_es_fps:")
    shortlist = (sortedlist[0:1000])
    for line in shortlist:
        print(', '.join(line))

with open(file) as csvfile:
    reader = csv.reader(csvfile)
    sortedlist = sorted(reader, key=lambda i: float(i[17]))
    print("desc_FPS_vdw_fps:")
    shortlist = (sortedlist[0:1000])
    for line in shortlist:
        print(', '.join(line))

with open(file) as csvfile:
    reader = csv.reader(csvfile)
    sortedlist = sorted(reader, key=lambda i: float(i[22]))
    print("Pharmacophore_Score:")
    shortlist = (sortedlist[0:1000])
    for line in shortlist:
        print(', '.join(line))

with open(file) as csvfile:
    reader = csv.reader(csvfile)
    sortedlist = sorted(reader, key=lambda i: float(i[36]))
    print("Hungarian_Matching_Similarity_Score:")
    shortlist = (sortedlist[0:1000])
    for line in shortlist:
        print(', '.join(line))

with open(file) as csvfile:
    reader = csv.reader(csvfile)
    sortedlist = sorted(reader, key=lambda i: float(i[40]))
    print("Property_Volume_Score:")
    shortlist = (sortedlist[0:1000])
    for line in shortlist:
        print(', '.join(line)) 
