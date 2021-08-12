import sys,os
import csv
import operator
import os.path
import glob
import time
from os import listdir
from os.path import isfile, join

start_time = time.time()
inputdirectory = sys.argv[1]
total_paths = [filenames_for_pre for filenames_for_pre in listdir(inputdirectory) if isfile(join(inputdirectory, filenames_for_pre))]

print("gathered file paths and names")
chunk_number = []
for score_filenames in total_paths:
    score_filenames = score_filenames.replace("_scored.mol2","")
    chunk_number.append(score_filenames)
print(chunk_number)

with open("chunk_list.txt", 'w') as chunk_file:
    for line in chunk_number:
        chunk_file.write(line)
        chunk_file.write("\n")
         

dict = {}
for i in chunk_number:
    dict[i] = []
print(dict)

def open_chunk (chunk_number):
    with open(inputdirectory + '/'+chunk_number+"_scored.mol2", 'r') as chunk_file:
        read_files = chunk_file.readlines()
        return read_files

dict_score_names = {}
for i in chunk_number:
    dict_score_names[i] = []
print(dict_score_names)

print(" adding new columns and extracting the descriptor scores")
for i in chunk_number:
    scorecopy = False
    for line in open_chunk(i):
        print(line)
        if len(line.split()) == 3 and "##########" in line.split()[0]:
            scorecopy = True
        if len(line.split()) == 3 and scorecopy == True:
            dict_score_names[i].append(line.split()[1])
        if "ROOT" in line:
            break
        scorecopy = False
    dict_score_names[i].append("Total_Score:")
    dict_score_names[i].append("Line_Start")
    dict_score_names[i].append("Line_End")
    dict_score_names[i].append("chunk_source")

#adding header
score_names_list = []
for score_names in dict_score_names[i]:
    score_names_list.append(score_names)

with open("score_name_list.txt", 'w') as score_names_file:
    for line in score_names_list:
        score_names_file.write(line)
        score_names_file.write("\n")

print("extracted just the descriptor scores for sorting...")

print("makeing all_molecule.csv")
total_list = []
temp_contin = 0
temp_fps_sum = 0
with open('all_molecules.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    for s in chunk_number:
        row = []
        for i,line in enumerate(open_chunk(s),1): 
            if "##########" not in line:
                pass
            if "##########" in line and str(line.split()[2]):
                row.append(line.split()[2])
            elif "##########" in line and int(line.split()[2]):
                row.append(float(line.split()[2]))
            if "##########" in line and "Continuous_Score:" in line:
                temp_contin = float(line.split()[2])
            if "##########" in line and "Footprint_Similarity_Score:" in line:
                temp_fps_sum = float(line.split()[2])
            if len(row) == (len(dict_score_names[s])-4):
                row.append(temp_contin + temp_fps_sum)
                temp_contin = 0
                temp_fps_sum = 0
            if len(row) == (len(dict_score_names[s])-3):
                row.append(float(i-(len(dict_score_names[s])-5)))
            if "ROOT" in line:
                row.append(float(i))
                row.append(s)
                writer.writerow(row)
                total_list.append(row)
                row = []



end_time = time.time()
print("end phase 1")
  
