#!/usr/bin/env python
import sys,os
import csv
import operator
import os.path
import glob
import time
from os import listdir
from os.path import isfile, join
from itertools import islice

start_time = time.time()
inputdirectory = sys.argv[1]
input_score = sys.argv[2]
amt_of_molecules = sys.argv[3]
scorelist_relevant_scores = ['Descriptor_Score:','Continuous_Score:','Pharmacophore_Score:','Hungarian_Matching_Similarity_Score:','Property_Volume_Score:','desc_FPS_vdw_fps:','desc_FPS_es_fps:','Footprint_Similarity_Score:','Total_Score:']

chunk_number = []
with open("chunk_list.txt",'r') as chunks:
    read_chunk = chunks.readlines()
    for line in read_chunk:
        chunk_number.append(line.rstrip())

temp_mol2 = []
with open(input_score+'_sorted_'+amt_of_molecules+'.csv','r') as file:
    for j,row_lines in enumerate(file,0):
        row_lines = row_lines.rstrip()
        for s in chunk_number:
            if s != row_lines.split(',')[48]:
                pass
            if s == row_lines.split(',')[48]:
                with open(inputdirectory+'/'+s+"_scored.mol2",'r') as f:
                    lines = islice(f, int(float(row_lines.split(',')[46]))-1, int(float(row_lines.split(',')[47])))
                    with open(input_score + "_top_"+amt_of_molecules+".mol2",'a') as initial_file:
                        for listitem in lines:  
                            initial_file.write(listitem)
                        initial_file.write("\n")




end_time = time.time()
print(end_time - start_time) 
print("end phase 3")
