#!/usr/bin/env python
import sys,os
import csv
import operator
import os.path
import glob
import time
import codecs
from os import listdir
from os.path import isfile, join

start_time = time.time()
amt_of_molecules = sys.argv[1]
scorelist_relevant_scores = ['Descriptor_Score:','Continuous_Score:','Pharmacophore_Score:','Hungarian_Matching_Similarity_Score:','Property_Volume_Score:','desc_FPS_vdw_fps:','desc_FPS_es_fps:','Footprint_Similarity_Score:','Total_Score:' ]


score_names_list = []
with open("score_name_list.txt", 'r') as score_names_file:
    read_score_names = score_names_file.readlines()
    for line in read_score_names:
        score_names_list.append(line.rstrip()) 

print(score_names_list)
with open('all_molecules.csv','r') as file_reader:
    csv_read_file = csv.reader(file_reader, delimiter = ',')

    csv_temp = []
    for csv_lines in csv_read_file:
        csv_temp.append(csv_lines)
    for i,element in enumerate(score_names_list,0):
       if element not in str(scorelist_relevant_scores):
           pass
       if element in str(scorelist_relevant_scores): 
           with open(element+'_sorted.csv', 'w', newline='') as file:
               sorted_list = []
               writer = csv.writer(file)
               writer.writerow(score_names_list)
               sorted_list = sorted(csv_temp,key=lambda x: float(x[i]), reverse = False)
               for line in sorted_list:
                   writer.writerow(line)


for element in scorelist_relevant_scores:
    with open(element+'_sorted.csv', 'r') as file_reader:
        csv_read_file = csv.reader(file_reader, delimiter = ',')
        with open(element+'_sorted_'+amt_of_molecules+'.csv', 'w', newline='') as file_writer:
           writer = csv.writer(file_writer)
           for i,line in enumerate(csv_read_file,0):
               if i > 0 and i <= int(amt_of_molecules):
                   writer.writerow(line)

end_time = time.time() 
print("end phase 2")
