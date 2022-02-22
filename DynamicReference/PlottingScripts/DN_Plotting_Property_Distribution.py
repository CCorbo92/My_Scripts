import sys,os
import csv
import matplotlib.mlab as mlab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import scipy.stats
from scipy.stats import norm
from scipy import stats

x_list =     []
y_list =     []
z_list =     []
xx_list =    []
yy_list =    []
zz_list =    []
xxx_list =   []
csv_total =  []

pdb = sys.argv[1]
desc_list = sys.argv[2]
N=sys.argv[3]

with open(desc_list,'r') as desc_list_read:
    descriptors = desc_list_read.readlines()
for descs in descriptors:

    x_list =     []
    y_list =     []
    z_list =     []
    xx_list =    []
    yy_list =    []
    zz_list =    []
    xxx_list =   []
    csv_total =  []

    desc=descs.rstrip()
    print(desc)
    file1="All_Results/"+pdb+"_"+desc+"_"+N+".txt"
    file2="All_Results_NonDynamic/"+pdb+"_"+desc+"_"+N+".txt"

    with open(file1,"r") as fileRead5:
        readFile5 = fileRead5.readlines()
    with open(file2,"r") as fileRead6:
        readFile6 = fileRead6.readlines()
    if (1==1):
        newFile5 = readFile5
        newFile6 = readFile6

        for line5 in newFile5:
            yy_list.append(float(line5))
        for line6 in newFile6:
            zz_list.append(float(line6))
    #    for line7 in readFile7:
    #        xxx_list.append(float(line7))
        fig, ax = plt.subplots(1, figsize=(8, 6))

        kdeyy = stats.gaussian_kde(yy_list,0.15)
        yy = np.linspace(min(yy_list), max(yy_list), 1000)
        kdezz = stats.gaussian_kde(zz_list,0.15)
        zz = np.linspace(min(zz_list), max(zz_list), 1000)

        ax.plot(zz, kdezz(zz) ,color='blue', label= "Standard")
        ax.plot(yy, kdeyy(yy) ,color='red', linestyle = ':', label= "Dynamic")

        plt.legend(loc=(0.62,0.850), title="Legend", frameon=True)
       
      
        plt.suptitle("Distribution of " + desc+" Scores: " + pdb ,fontsize=18)
        plt.xlabel("Top "+N+" Scoring Mol: Associated "+desc, fontsize=12) 
        plt.ylabel('Frequency (% of # of scores/total # of scores)') 
        #plt.show()
        plt.savefig(pdb+"_"+N+"_"+desc+".png") 
