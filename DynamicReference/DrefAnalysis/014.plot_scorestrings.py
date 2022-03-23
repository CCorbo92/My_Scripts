import matplotlib.pyplot as plt
import numpy as np

x=[1,2,3,4,5,6,7,8,9,10,11]

with open("scorestrs.csv","r") as fileReader:
    readFile = fileReader.readlines()
    
    string_dyn= []
    stringsplit_flt=[]
    for line in readFile:
        string=line.rstrip()
        stringsplit = string.split(",")

        length = len(stringsplit)
        add = 11 - length
        #print(stringsplit)
        for num in stringsplit:
            stringsplit_flt.append(float(num))
        for i in range(add):
            stringsplit_flt.append(0)
        ypoints = np.array(stringsplit_flt)
        ypoints[ypoints==0] = np.NaN
        print(x)
        print(ypoints)
        plt.plot(x,ypoints, linestyle = 'solid',color="red")

        stringsplit_flt=[]

with open("scorestrs_dyn.csv","r") as fileReader:
    readFile = fileReader.readlines()
    
    string_dyn= []
    stringsplit_flt=[]
    for line in readFile:
        string=line.rstrip()
        stringsplit = string.split(",")

        length = len(stringsplit)
        add = 11 - length
        #print(stringsplit)
        for num in stringsplit:
            stringsplit_flt.append(float(num))
        for i in range(add):
            stringsplit_flt.append(0)
        ypoints = np.array(stringsplit_flt)
        ypoints[ypoints==0] = np.NaN
        print(x)
        print(ypoints)
        plt.plot(x,ypoints, linestyle = 'solid',color="green")

        stringsplit_flt=[]

plt.xlabel("Added Segments")
plt.ylabel("Score: 25*HMS + 1*GRD")
plt.title("Score Strings For 2GQG on Full Reference")
plt.show()
