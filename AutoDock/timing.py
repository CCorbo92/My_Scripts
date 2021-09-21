import sys

file= sys.argv[1]
file1 = open(file,'r')
lines  =  file1.readlines()

file1.close()

for line in lines:
    line = line.rstrip()
    if "CPU=" in line:
        line = line.split(",")

        if "m" in line[0]:
           line = line[0].split()
           min=line[1]
           sec=line[2]
           min = (min[0])
           sec = (sec[0:2])
           total = (int(min) * 60) + int(sec)
           print(total)

        else:
           line = line[0].split()
           sec=line[1]
           sec = (sec[0:2])
           print(sec)
