import sys
file=sys.argv[1]
with open(file) as file_object:
    lines = file_object.readlines()
    for line in lines:
        stripline = line.rstrip()
        print(stripline)
