import sys

file = sys.argv[1]

file_open = open(file,'r')
lines_raw=file_open.readlines()

file_open.close()
sum = 0
for i in range(len(lines_raw)):
	sum = 0
	line_split = lines_raw[i].split()
	for j in range(len(line_split)):
		sum = float(line_split[j]) + sum
number_round = round(sum,1)
string_num = str(number_round)
for i in range(len(string_num)):
	if (string_num[i] == "."):
		index = i

if (string_num[index + 1] != "0"):
	print("  Not Whole")

