#This script assumes the full reference is the first molecule in the multimol2

import sys
# This is the number of molecules in the multimol2
molsize=int(sys.argv[2])

print("open " + sys.argv[1])
#Hide everything and then set up the reference which should be first in file
print("background solid white\n~disp\ncolor pink\ncolor byhet\ncolor dark gray #0.1\ncolor byhet #0.1\nsetattr m stickScale 0.2 #0.1\ndisp #0.1\ndelete HC")

#Dynamic Reference at each stage
print("color orange #0.3 ", end = ' ')

for i in range(8,molsize,4):
	print("#0." + str(i) + " ", end = '')

print("color byhet #0.3 ", end = ' ')

for i in range(8,molsize,4):
        print("#0." + str(i) + " ", end = '')

print("setattr m stickScale 0.5 #0.3 ", end = ' ')

for i in range(8,molsize,4):
        print("#0." + str(i) + " ", end = '')

#Make the heckin movie
print("movie record\nwait 55")

for i in range(2,molsize,4):
	print("disp #0." + str(i))
	print("~select #0." + str(i))
	print("wait 55")
	print("turn y 2.5 144")
	print("turn x 2.5 144")
	print("\n")

	print("disp #0." + str(i + 1))
	print("~select #0." + str(i + 1))
	print("wait 55")
	print("turn y 2.5 144")
	print("turn x 2.5 144")
	print("\n")

	print("~disp #0." + str(i))
	print("disp #0." + str(i + 2))
	print("~select #0." + str(i + 2))
	print("wait 75")
	print("turn y 2.5 144")
	print("turn x 2.5 144")
	print("~disp #0." + str(i + 1) + " #0." + str(i + 2))
	print("\n")

print("movie encode output ~/Downloads/" + sys.argv[3] + "_" + sys.argv[4] + "_movie.mp4")
	
