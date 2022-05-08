# This script aligns a multi FASTA file and prints out the alignment. It reads in the multi FASTA as a command line variable 
# Written by: Christopher Corbo
# Affiliation: Rizzo Lab, Stony Brook University
# Last Edit Date: 05/2022
# Last Edit by: Christopher Corbo
import sys

family=sys.argv[1]

from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
#Align sequences with MUSCLE (using parameters to make the alignment
#process as fast as possible)
muscle_cline = MuscleCommandline(input=(family+".fasta"), 
                                 out=(family+"_aligned.fasta"), 
                                 diags = True, 
                                 maxiters = 1, 
                                 log="Test_align_log.txt")
muscle_cline()

