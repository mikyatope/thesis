import sys              #needed for argparse
import argparse
import re               #regexp
from Bio import SeqIO

# recoded2tabulated.py.py
# 	by Miquel Ramia miquel.ramia@gmail.cat 
# 	Universitat Autonoma de Barcelona
#   29/01/2014

# Read recoded sequence, output the same but tabulated format (1 position per line)

# usage: 
# recoded2tabulated.py <inFile> <chr> 

# Input:
#   - recoded sequence
	
# Output:
# 	- STDOUT (tabulated format)


#####
## MAIN
#####

parser = argparse.ArgumentParser()
parser.add_argument('inFile', help="recoded sequence") # read stdin if file not provided
parser.add_argument('label', help="Label to 1st column output (chr)") # read stdin if file not provided
args = parser.parse_args()


# Match gaps and get region coordinates in Dmel


with open(args.inFile, 'r') as content_file:
    
	recoded = content_file.read()
	position = 1

	for pos in recoded:
		print args.label+"\t"+str(position)+"\t"+pos
		position = position+1