import sys              #needed for argparse
import argparse
import re               #regexp
from Bio import SeqIO

# callFixedINDELs.py
# 	by Miquel Ramia miquel.ramia@gmail.cat 
# 	Universitat Autonoma de Barcelona
#   20/08/2013

# Read (one) blast XML. Output relevant data for 1st HIT and 1st HSP.

# usage: 
# callFixedINDELs.py inFile

# Input:
#   - multi fasta file
	
# Output:
# 	- STDOUT


#####
## MAIN
#####

parser = argparse.ArgumentParser()
parser.add_argument('inFile', help="multiple alignment file in fasta") # read stdin if file not provided
parser.add_argument('label', help="Label to 1st column output") # read stdin if file not provided
args = parser.parse_args()


# Match gaps and get region coordinates in Dmel


with open(args.inFile, 'r') as content_file:
    
	content = content_file.read()
	positions  = re.finditer(r"4", content)

	for pos in positions:
		print args.label+"\t"+str(pos.start()+1)