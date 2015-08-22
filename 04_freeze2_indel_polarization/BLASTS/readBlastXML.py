
import sys              #needed for argparse
import argparse
import re               #regexp
from Bio import SearchIO
from os.path import basename, splitext

# readBlastXML.py
# 	by Miquel Ramia , miquel.ramia@gmail.cat
# 	Universitat Autonoma de Barcelona

# Read (one) blast XML. Output relevant data for 1st HIT and 1st HSP.

# Input:
#   - blast XML result
	
# Output:
# 	- STDOUT




#####
## MAIN
#####
parser = argparse.ArgumentParser()
parser.add_argument('xml', help="Reference genome in fasta") # read stdin if file not provided
args = parser.parse_args()

blast_qresult = SearchIO.read(args.xml, 'blast-xml')

XMLname    = splitext(basename(args.xml))[0]
nameParts  = XMLname.split("_")

# INDEL related data
chrom      = nameParts[0]
start      = nameParts[1]
len_indel  = nameParts[4][3:]


if (len(blast_qresult)):

	num_hits   = len(blast_qresult)

	# HIT related data
	num_hsps         = len(blast_qresult[0])     # num hsps of first hit

	# if (re.search('\[chromosome \= (.*)\]', blast_qresult[0].description).group(1)) is not None:
	# 	outgroup_chr   = re.search('\[chromosome \= (.*)\]', blast_qresult[0].description).group(1)
	# 	outgroup_start = blast_qresult[0][0].hit_start   # it's the XML coordinate -1 (it really is HSP data)

	# else:
	# 	outgroup_chr   = "NA"
	# 	outgroup_start = "NA"
	
	# HSP related data
	evalue    = blast_qresult[0][0].evalue
	qstart    = blast_qresult[0][0].query_start+1 # it's the XML coordinate -1
	qend      = blast_qresult[0][0].query_end
	gaps      = blast_qresult[0][0].gap_num
	ident     = blast_qresult[0][0].ident_num
	posit     = blast_qresult[0][0].pos_num
	total_aln = blast_qresult[0][0].aln_span

	print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrom, start, len_indel, num_hits, num_hsps, evalue, qstart, qend, gaps, ident, posit, total_aln)

else:
	print "{}\t{}\t{}\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA".format(chrom, start, len_indel)