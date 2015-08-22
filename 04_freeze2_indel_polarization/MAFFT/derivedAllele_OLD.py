# derivedAllele.py
# 	by Miquel Ramia , miquel.ramia@gmail.cat
# 	Universitat Autonoma de Barcelona

# (no accents to ensure encoding)

# Infer derived/ancestral allele status for INDELS using one outgroup (right now, harcoded as the blast DB used!)

##### STEPS
# 1. Get indel coordinates from VCF, get a number of up and downstream flanking nucleotides.
# 2. Blast is done for each REF and ALT. Select longest hit sequence as outgroup sequence
#	(add blast criteria for valid hits!)
# 3. Perform multiple alignment for REF, ALT and Outgroup sequences (ALT and REF using the added flanking sequences)
#   (add alignment criteria to infer derived/ancestral state!)

### NOTES for future-Me or future-Another:
# open() already opens in batch mode, so no problems are expected with hughe files 
# (15min for 617454 indels and 50nuc flanking regions without parallelization)
# some stuff could be parallelized using this approach? --> http://stackoverflow.com/questions/2359253/solving-embarassingly-parallel-problems-using-python-multiprocessing

# Input:
#   - Refrence fasta genome (All chromosomes in single file)
# 	- VCF file -> at least first 5 columns. ID column (column[2], 3rd column) NOT used but maintains order
#   - Number of flanking nucleotides
#	- Outgroup identifier/name
	
# Output:
# 	- CSV file 




######## IMPORTS
import sys                                            # needed for argparse
import argparse                                   
import os                                             # needed for dir creation function
import errno                                          # needed for dir creation function
import subprocess                                     # needed for command line
from Bio import SeqIO                                 # to read fastas
from Bio import SearchIO                              # to parse BLAST output
from Bio import AlignIO                               # to parse alignments
from Bio.Align.Applications import MafftCommandline
from StringIO import StringIO                         # to create fake file objects from strings
import re                                             # regular expressions
from os.path import basename, splitext

######## Binary routes
blast_route     = "/home/mike/SV_analysis/blasts/ncbi-blast/bin/blastn"
blast_DBs_route = "/home/mike/SV_analysis/blasts/DBs/dsimV2"
mafft_route     = "/home/mike/SV_analysis/blasts/mafft/bin/mafft"



######## FUNCTIONS
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def parseBLAST(blastXML):

	blast_qresult = SearchIO.read(StringIO(blastXML), 'blast-xml')

	# initialize output variables
	num_hits  = 0
	num_hsps =evalue = qstart = qend = hit_seq = "NA"


	if (len(blast_qresult)):

		num_hits   = len(blast_qresult)

		# HIT related data 
		num_hsps         = len(blast_qresult[0])       # blast_qresult[0] -> num hsps of first hit
		
		# HSP related data ( blast_qresult[0][0] -> first HSP from first HIT)
		evalue    = blast_qresult[0][0].evalue
		qstart    = blast_qresult[0][0].query_start+1  # it's the XML coordinate -1
		qend      = blast_qresult[0][0].query_end
		hit_seq   = str(blast_qresult[0][0][0].hit.seq)     # this is a SeqRecord object converted to string (blast_qresult[0][0][0] -> third level, fragments of the HSP)

	return(num_hits, num_hsps, evalue, qstart, qend, hit_seq)


# Create lists of gap spans
def getGapCoords(gaps, spans):
	for m in gaps:
		spans.append(m.span())


# check if ranges in A intersects in B.
def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


######## MAIN PROGRAM
#####
## Get Arguments
#####
parser = argparse.ArgumentParser()
parser.add_argument('ref', help="Reference genome in fasta") #read stdin if file not provided
parser.add_argument('vcf', help="INDELs in VCF")             #read stdin if file not provided
parser.add_argument('flank', type=int, default=100, help="Number of flanking nucleotides per INDEL (default=100)")
parser.add_argument('outgroup', default="outgroup", help="Outgroup identifier")             
args = parser.parse_args()


### Read Fasta with ALL chromosomes. Index sequences with to_dict.
ref_handle = open(args.ref, 'r')
ref = SeqIO.to_dict(SeqIO.parse(ref_handle, "fasta"))
ref_handle.close()

### Read VCF (opened from argparse)
with open(args.vcf, 'r') as vcf:
	for line in vcf:

		# preinitialize MAIN OUTPUT, The derived allele state, and the sllices of aligned sequences and values
		derived  = "NA"
		sliceREF = "NA"
		sliceALT = "NA"
		sliceOUT = "NA"

		REF_num_hits=REF_num_hsps=REF_evalue=REF_qstart=REF_qend="NA"
		ALT_num_hits=ALT_num_hsps=ALT_evalue=ALT_qstart=ALT_qend="NA"


		col = line.split('\t') # split VCF columns


		### STEP 1: Flanking sequences for REF and ALT INDEL alleles
		# 'before' flanking sequence
		# seq is 0 based in Biopython!
		before_start = int(col[1])-args.flank-1   # before sequence start coord
		before_end   = int(col[1])-1              # before sequence end coord

		# 'after' flanking sequence
		after_start = int(col[1])+len(col[3])-1   # after sequence start coord
		after_end   = after_start+args.flank      # after sequence end coord

		# Get flanking sequences
		before_seq = ref[col[0]].seq[before_start:before_end]
		after_seq  = ref[col[0]].seq[after_start:after_end]

		# Create REF allele region fastqa
		REF_id    = "REF"
		REF_fasta = ">"+REF_id+"\n"+before_seq+col[3]+after_seq

		# Create ALT allele region fasta
		ALT_id    = "ALT"   
		ALT_fasta = ">"+ALT_id+"\n"+before_seq+col[4]+after_seq


		### STEP 2: BLASTs
		# Do BLASTs in command line and get output(blast output is in XML mode (5))    
		commandREF = "echo '" + str(REF_fasta) + "' | " + blast_route + " -db " + blast_DBs_route + " -outfmt 5"
		REF_BLAST_output = subprocess.check_output(commandREF, shell=True)

		commandALT = "echo '" + str(ALT_fasta) + "' | " + blast_route + " -db " + blast_DBs_route + " -outfmt 5"
		ALT_BLAST_output = subprocess.check_output(commandALT, shell=True)

		# Recover BLAST results
		REF_num_hits, REF_num_hsps, REF_evalue, REF_qstart, REF_qend, REF_hit_seq = parseBLAST(REF_BLAST_output)
		ALT_num_hits, ALT_num_hsps, ALT_evalue, ALT_qstart, ALT_qend, ALT_hit_seq = parseBLAST(ALT_BLAST_output)



		# NEXT INDEL if bad BLAST ('continue' works as 'next' in Perl)
		if ((REF_hit_seq == ALT_hit_seq == "NA") or (REF_num_hits > 1) or (ALT_num_hits > 1) or (REF_num_hits==ALT_num_hits==0)):
			print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t|{}|{}|{}\t{}\t{}".format(col[0], col[1], derived, REF_num_hits, REF_num_hsps, REF_evalue, REF_qstart, REF_qend, ALT_num_hits, ALT_num_hsps, ALT_evalue, ALT_qstart, ALT_qend, sliceREF, sliceALT, sliceOUT, col[3], col[4])
			continue 


		### STEP 3: Multiple Alignment
		# Create multi FASTA
		REF_hit_seq = REF_hit_seq.replace('-', '')
		ALT_hit_seq = ALT_hit_seq.replace('-', '')

		# Select longest outgroup sequence
		OUT_seq = REF_hit_seq
		if len(REF_hit_seq) < len(ALT_hit_seq):
			OUT_seq = ALT_hit_seq

		# Write multifasta file
		OUT_fasta = ">"+args.outgroup+"\n"+OUT_seq
		multiFasta = REF_fasta + "\n" + ALT_fasta + "\n" + OUT_fasta
		fileName     = "TEMP_"+splitext(basename(args.vcf))[0]+".mfa"
		f_multiFasta = open(fileName, 'w')
		f_multiFasta.write(str(multiFasta))
		f_multiFasta.close()

		# Execute MAFFT and get output
		commandMAFFT = mafft_route + " --globalpair --maxiterate 10 --inputorder --quiet "+fileName
		MAFFT_output = subprocess.check_output(commandMAFFT, shell=True)

		os.remove(fileName)

		alignment = AlignIO.read(StringIO(MAFFT_output), "fasta")

		# Get aligned sequences
		REF_aln = alignment[0].seq
		ALT_aln = alignment[1].seq
		OUT_aln = alignment[2].seq


		# GET INDEL region inside alignment 
		# start is variable 'before_end'(end of first flanking seq.) and end is 'after_start'(start of second flanking seq.)
		lenINDEL = abs(len(col[3])-len(col[4]))

		sliceREF = str(REF_aln[75:(125+lenINDEL)])
		sliceALT = str(ALT_aln[75:(125+lenINDEL)])
		sliceOUT = str(OUT_aln[75:(125+lenINDEL)])

		# complete gaps (whole gap inside slice)
		gapsREF  = re.finditer(r"(?<=[\w])-+(?=[\w])", sliceREF) # discard gaps at start or end of alignment region (onli "complete" gaps)
		gapsALT  = re.finditer(r"(?<=[\w])-+(?=[\w])", sliceALT) # discard gaps at start or end of alignment region (onli "complete" gaps)
		gapsOUT  = re.finditer(r"(?<=[\w])-+(?=[\w])", sliceOUT) # discard gaps at start or end of alignment region (onli "complete" gaps)
		spansREF=[]
		spansALT=[]
		spansOUT=[]
		getGapCoords(gapsREF, spansREF)
		getGapCoords(gapsALT, spansALT)
		getGapCoords(gapsOUT, spansOUT)

		# All gaps (also consider gaps that end outside the slice)
		allgapsREF  = re.finditer(r"-+", sliceREF) 
		allgapsALT  = re.finditer(r"-+", sliceALT) 
		allgapsOUT  = re.finditer(r"-+", sliceOUT) 
		allspansREF=[]
		allspansALT=[]
		allspansOUT=[]
		getGapCoords(allgapsREF, allspansREF)
		getGapCoords(allgapsALT, allspansALT)
		getGapCoords(allgapsOUT, allspansOUT)


		# get exact gaps between REF-OUT
		exactGapsREF_OUT = 0
		for gapREF in spansREF:
			for gapOUT in spansOUT:
				if gapREF == gapOUT:
					exactGapsREF_OUT += 1

		# get exact gaps between ALT-OUT
		exactGapsALT_OUT = 0
		for gapALT in spansALT:
			for gapOUT in spansOUT:
				if gapALT == gapOUT:
					exactGapsALT_OUT += 1



		# Cases with 2 gaps (REF or ALT with OUT). Derived is always sequence WITHOUT GAP
		if (len(spansREF) >= 1 or len(spansALT) >= 1):
			if (exactGapsREF_OUT > 0 or exactGapsALT_OUT > 0):
				if (exactGapsREF_OUT == 1 and exactGapsALT_OUT == 0):
					derived = "ALT"
				if (exactGapsALT_OUT == 1 and exactGapsREF_OUT == 0):
					derived = "REF"
			
			# Cases with 1 gap (no same gap in OUT)
			elif(exactGapsREF_OUT == 0 and exactGapsALT_OUT == 0):
				if len(spansREF) > 0:
					for currentGapREF in spansREF:                       # loop each REF gap INSIDE slide
						overlapsALT = overlapsOUT = 0
						
						for gapALT in allspansALT:
							if getOverlap(currentGapREF, gapALT) > 0:
								overlapsALT += 1

						for gapOUT in allspansOUT:
							if getOverlap(currentGapREF, gapOUT) > 0:
								overlapsOUT += 1

						if (overlapsALT+overlapsOUT==0):
							derived = "REF"

				if len(spansALT) > 0:
					for currentGapALT in spansALT:                       # loop each ALT gap INSIDE slide
						overlapsREF = overlapsOUT = 0
						
						for gapREF in allspansREF:
							if getOverlap(currentGapALT, gapREF) > 0:
								overlapsREF += 1
			
						for gapOUT in allspansOUT:
							if getOverlap(currentGapALT, gapOUT) > 0:
								overlapsOUT += 1

						if (overlapsREF+overlapsOUT==0):
							derived = "ALT"


		print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t|{}|{}|{}\t{}\t{}".format(col[0], col[1], derived, REF_num_hits, REF_num_hsps, REF_evalue, REF_qstart, REF_qend, ALT_num_hits, ALT_num_hsps, ALT_evalue, ALT_qstart, ALT_qend, sliceREF, sliceALT, sliceOUT, col[3], col[4])

vcf.close




##### Example write document
#make_sure_path_exists('fastas_REF')                            # create dir
#f_REF   = open("fastas_REF/"+REF_id+".fa", 'w')                # open file
#f_REF.writelines(">"+REF_id+"\n")                              # write id
#f_REF.writelines(before_seq+col[3]+after_seq)                  # write sequence
#f_REF.close()

