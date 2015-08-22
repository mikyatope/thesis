
import sys              #needed for argparse
import argparse
import os               #needed for dir creation function
import errno            #needed for dir creation function
from Bio import SeqIO

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

# createFastas.py
# 	by Miquel Ramia , miquel.ramia@gmail.cat
# 	Universitat Autonoma de Barcelona

# Get indel coordinates from VCF, get a number of up and downstream flanking 
# nucleotides and output a fasta file for REF and ALT.

# open() already opens in batch mode, so no problems are expected with hughe files 
# (15min for 617454 indels and 50nuc flanking regions without parallelization)
# some stuff could be parallelizes with this approach? --> http://stackoverflow.com/questions/2359253/solving-embarassingly-parallel-problems-using-python-multiprocessing

# Input:
#   - Refrence fasta genome
# 	- VCF file (at least first 5 columns. ID not used but maintains order)
#   - Number of flanking nucleotides
	
# Output:
# 	- (automated fasta files for REF and ALT)




#####
## Get Arguments
#####
parser = argparse.ArgumentParser()
parser.add_argument('ref', help="Reference genome in fasta") #read stdin if file not provided
parser.add_argument('vcf', help="INDELs in VCF") #read stdin if file not provided
parser.add_argument('flank', type=int, default=100, help="Number of flanking nucleotides per INDEL (default=100)")
args = parser.parse_args()


### Read Fasta. Index sequences with to_dict.
ref_handle = open(args.ref, 'r')
ref = SeqIO.to_dict(SeqIO.parse(ref_handle, "fasta"))
ref_handle.close()

### Read VCF (opened from argparse)
with open(args.vcf, 'r') as vcf:
	for line in vcf:
	    col = line.split('\t') # split VCF columns
	    

	    # Downstream flanking sequence
	    # seq is 0 based in Biopython!
	    dfl_start = int(col[1])-args.flank-1   # downstream Flank start
	    dfl_end   = int(col[1])-1              # downstream Flank end

	    # Upstream flanking sequence
	    ufl_start = int(col[1])+len(col[3])-1   # upstream Flank start
	    ufl_end   = ufl_start+args.flank    # upstream Flank end

	    # Get flanking sequences
	    dfl_seq = ref[col[0]].seq[dfl_start:dfl_end]
	    ufl_seq = ref[col[0]].seq[ufl_start:ufl_end]

	    # Create REF allele region fasta
	    make_sure_path_exists('fastas_REF')                            # create dir
	    REF_id  = col[0]+"."+str(int(col[1]))+".REF."+str(args.flank)  # id string
	    REF_id += "_len"+str(len(col[3]))
	    f_REF   = open("fastas_REF/"+REF_id+".fa", 'w')                # open file
	    f_REF.writelines(">"+REF_id+"\n")                              # write id
	    f_REF.writelines(dfl_seq+col[3]+ufl_seq)                       # write sequence
	    f_REF.close()


	    # Create ALT allele region fasta
	    make_sure_path_exists('fastas_ALT')                            # create dir
	    ALT_id  = col[0]+"."+str(int(col[1]))+".ALT."+str(args.flank)  # id string
	    ALT_id += "_len"+str(len(col[4]))
	    f_ALT   = open("fastas_ALT/"+ALT_id+".fa", 'w')                # open file
	    f_ALT.writelines(">"+ALT_id+"\n")                              # write id
	    f_ALT.writelines(dfl_seq+col[4]+ufl_seq)                       # write sequence
	    f_ALT.close()

	    # print ">"+col[0]+"_"+str(int(col[1]))+"_ALT("+col[4][:10]+")"
	    # print dfl_seq+"-"+col[4]+"-"+ufl_seq

	    # print ref[col[0]].seq[dfl_start:ufl_end]
	    # print ""

vcf.close

# ### Get Reference
# for ref in SeqIO.parse(args.ref, "fasta"):
#     print ref.id
#     print ref.seq[0:3]
#     print len(ref)


