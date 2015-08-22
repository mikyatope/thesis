import sys              #needed for argparse
import argparse
import re               #regexp
from Bio import AlignIO

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
args = parser.parse_args()


# Match gaps and get region coordinates in Dmel
chrom = ""
start = 0
block = 0

alignment = AlignIO.read(args.inFile, "fasta")

for record in alignment :
    if re.match("D_mel", record.id):
    	gapsDmel    = re.finditer(r"(?<=[\w])-+(?=[\w])", str(record.seq)) # discard gaps at start or end of alignment
        gapsDmelAll = re.finditer(r"(-+)", str(record.seq))                # include gaps at start or end of alignment
        coords      = re.search('_(\w+):(\d+)-(\d+)_', record.id)
        chrom       = coords.group(1)
        start       = int(coords.group(2))
        block       = coords.group(0)

    if re.match("D_sim", record.id):
    	gapsDsim    = re.finditer(r"(?<=[\w])-+(?=[\w])", str(record.seq))
        gapsDsimAll = re.finditer(r"(-+)", str(record.seq))                # include gaps at start or end of alignment


    if re.match("D_yak", record.id):
    	gapsDyak    = re.finditer(r"(?<=[\w])-+(?=[\w])", str(record.seq))
        gapsDyakAll = re.finditer(r"(-+)", str(record.seq))                # include gaps at start or end of alignment


try:
    gapsDmel
except:
    sys.exit() 

try:
    gapsDsim
except:
    sys.exit() 

try:
    gapsDyak
except:
    sys.exit() 


# Create lists of gap spans
def getGapCoords(gaps, spans):
    for m in gaps:
        spans.append(m.span())

spansDmel=[]
spansDsim=[]
spansDyak=[]

spansDmelAll=[]
spansDsimAll=[]
spansDyakAll=[]

getGapCoords(gapsDmel, spansDmel)
getGapCoords(gapsDsim, spansDsim)
getGapCoords(gapsDyak, spansDyak)

getGapCoords(gapsDmelAll, spansDmelAll)
getGapCoords(gapsDsimAll, spansDsimAll)
getGapCoords(gapsDyakAll, spansDyakAll)



# allways add blank line (avoid problems with gnu-parallel mixing rows)
print "\n"

# check if ranges in A intersects in B.
def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


##### Polarized Dmel

# get exact gaps between Dyak Dsim
exactGapsDsimDyak = []
for gapDsim in spansDsim:
    for gapDyak in spansDyak:
        if gapDsim == gapDyak:
            exactGapsDsimDyak.append(gapDsim)

# check Deletions IN dmel (not present in Dmel)
for gapDmel in spansDmel:
    dontIntersect = True
    for gapDsim in spansDsimAll:
        if getOverlap(gapDmel, gapDsim) > 0:
           dontIntersect = False

    for gapDyak in spansDyakAll:
        if getOverlap(gapDmel, gapDyak) > 0:
           dontIntersect = False

    if dontIntersect:

        # get Absolute Dmel reference coordinates.
        # count gaps before indel
        gapsBefore = 0
        for m in spansDmelAll:
            if m < gapDmel:
                gapsBefore += m[1]-m[0]+1


        gapAbsStart = str((gapDmel[0]+start)-gapsBefore)
        gapAbsEnd   = str(int(gapAbsStart)+1)

        print chrom+"\t"+gapAbsStart+"\t"+gapAbsEnd+"\t"+"DEL\tpolDmel"+"\t"+block+"\t"+str(gapsBefore)+"\t"+args.inFile+"\t"+str(gapDmel[0])+"\t"+str(gapDmel[1]) # for deletions, as they don't 'exist' in Dmel, before and after nucleotide coordinates are provided

# check Insertions IN dmel (not present in Dyak or Dsim)
for gapDsimDyak in exactGapsDsimDyak:
    dontIntersect = True
    for gapDmel in spansDmelAll:
        if getOverlap(gapDsimDyak, gapDmel) > 0:
           dontIntersect = False

    if dontIntersect:
        gapsBefore = 0
        for m in spansDmelAll:
            if m < gapDmel:
                gapsBefore += m[1]-m[0]+1

        gapAbsStart = str((gapDsimDyak[0]+1+start)-gapsBefore)
        gapAbsEnd   = str((gapDsimDyak[1]+start)-gapsBefore)
        print chrom+"\t"+gapAbsStart+"\t"+gapAbsEnd+"\t"+"IN\tpolDmel"+"\t"+block+"\t"+str(gapsBefore)+"\t"+args.inFile+"\t"+str(gapDsimDyak[0])+"\t"+str(gapDsimDyak[1]) # start and end is 1based





##### Polarized Dyak
## coords in reference Dmel!!

# get exact gaps between Dmel Dsim
exactGapsDmelDsim = []
for gapDmel in spansDmel:
    for gapDsim in spansDsim:
        if gapDmel == gapDsim:
            exactGapsDmelDsim.append(gapDmel)

# check Deletions IN Dyak (not present in Dyak)
for gapDyak in spansDyak:
    dontIntersect = True

    for gapDmel in spansDmelAll:
        if getOverlap(gapDyak, gapDmel) > 0:
           dontIntersect = False

    for gapDsim in spansDsimAll:
        if getOverlap(gapDyak, gapDsim) > 0:
           dontIntersect = False

    if dontIntersect:

        # get Absolute Dmel reference coordinates.
        # count gaps before indel
        gapsBefore = 0
        for m in spansDmelAll:
            if m < gapDyak:                 # gaps in downstream Dmel seq but from Dyak gap!
                gapsBefore += m[1]-m[0]+1


        gapAbsStart = str((gapDyak[0]+start)-gapsBefore)
        gapAbsEnd   = str(int(gapAbsStart)+1)

        print chrom+"\t"+gapAbsStart+"\t"+gapAbsEnd+"\t"+"DEL\tpolDyak"+"\t"+block+"\t"+str(gapsBefore)+"\t"+args.inFile+"\t"+str(gapDyak[0])+"\t"+str(gapDyak[1]) # for deletions, as they don't 'exist' in Dmel, before and after nucleotide coordinates are provided


# check Insertions IN dyak (not present in Dmel or Dsim)
for gapDmelDsim in exactGapsDmelDsim:
    dontIntersect = True
    for gapDyak in spansDyakAll:
        if getOverlap(gapDmelDsim, gapDyak) > 0:
           dontIntersect = False

    if dontIntersect:
        gapsBefore = 0
        for m in spansDmelAll:
            if m < gapDmelDsim:                 # gaps in downstream Dmel seq but from DmelDsim gap!
                gapsBefore += m[1]-m[0]+1

        gapAbsStart = str((gapDmelDsim[0]+1+start)-gapsBefore)
        gapAbsEnd   = str((gapDmelDsim[1]+start)-gapsBefore)
        print chrom+"\t"+gapAbsStart+"\t"+gapAbsEnd+"\t"+"IN\tpolDyak"+"\t"+block+"\t"+str(gapsBefore)+"\t"+args.inFile+"\t"+str(gapDmelDsim[0])+"\t"+str(gapDmelDsim[1]) # start and end is 1based





##### Polarized Dsim
## coords in reference Dmel!!

# get exact gaps between Dmel Dyak
exactGapsDmelDyak = []
for gapDmel in spansDmel:
    for gapDyak in spansDyak:
        if gapDmel == gapDyak:
            exactGapsDmelDyak.append(gapDmel)

# check Deletions IN Dsim (not present in Dsim)
for gapDsim in spansDsim:
    dontIntersect = True

    for gapDmel in spansDmelAll:
        if getOverlap(gapDsim, gapDmel) > 0:
           dontIntersect = False

    for gapDyak in spansDyakAll:
        if getOverlap(gapDsim, gapDyak) > 0:
           dontIntersect = False

    if dontIntersect:

        # get Absolute Dmel reference coordinates.
        # count gaps before indel
        gapsBefore = 0
        for m in spansDmelAll:
            if m < gapDsim:                 # gaps in downstream Dmel seq but from Dyak gap!
                gapsBefore += m[1]-m[0]+1


        gapAbsStart = str((gapDsim[0]+start)-gapsBefore)
        gapAbsEnd   = str(int(gapAbsStart)+1)

        print chrom+"\t"+gapAbsStart+"\t"+gapAbsEnd+"\t"+"DEL\tpolDsim"+"\t"+block+"\t"+str(gapsBefore)+"\t"+args.inFile+"\t"+str(gapDsim[0])+"\t"+str(gapDsim[1]) # for deletions, as they don't 'exist' in Dmel, before and after nucleotide coordinates are provided


# check Insertions IN dsim (not present in Dyak or Dmel)
for gapDmelDyak in exactGapsDmelDyak:
    dontIntersect = True
    for gapDsim in spansDsimAll:
        if getOverlap(gapDmelDyak, gapDsim) > 0:
           dontIntersect = False

    if dontIntersect:
        gapsBefore = 0
        for m in spansDmelAll:
            if m < gapDmelDyak:                 # gaps in downstream Dmel seq but from DmelDyak gap!
                gapsBefore += m[1]-m[0]+1

        gapAbsStart = str((gapDmelDyak[0]+1+start)-gapsBefore)
        gapAbsEnd   = str((gapDmelDyak[1]+start)-gapsBefore)
        print chrom+"\t"+gapAbsStart+"\t"+gapAbsEnd+"\t"+"IN\tpolDsim"+"\t"+block+"\t"+str(gapsBefore)+"\t"+args.inFile+"\t"+str(gapDmelDyak[0])+"\t"+str(gapDmelDyak[1]) # start and end is 1based
