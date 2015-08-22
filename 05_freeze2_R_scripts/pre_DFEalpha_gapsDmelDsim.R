# Pre DFE-alpha computations


#### Annotate GAPS between Dmel-Dsim

GAPS_DmelDsim <- read.delim("~/SV_analysis/selection/DmelDsimGAPS.ALL.csv", header=F, quote="")
names(GAPS_DmelDsim) <- c("chr", "start")
#sort
GAPS_DmelDsim <- GAPS_DmelDsim[with(GAPS_DmelDsim, order(chr, start)), ]
row.names(GAPS_DmelDsim)<-NULL #reset row.names


## Annotate GAPS
## NOTE (jan 2014): Huge reused code, should create function/script some day
library("GenomicRanges")

##### --> previous GenomicRange objects for dmelAnnotations should exist

# create GenomicRange object for fixed divergence indels
gr.GAPS_DmelDsim <- with(GAPS_DmelDsim, GRanges(chr, IRanges(start-1, start))) ### GAPS as SNPs! 0 based!! 

foo <- findOverlaps(gr.GAPS_DmelDsim, dmelAnnotations.gr, type="within") 


# get Annotation names merging by subjectHit (subjectHit == annotation row.name)
foo_merged <- merge(as.data.frame(foo), dmelAnnotations, by.x="subjectHits", by.y="row.names", all.y=FALSE)[c(2,5:7)]

# get rownames for future merging (cross-reference the consensus annotation with its original row.name)
foo_merged$rn <- row.names(foo_merged)

### Use data.table to select consensus annotation for INDELs with multiple annotations
foo.DT <- as.data.table(foo_merged)

##### Convert annotations to numbers
# Get 'classes' vector  
vec <- as.character(foo.DT[,V3])
# recode vec to numbers
vec[vec!="three_prime_UTR" & vec!="five_prime_UTR" & vec!="CDS" & vec!="intron" & vec!="intergenic"]<-99
vec[vec=="CDS"]<-0
vec[vec=="five_prime_UTR" | vec=="three_prime_UTR"]<-1
vec[vec=="intron"]<-2
vec[vec=="intergenic"]<-3   
vec <- as.numeric(vec)
# Add recoded classes to table
foo.DT[,class:=vec]

# ORDER step (done automatically assigning the correct keys)
setkey(foo.DT, queryHits, class)  # orders by class!! (and the others)
foo.DT.ord    <- foo.DT[, list(classCode=class[1], 
                               rn=rn[1],                        
                               count=length(class), 
                               countCDS=length(class[class==0]),
                               countUTR=length(class[class==1]),
                               countIntron=length(class[class==2]),
                               countIntergenic=length(class[class==3]),
                               countOthers=length(class[class==99])), by=queryHits]


# merge with original annotations dataset to recover the annotation class name
GAPS_DmelDsimCl <- merge(as.data.frame(foo.DT.ord)[,c("rn", "queryHits", "classCode")],foo_merged)[,c("queryHits", "classCode","V3","V4","V5")]
names(GAPS_DmelDsimCl) <- c("queryHits", "classCode","class", "classStart","classEnd")
row.names(GAPS_DmelDsimCl)<-NULL #reset row.names

# get INDEL info names merging by queryHit (queryHit == annotation row.name)
GAPS_DmelDsimCl <- merge(GAPS_DmelDsim, GAPS_DmelDsimCl, by.x="row.names", by.y="queryHits", all.x=FALSE)

#sort
GAPS_DmelDsimCl <- GAPS_DmelDsimCl[with(GAPS_DmelDsimCl, order(chr, start)), ]
row.names(GAPS_DmelDsimCl)<-NULL #reset row.names
GAPS_DmelDsimCl <- subset(GAPS_DmelDsimCl, select=-c(Row.names))

rm(foo, foo_merged, foo.DT, foo.DT.ord)


##### add classes

#annotate Small/Long introns
GAPS_DmelDsimCl$classSize <- (GAPS_DmelDsimCl$classEnd-GAPS_DmelDsimCl$classStart)+1
GAPS_DmelDsimCl$class <- as.character(GAPS_DmelDsimCl$class) # make sure is not factor!

GAPS_DmelDsimCl$class2[GAPS_DmelDsimCl$class != "intron"] <- GAPS_DmelDsimCl$class[GAPS_DmelDsimCl$class != "intron"]
GAPS_DmelDsimCl$class2[GAPS_DmelDsimCl$class == "intron"] <- ifelse(GAPS_DmelDsimCl$classSize[GAPS_DmelDsimCl$class == "intron"] <= 100, "small_intron", "long_intron")

# Annotate X/Autosome in new column
GAPS_DmelDsimCl$aut_or_X <- ifelse(GAPS_DmelDsimCl$chr=="X", "X", "autosome")
row.names(GAPS_DmelDsimCl)<-NULL #reset row.names