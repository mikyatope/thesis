

###### PLOT INDEL Divergences

divergence <- read.delim("~/SV_analysis/divergence/ALL.fixed_indels.csv", header=F, quote="")
names(divergence) <- c("chr", "start", "stop", "in_or_del", "sp", "aln_block", "previous_gaps", "aln_file", "start_infile", "end_infile")
#sort
divergence <- divergence[with(divergence, order(chr, start)), ]
row.names(divergence)<-NULL #reset row.names

INDELs.freeze2 <- read.delim("../../share/dgrp.freeze2/ALL_SV_SNP/SV.noIV.uniq.onlyStart.csv", header=F, quote="")
names(INDELs.freeze2) <- c("chr", "start")

# Filter "false" divergence (INDELS polymorphic in melano)
gr.Divergenge <- with(divergence, GRanges(chr, IRanges(start, stop))) # from custom file, start ist stop-1 for deletions, insertions have its start and stop!
gr.INDELs     <- with(INDELs.freeze2, GRanges(chr, IRanges(start-1, start)))  # from vcf variants have only start coordinate

Div_INDELS_overlaps <- as.data.frame(findOverlaps(gr.Divergenge, gr.INDELs))

fdivergence <- divergence[!rownames(divergence) %in% Div_INDELS_overlaps$queryHits,]




################################
########## BY CLASSES ########## 
################################

## NOTE (jan 2014): Huge reused code, should create function/script some day


library("GenomicRanges")

##### --> previous GenomicRange objects for dmelAnnotations should exist

# create GenomicRange object for fixed divergence indels
gr.fDivergenge <- with(fdivergence, GRanges(chr, IRanges(start, stop)))

foo <- findOverlaps(gr.fDivergenge, dmelAnnotations.gr, type="within") 


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
fdivergenceCl <- merge(as.data.frame(foo.DT.ord)[,c("rn", "queryHits", "classCode")],foo_merged)[,c("queryHits", "classCode","V3","V4","V5")]
names(fdivergenceCl) <- c("queryHits", "classCode","class", "classStart","classEnd")
row.names(fdivergenceCl)<-NULL #reset row.names

# get INDEL info names merging by queryHit (queryHit == annotation row.name)
fdivergenceCl <- merge(fdivergence, fdivergenceCl, by.x="row.names", by.y="queryHits", all.x=FALSE)

#sort
fdivergenceCl <- fdivergenceCl[with(fdivergenceCl, order(chr, start)), ]
row.names(fdivergenceCl)<-NULL #reset row.names
fdivergenceCl <- subset(fdivergenceCl, select=-c(Row.names))

rm(foo, foo_merged, foo.DT, foo.DT.ord)


##### add classes

#annotate Small/Long introns
fdivergenceCl$classSize <- (fdivergenceCl$classEnd-fdivergenceCl$classStart)+1
fdivergenceCl$class <- as.character(fdivergenceCl$class) # make sure is not factor!

fdivergenceCl$class2[fdivergenceCl$class != "intron"] <- fdivergenceCl$class[fdivergenceCl$class != "intron"]
fdivergenceCl$class2[fdivergenceCl$class == "intron"] <- ifelse(fdivergenceCl$classSize[fdivergenceCl$class == "intron"] <= 100, "small_intron", "long_intron")

### Annotate frameshift indels
fdivergenceCl$class3 <- fdivergenceCl$class2
fdivergenceCl$size   <- (fdivergenceCl$end_infile-fdivergenceCl$start_infile) # coordinates zero based from python script
fdivergenceCl$class3[fdivergenceCl$size %% 3 == 0 & fdivergenceCl$class2 == "CDS"] <- "CDS_non_frameshift"
fdivergenceCl$class3[fdivergenceCl$size %% 3 != 0 & fdivergenceCl$class2 == "CDS"] <- "CDS_frameshift"

# Annotate X/Autosome in new column
fdivergenceCl$aut_or_X <- ifelse(fdivergenceCl$chr=="X", "X", "autosome")
row.names(fdivergenceCl)<-NULL #reset row.names








############
## GRAPHS ##
############

swDivergence <- function(slide, dataset) {
    
    Xchr   <- slide$chr
    Xstart <- slide$start
    Xend   <- slide$end
    
    wsize  <- Xend-Xstart+1
    slide$wsize <- wsize
    
    slideData<- subset(dataset, chr == Xchr & start >= Xstart & start <= Xend)

    slide$kdmel <- nrow(subset(slideData, sp == "polDmel"))/wsize
    slide$kdsim <- nrow(subset(slideData, sp == "polDsim"))/wsize
    slide$kdyak <- nrow(subset(slideData, sp == "polDyak"))/wsize
    
    return(slide)
}


#get slidding windows
SWs <- read.delim("~/SV_analysis/slidingWindows_100k.bed", header=F, quote="")
cols <- c("chr", "start", "end")
names(SWs) <- cols; rm(cols)
#testSlide <- SWs[c(1:5),]


## Coordinate slides loop
system.time( 
    DIV <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), swDivergence, dataset=fdivergence ))) 
)


ggplot(DIV, aes(x=start)) + 
    geom_point(size=2.2, aes(y=kdmel)) + 
    geom_point(size=2.2, aes(y=kdsim), colour="green") +
    geom_point(size=2.2, aes(y=kdyak), colour="blue") +
    facet_grid(. ~ chr, scales="free_x", space="free_x") +
    #coord_cartesian(ylim = c(0, 0.001)) +
    ylab("Divergence") +
    theme(legend.position = "none",axis.text.x=element_blank(), 
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(), 
          plot.margin=unit(rep(.5, 4), "lines")) #+
    #labs(title = " ") #


ggplot(cbind(PI, DIV$kdmel, correls_PI_adim$rec_comeron), aes(x=start, y=(PIindel_adim/DIV$kdmel)*correls_PI_adim$rec_comeron)) + 
    geom_point(size=2.2, aes(colour = "blue")) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +
    #coord_cartesian(ylim = c(0, 0.001)) +
    ylab("(Diversity/Divergence)*recombination") +
    theme(legend.position = "none",axis.text.x=element_blank(), 
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(), 
          plot.margin=unit(rep(.5, 4), "lines"))



ggplot(cbind(PI, DIV$kdmel, correls_PI_adim$rec_comeron), aes(x=start, y=PIindel_adim/DIV$kdmel)) + 
    geom_point(size=2.2, aes(colour = "blue")) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +
    #coord_cartesian(ylim = c(0, 0.001)) +
    ylab("Diversity/Divergence") +
    theme(legend.position = "none",axis.text.x=element_blank(), 
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(), 
          plot.margin=unit(rep(.5, 4), "lines"))


