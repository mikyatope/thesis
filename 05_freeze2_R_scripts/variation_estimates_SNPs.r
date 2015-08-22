# LOAD LIBS
#custom libs
source("~/SV_analysis/R_analisis/lib_mike_SV.r")

#graphics
library(ggplot2)
library(reshape2)
library(gridExtra)

#ddply
library(plyr)    


#########################
########## PI SNPS
########################

### GET DATA
#SNPdata <- read.delim("~/SV_analysis/freeze2.SNP.UNIQUE.pi.het.vcf", header=F, quote="")
system.time( SNPdata <- getAndFilterFreeze2SNP("~/SV_analysis/freeze2.SNP.UNIQUE.pi.het.vcf") )

#### total  
fSNPdata <- subset(SNPdata, MISS <= 20 & RC > 0 & FAC > 0)  # get indels with less than 10% missing lines and at least 1 line with the reference allele in homozigosis
#### total 4276393


### polarized SNPs
SNPs_polarized <- read.delim("~/SV_analysis/SNPs_polarized.csv", header=F)
names(SNPs_polarized) <- c("chr", "start", "REF", "ALT", "INFO", "dsim","class", "derived")


### MERGE with polarized
fSNPdata <- merge(fSNPdata, subset(SNPs_polarized, select=-c(REF,ALT,INFO,dsim)))
#sort
fSNPdata <- fSNPdata[with(fSNPdata, order(chr, start)), ]
row.names(fSNPdata)<-NULL #reset row.names



### GET CONSENSUS ANNOTATIONS.
### Give weight to each annotation CDS[0] > UTR[1] > Intron[2] > Intergenic > [3]
SNPs_class_polarized <- read.delim("~/SV_analysis/SNPs_crossClasses.csv", header=F)
names(SNPs_class_polarized) <- c("chr", "start",  "derived", "class", "classSize")

# get rownames for future merging
SNPs_class_polarized$rn <- row.names(SNPs_class_polarized)

# Create data.table for efficiency
library(data.table)
SNPs_class_polarized.dt <- data.table(SNPs_class_polarized[c("chr", "start", "class")], keep.rownames=T) # use data.table for efficiency inside loop || drop classInfo column

# Get 'classes' vector  
vec <-  as.character(SNPs_class_polarized.dt[,class])
# recode vec to numbers
vec[vec=="CDS"]<-0
vec[vec=="five_prime_UTR" | vec=="three_prime_UTR"]<-1
vec[vec=="intron"]<-2
vec[vec=="intergenic"]<-3   
vec <- as.numeric(vec)
# Add recoded classes to table, order by it.
SNPs_class_polarized.dt[,class:=vec]

### Get consensus annotation
setkey(SNPs_class_polarized.dt, chr, start, class)  # orders by class!! (and the others)

# ordered by the keys, class[1] is the class lowest value
SNPs_class_polarized.consensus.dt <- SNPs_class_polarized.dt[, list(classCode=class[1], 
                                                                    rn=rn[1], 
                                                                    count=length(class), 
                                                                    countCDS=length(class[class==0]),
                                                                    countUTR=length(class[class==1]),
                                                                    countIntron=length(class[class==2]),
                                                                    countIntergenic=length(class[class==3])), by=list(chr, start)]

# convert to data.frame
SNPs_class_polarized.consensus <- as.data.frame(SNPs_class_polarized.consensus.dt)
SNPs_class_polarized.consensus <- subset(SNPs_class_polarized.consensus, select=c("chr", "start", "rn"))

# merge with original data.frame, discard empty!
SNPs_class_polarized <- merge(SNPs_class_polarized, SNPs_class_polarized.consensus)

#remove rm column
SNPs_class_polarized <- subset(SNPs_class_polarized, select=-rn)

#sort
SNPs_class_polarized <- SNPs_class_polarized[with(SNPs_class_polarized, order(chr, start)), ]
row.names(SNPs_class_polarized)<-NULL #reset row.names

# remove variables
#rm(SNPs_class_polarized.consensus, SNPs_class_polarized.consensus.dt, vec)



###### MERGE WITH WORKING DATASET!
fSNPdata <- merge(fSNPdata, SNPs_class_polarized, by=c("chr", "start"))
fSNPdata <- subset(fSNPdata, select=-c(derived.y))
colnames(fSNPdata)[16] <- "class"   #from class.x
colnames(fSNPdata)[17] <- "derived" #from derived.x
colnames(fSNPdata)[18] <- "class2"  #from class.y

#sort
fSNPdata <- fSNPdata[with(fSNPdata, order(chr, start)), ]
row.names(fSNPdata)<-NULL #reset row.names




### ******Filter JGIL SNPs******** (this step had to be done before :( )

# # old.fSNPdata <- fSNPdata

jgilSNPs <- read.delim("../../share/dgrp.freeze2/freeze2.jgil.only_SNPs.csv", header=T, quote="", sep="") # separator is "white space"
names(jgilSNPs) <- c("chr", "start", "id", "refc", "altc", "qual", "cov")

fSNPdata <- merge(fSNPdata, jgilSNPs, all.y=FALSE)
fSNPdata <- fSNPdata[with(fSNPdata, order(chr, start)), ]
row.names(fSNPdata)<-NULL #reset row.names

jgilPI_SNPs <- 2*(fSNPdata$altc*fSNPdata$refc)/(((fSNPdata$altc+fSNPdata$refc)*(fSNPdata$altc+fSNPdata$refc-1)))
fSNPdata$jgilPI <- jgilPI_SNPs
rm(jgilPI_SNPs)



system.time( 
    PI_SNP  <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), sumAndMeansSNP, dataset=fSNPdata ) )) 
)


################# OLD
#####################

ggplot(PI_SNP [c("chr", "start","PI_SNP_adim")], aes(start, PI_SNP_adim, colour=chr)) + 
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +  #coord_cartesian(ylim = c(0, 0.01)) +
    theme(legend.position = "none",axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks=element_blank()) +
    labs(title = "Pi SNPs") 


