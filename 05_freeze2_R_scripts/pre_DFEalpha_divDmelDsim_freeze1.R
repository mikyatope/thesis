# Pre DFE-alpha computations


#### Annotate GAPS between Dmel-Dsim
DIV_DmelDsim_freeze1SNPS <- read.delim("~/SV_analysis/selection/divergence.dsim-dmel.ALL.csv", header=F, quote="")
names(DIV_DmelDsim_freeze1SNPS) <- c("chr", "start")
#sort
DIV_DmelDsim_freeze1SNPS <- DIV_DmelDsim_freeze1SNPS[with(DIV_DmelDsim_freeze1SNPS, order(chr, start)), ]
row.names(DIV_DmelDsim_freeze1SNPS)<-NULL #reset row.names


# remove duplicates (remove segregating sites not detected in freeze 1) --> Literally, display rows in df1 not present in df2 by chr AND start
tempSNPs <- subset(fSNPdata, select=c(chr,start))
DIV_DmelDsim_freeze1SNPS <- DIV_DmelDsim_freeze1SNPS[!(DIV_DmelDsim_freeze1SNPS$chr %in% tempSNPs$chr & DIV_DmelDsim_freeze1SNPS$start %in% tempSNPs$start),]
rm(tempSNPs)

### output
outdf <- DIV_DmelDsim_freeze1SNPS
outdf$end <- outdf$start
outdf$start <- outdf$start-1
write.table(outdf, file="selection/divergence.dsim-dmel.ALL.NOSEG.csv", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
rm(outdf,DIV_DmelDsim_freeze1SNPS)


#######
## DO ANNOTATION OUTSIDE WITH INTERSECT BED
#######
# recover file
DIV_DmelDsim_freeze1SNPSCl <- read.delim("~/SV_analysis/selection/divergence.MERGED.dmelAnnotations.csv", header=F, quote="")
names(DIV_DmelDsim_freeze1SNPSCl) <- c("chr", "start", "class", "classStart", "classEnd")

# Create classSize
DIV_DmelDsim_freeze1SNPSCl$classSize <- (DIV_DmelDsim_freeze1SNPSCl$classEnd-DIV_DmelDsim_freeze1SNPSCl$classStart)+1
#DIV_DmelDsim_freeze1SNPSCl$class <- as.character(DIV_DmelDsim_freeze1SNPSCl$class) # make sure is not factor!

#sort
DIV_DmelDsim_freeze1SNPSCl <- DIV_DmelDsim_freeze1SNPSCl[with(DIV_DmelDsim_freeze1SNPSCl, order(chr, start)), ]
row.names(DIV_DmelDsim_freeze1SNPSCl)<-NULL #reset row.names

#### OLD ANNOTATION CONSENSUS CODE

# get rownames for future merging
DIV_DmelDsim_freeze1SNPSCl$rn <- row.names(DIV_DmelDsim_freeze1SNPSCl)
#### total 

library(data.table)
DIV_DmelDsim_freeze1SNPSCl.dt <- data.table(DIV_DmelDsim_freeze1SNPSCl[c("chr", "start", "class")], keep.rownames=T) # use data.table for efficiency inside loop || drop classInfo column

### GET CONSENSUS ANNOTATIONS.
### Give weight to each annotation CDS[0] > UTR[1] > Intron[2] > Intergenic > [3]
# Get 'classes' vector  
vec <- DIV_DmelDsim_freeze1SNPSCl.dt[,class]
# recode vec to numbers
vec[vec!="three_prime_UTR" & vec!="five_prime_UTR" & vec!="CDS" & vec!="intron" & vec!="intergenic"]<-99
vec[vec=="CDS"]<-0
vec[vec=="five_prime_UTR" | vec=="three_prime_UTR"]<-1
vec[vec=="intron"]<-2
vec[vec=="intergenic"]<-3   
vec <- as.numeric(vec)
# Add recoded classes to table, order by it.
DIV_DmelDsim_freeze1SNPSCl.dt[,class:=vec]

### Get consensus annotation
setkey(DIV_DmelDsim_freeze1SNPSCl.dt, chr, start, class)  # orders by class!! (and the others)

fDIV_DmelDsim_freeze1SNPSCl.dt <- DIV_DmelDsim_freeze1SNPSCl.dt[, list(classCode=class[1], 
                               rn=rn[1], 
                               count=length(class), 
                               countCDS=length(class[class==0]),
                               countUTR=length(class[class==1]),
                               countIntron=length(class[class==2]),
                               countIntergenic=length(class[class==3]),
                               countOthers=length(class[class==99])), by=list(chr, start)]

# convert to data.frame
fDIV_DmelDsim_freeze1SNPSCl.df <- as.data.frame(fDIV_DmelDsim_freeze1SNPSCl.dt)
# merge with original data.frame, discard empty!
fDIV_DmelDsim_freeze1SNPSCl <- merge(DIV_DmelDsim_freeze1SNPSCl, fDIV_DmelDsim_freeze1SNPSCl.df)
#sort
fDIV_DmelDsim_freeze1SNPSCl <- fDIV_DmelDsim_freeze1SNPSCl[with(fDIV_DmelDsim_freeze1SNPSCl, order(chr, start)), ]
row.names(fDIV_DmelDsim_freeze1SNPSCl)<-NULL #reset row.names
# remove variables
rm(fDIV_DmelDsim_freeze1SNPSCl.df, fDIV_DmelDsim_freeze1SNPSCl.dt, DIV_DmelDsim_freeze1SNPSCl.dt, vec)

#annotate Small/Long introns
fDIV_DmelDsim_freeze1SNPSCl$class2[fDIV_DmelDsim_freeze1SNPSCl$class != "intron"] <- fDIV_DmelDsim_freeze1SNPSCl$class[fDIV_DmelDsim_freeze1SNPSCl$class != "intron"] # the same if annotation is not 'intron'
fDIV_DmelDsim_freeze1SNPSCl$class2[fDIV_DmelDsim_freeze1SNPSCl$class == "intron"] <- ifelse(fDIV_DmelDsim_freeze1SNPSCl$classSize[fDIV_DmelDsim_freeze1SNPSCl$class == "intron"] <= 100, "small_intron", "long_intron")

#remove rn column
fDIV_DmelDsim_freeze1SNPSCl <- subset(fDIV_DmelDsim_freeze1SNPSCl, select=-rn)

#remove unfiltered dataset
rm(DIV_DmelDsim_freeze1SNPSCl)








