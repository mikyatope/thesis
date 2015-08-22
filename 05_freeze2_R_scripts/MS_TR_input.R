# LOAD LIBS
#custom libs
source("~/SV_analysis/R_analisis/lib_mike_SV.r")

#graphics
library(ggplot2)
library(reshape2)
library(gridExtra)

#ddply
library(plyr)    



input_duplicates_filter <- function(path) {
    ### GET DATA
    data <- getAndFilterFreeze2(path)
    #### rows 596640 
    
    ## remove duplicate/non-biallelic indels (multiple rows with same chr and start/coordinate)
    fdata <- removeDuplicatedChrStart(data) 
    
    # QC filters
    fdata <- subset(fdata, MISS < 102 & RC > 0 & FAC > 0 )  # get indels with less than half missing lines and at least 1 line with the reference allele in homozigosis
    row.names(fdata)<-NULL #reset row.names
    
    #return
    return(fdata)
}

fdata_MS <- input_duplicates_filter("~/SV_analysis/freeze2.MS.ALL.pi.het.vcf")
fdata_TR <- input_duplicates_filter("~/SV_analysis/freeze2.TR.ALL.pi.het.vcf")


# write files to create fastas and do blasts
outdf<-as.data.frame(append(fdata_MS[1:4], list(ID = "."), after = 2)); outdf$QUAL<- NA
write.table(outdf, "~/SV_analysis/blasts/results/MS_indels/MS.semi.vcf", sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
outdf<-as.data.frame(append(fdata_TR[1:4], list(ID = "."), after = 2)); outdf$QUAL<- NA
write.table(outdf, "~/SV_analysis/blasts/results/TR_indels/TR.semi.vcf", sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
rm(outdf)




### Merge again with blast results

### Microsatelites
fdata_MS <- subset(fdata_MS, select=-c(size, REF, ALT)) # remove size column, conflict with blasts dataset. Remove sequences REF/ALT (not really used)
fdata_MS$indelType <- "MS"

#### with Dsim
fdata_MS_dsim <- merge(fdata_MS, blasts_MS_dsim)
#sort
fdata_MS_dsim <- fdata_MS_dsim[with(fdata_MS_dsim, order(chr, start)), ]
row.names(fdata_MS_dsim)<-NULL #reset row.names

#### with Dyak
fdata_MS_dyak <- merge(fdata_MS, blasts_MS_dyak)
#sort
fdata_MS_dyak <- fdata_MS_dyak[with(fdata_MS_dyak, order(chr, start)), ]
row.names(fdata_MS_dyak)<-NULL #reset row.names


### Tandem repeats
fdata_TR <- subset(fdata_TR, select=-c(size, REF, ALT)) # remove size column, conflict with blasts dataset. Remove sequences REF/ALT (not really used)
fdata_TR$indelType <- "TR"

#### with Dsim
fdata_TR_dsim <- merge(fdata_TR, blasts_TR_dsim)
#sort
fdata_TR_dsim <- fdata_TR_dsim[with(fdata_TR_dsim, order(chr, start)), ]
row.names(fdata_TR_dsim)<-NULL #reset row.names

#### with Dyak
fdata_TR_dyak <- merge(fdata_TR, blasts_TR_dyak)
#sort
fdata_TR_dyak <- fdata_TR_dyak[with(fdata_TR_dyak, order(chr, start)), ]
row.names(fdata_TR_dyak)<-NULL #reset row.names
