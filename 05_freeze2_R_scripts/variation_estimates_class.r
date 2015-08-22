# LOAD LIBS
#custom libs
source("~/SV_analysis/R_analisis/lib_mike_SV.r")

#graphics
library(ggplot2)
library(reshape2)
library(gridExtra)

#ddply
library(plyr)    



### GET DATA
#dataCl <- getAndFilterFreeze2Class("~/SV_analysis/freeze2.Indels.UNIQUE.crossClass.csv")

#1coord annotations
dataCl_1coord <- getAndFilterFreeze2Class("~/SV_analysis/freeze2.Indels.crossClass.1coord.dsim.csv")

#2coord annotations
dataCl_2coord <- getAndFilterFreeze2Class("~/SV_analysis/freeze2.Indels.crossClass.2coord.dsim.csv")
dataCl_2coord$start <- dataCl_2coord$start +1 #correct from 0based to 1based

#rbind and sort
dataCl <- rbind(dataCl_1coord, dataCl_2coord)
fdataCl <- dataCl[with(dataCl, order(chr, start)), ]
row.names(dataCl)<-NULL #reset row.names



# get rownames for future merging
dataCl$rn <- row.names(dataCl)
#### total 

library(data.table)
dataCl.dt <- data.table(dataCl[c("chr", "start", "class")], keep.rownames=T) # use data.table for efficiency inside loop || drop classInfo column

### GET CONSENSUS ANNOTATIONS.
### Give weight to each annotation CDS[0] > UTR[1] > Intron[2] > Intergenic > [3]
# Get 'classes' vector  
vec <- dataCl.dt[,class]
# recode vec to numbers
vec[vec!="three_prime_UTR" & vec!="five_prime_UTR" & vec!="CDS" & vec!="intron" & vec!="intergenic"]<-99
vec[vec=="CDS"]<-0
vec[vec=="five_prime_UTR" | vec=="three_prime_UTR"]<-1
vec[vec=="intron"]<-2
vec[vec=="intergenic"]<-3   
vec <- as.numeric(vec)
# Add recoded classes to table, order by it.
dataCl.dt[,class:=vec]

### Get consensus annotation
setkey(dataCl.dt, chr, start, class)  # orders by class!! (and the others)

fdataCl.dt <- dataCl.dt[, list(classCode=class[1], 
                               rn=rn[1], 
                               count=length(class), 
                               countCDS=length(class[class==0]),
                               countUTR=length(class[class==1]),
                               countIntron=length(class[class==2]),
                               countIntergenic=length(class[class==3]),
                               countOthers=length(class[class==99])), by=list(chr, start)]

# convert to data.frame
fdataCl.df <- as.data.frame(fdataCl.dt)
# merge with original data.frame, discard empty!
fdataCl <- merge(dataCl, fdataCl.df)
#sort
fdataCl <- fdataCl[with(fdataCl, order(chr, start)), ]
row.names(fdataCl)<-NULL #reset row.names
# remove variables
rm(fdataCl.df, fdataCl.dt, dataCl.dt, vec)

#annotate Small/Long introns
fdataCl$class2[fdataCl$class != "intron"] <- fdataCl$class[fdataCl$class != "intron"]
fdataCl$class2[fdataCl$class == "intron"] <- ifelse(fdataCl$classSize[fdataCl$class == "intron"] <= 100, "small_intron", "long_intron")

#remove rn column
fdataCl <- subset(fdataCl, select=-rn)

# Merge with original data with polarization
#### with Dsim
fdataCl_dsim <- merge(fdata_dsim, fdataCl)
#sort
fdataCl_dsim <- fdataCl_dsim[with(fdataCl_dsim, order(chr, start)), ]
row.names(fdataCl_dsim)<-NULL #reset row.names

# #### with Dsim
# fdataCl_dyak<- merge(fdata_dyak, fdataCl)
# #sort
# fdataCl_dyak <- fdataCl_dyak[with(fdataCl_dyak, order(chr, start)), ]
# row.names(fdataCl_dyak)<-NULL #reset row.names








# GRAPHS pre-dataframes
#remember to change datasets!!!!

newPI <- 2*(fdataCl_dsim$FAC*fdataCl_dsim$RC)/(((fdataCl_dsim$FAC+fdataCl_dsim$RC)*(fdataCl_dsim$FAC+fdataCl_dsim$RC-1)))
fdataCl_dsim$newPI <- newPI
rm(newPI)


fdataCl_graphs <- subset(fdataCl_dsim, blast_stat == "OK")

fdataCl_graphs$class3 <- fdataCl_graphs$class2
fdataCl_graphs$class3[fdataCl_graphs$size %% 3 == 0 & fdataCl_graphs$class2 == "CDS"] <- "CDS_non_frameshift"
fdataCl_graphs$class3[fdataCl_graphs$size %% 3 != 0 & fdataCl_graphs$class2 == "CDS"] <- "CDS_frameshift"
# Annotate X/Autosome in new column
fdataCl_graphs$aut_or_X <- ifelse(fdataCl_graphs$chr=="X", "X", "autosome")
row.names(fdataCl_graphs)<-NULL #reset row.names

# Compare ref IN/DEL to polarized IN/DEL
fdataCl_graphs$ref_in_or_del <- ifelse(fdataCl_graphs$RLN > 1, "DEL", "IN")
fdataCl_graphs$concordant_ref_pol <- ifelse(fdataCl_graphs$in_or_del == fdataCl_graphs$ref_in_or_del , TRUE, FALSE)





### ******Filter JGIL Indels******** (this step had to be done before :( )

# # old.fdataCl_graphs <- fdataCl_graphs

# jgil <- read.delim("../../share/dgrp.freeze2/freeze2.jgil.NO_SNPs_MNPs.csv", header=T, quote="", sep="") # separator is "white space"
# names(jgil) <- c("chr", "start", "id", "refc", "altc", "qual", "cov")
# 
# fdataCl_graphs <- merge(fdataCl_graphs, jgil, all.y=FALSE)
# fdataCl_graphs <- fdataCl_graphs[with(fdataCl_graphs, order(chr, start)), ]
# row.names(fdataCl_graphs)<-NULL #reset row.names

# jgilPI <- 2*(fdataCl_graphs$altc*fdataCl_graphs$refc)/(((fdataCl_graphs$altc+fdataCl_graphs$refc)*(fdataCl_graphs$altc+fdataCl_graphs$refc-1)))
# fdataCl_graphs$jgilPI <- jgilPI
# rm(jgilPI)



# counts
prop.table(
    table(
        subset(fdataCl_graphs, blast_stat == "OK" )$in_or_del,
        subset(fdataCl_graphs, blast_stat == "OK" )$class3
    )
)*100
#     CDS_frameshift CDS_non_frameshift five_prime_UTR intergenic long_intron small_intron three_prime_UTR
# DEL           2205               5328           7610      74035      110667         3659           13069
# IN             911               2053           4339      54755       86383         3664            8917

# DEL      0.5839590          1.4110356      2.0153869 19.6069863  29.3083860    0.9690277       3.4611157
# IN       0.2412638          0.5437042      1.1491148 14.5009865  22.8771567    0.9703518       2.3615249




### distribution by class
ggplot(subset(fdataCl_graphs, size <= 50), aes(x=size, fill=in_or_del )) + 
    geom_histogram(binwidth=1, origin=-0.5) + # origin centers x label tick
    facet_grid(class2 ~ in_or_del , scales="free_y") +
    scale_y_continuous(expand = c(0,0)) + # x axis starts at y zero
    scale_x_continuous(expand = c(0,0)) + # y axis starts at x zero
    theme(#panel.grid.major = element_line(colour = "grey", size = 0.7),
          #panel.grid.minor = element_line(colour = "grey"),
          legend.position="bottom")



### Size by class
ggplot(fdataCl_graphs, aes(class3, size, fill=in_or_del)) + 
    geom_boxplot(notch = TRUE, notchwidth = 0.8) +
    #coord_cartesian(ylim = ylim1*1.05) +
    coord_cartesian(ylim = c(0,25)) +
    facet_grid(. ~ in_or_del ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) # xlabels vertical, centered & aligned

### PI by class
ggplot(fdataCl_graphs, aes(class3, newPI, fill=in_or_del)) + 
    geom_boxplot(notch = TRUE, notchwidth = 0.8) +
    #coord_cartesian(ylim = c(0,25)) +
    facet_grid(. ~ in_or_del ) +
    ylab("PI indels") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) # xlabels vertical, centered & aligned







########### OLD
####################



#compute PI per indel
### manual PI (normal k per indel: FAC*RC/(n over 2)
newPI <- 2*(fdataCl$FAC*fdataCl$RC)/(((fdataCl$FAC+fdataCl$RC)*(fdataCl$FAC+fdataCl$RC-1)))
fdataCl$newPI <- newPI

### (pre) weighted by size PI ( size * k )
wsizePI <- fdataCl$newPI * abs(fdataCl$size)
fdataCl$wsizePI <- wsizePI

rm(newPI,wsizePI)

### Counts
ddply(fdataCl, .(class2), summarize, meanSize = mean(abs(size)), medianSize = median(abs(size)), maxSize = max(abs(size)) , numRow = length(abs(size)))
ddply(fdataCl, .(class2, chr), summarize, meanSize = mean(abs(size)), medianSize = median(abs(size)), maxSize = max(abs(size)) , numRow = length(abs(size)))
ddply(fdataCl, .(class2, isINDEL), summarize, meanSize = mean(abs(size)), medianSize = median(abs(size)), maxSize = max(abs(size)) , numRow = length(abs(size)))
ddply(fdataCl, .(class2, chr, isINDEL), summarize, meanSize = mean(abs(size)), medianSize = median(abs(size)), maxSize = max(abs(size)) , numRow = length(abs(size)))
## alternative --> aggregate(fdataCl["size"], by=fdataCl[c("class","chr")], FUN=length)


### size by class
#compute outliers

fdataCl_graphs <- subset(fdataCl, abs(size) <= 1000)

### Counts
ddply(fdataCl_graphs, .(class2), summarize, meanSize = mean(abs(size)), medianSize = median(abs(size)), maxSize = max(abs(size)) , numRow = length(abs(size)))
ddply(fdataCl_graphs, .(class2, chr), summarize, meanSize = mean(abs(size)), medianSize = median(abs(size)), maxSize = max(abs(size)) , numRow = length(abs(size)))
ddply(fdataCl_graphs, .(class2, isINDEL), summarize, meanSize = mean(abs(size)), medianSize = median(abs(size)), maxSize = max(abs(size)) , numRow = length(abs(size)))
ddply(fdataCl_graphs, .(class2, chr, isINDEL), summarize, meanSize = mean(abs(size)), medianSize = median(abs(size)), maxSize = max(abs(size)) , numRow = length(abs(size)))


### size boxplot
p <- list()
for (i in c("in", "del")) {
    grData <- subset(fdataCl_graphs, isINDEL == i)
    #ylim1 = boxplot.stats(abs(grData$size))$stats[c(1, 5)]
    
    p[[i]] <- ggplot(grData, aes(class2, abs(size))) + 
        geom_boxplot(notch = TRUE, notchwidth = 0.8) +
        #coord_cartesian(ylim = ylim1*1.05) +
        coord_cartesian(ylim = c(0,50)) +
        labs(title = i)
}
do.call(grid.arrange, c(p, nrow=1))  

### Size weighted boxplot
p <- list()
for (i in c("in", "del")) {
    grData <- subset(fdataCl_graphs, isINDEL == i)
    #ylim1 = boxplot.stats(abs(grData$size))$stats[c(1, 5)]
    
    p[[i]] <- ggplot(grData, aes(class2, abs(size)/classSize)) + 
        geom_boxplot(notch = TRUE, notchwidth = 0.8) +
        #coord_cartesian(ylim = ylim1*1.05) +
        coord_cartesian(ylim = c(0,0.15)) +
        labs(title = i)
}
do.call(grid.arrange, c(p, nrow=1))  


### Pi_indel adimensional boxplot
p <- list()
for (i in c("in", "del")) {
    grData <- subset(fdataCl_graphs, isINDEL == i)
    #ylim1 = boxplot.stats(abs(grData$size))$stats[c(1, 5)]
    
    p[[i]] <- ggplot(grData, aes(class2, newPI)) + 
        geom_boxplot(notch = TRUE, notchwidth = 0.8) +
        #coord_cartesian(ylim = ylim1*1.05) +
        #coord_cartesian(ylim = c(0,10)) +
        labs(title = i)
}
do.call(grid.arrange, c(p, nrow=1))  



### distribution by class
ggplot(subset(fdataCl_graphs, abs(size) <= 50), aes(x=size, fill=isINDEL )) + 
    geom_histogram(binwidth=1) + 
    facet_grid(class2 ~ . , scales="free_y") +
    scale_y_continuous(expand = c(0,0)) + #remove blank space below zero
    theme(panel.grid.major = element_line(colour = "grey", size = 0.7),
        panel.grid.minor = element_line(colour = "grey"))


### distribution by class and chr
uniqueClasses <-unique(fdataCl_graphs$class2)
p <- list()
count <- 1
for (i in uniqueClasses) {
    
    grData <- subset(fdataCl_graphs, class2==i)
    preGraph <- ggplot(subset(grData, abs(size) <= 50), aes(x=size, fill=isINDEL )) + 
        geom_histogram(binwidth=1) + 
        facet_grid(class ~ chr, scales="free_y") + 
        theme(legend.position="none", 
              axis.title=element_blank(),
              axis.text.x=element_blank(), 
              axis.ticks.x=element_blank()
              )
    if (count == length(uniqueClasses)) {
        preGraph <- preGraph + theme(axis.text.x=element_text(), axis.ticks.x=element_line() )
    }
    
    p[[i]] <- preGraph
    count <- count+1
}
do.call(grid.arrange, c(p, ncol=1))  

