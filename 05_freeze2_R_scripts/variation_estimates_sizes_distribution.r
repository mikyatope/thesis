# LOAD LIBS
#custom libs
# source("~/SV_analysis/R_analisis/lib_mike_SV.r")
# 
# #graphics
# library(ggplot2)
# library(reshape2)
# library(gridExtra)
# 
# #ddply
# library(plyr)    
# 
# 
# ### GET DATA
# data <- getAndFilterFreeze2("~/SV_analysis/freeze2.Indels.pi.het.vcf")
# #### total 596941 
# fdata <- subset(data, MISS < 102 & RC > 0 & FAC > 0)  # get indels with less than half missing lines and at least 1 line with the reference allele in homozigosis
# #### total 559486
# 
# #get slidding windows
# SWs <- read.delim("~/SV_analysis/slidingWindows_100k.bed", header=F, quote="")
# cols <- c("chr", "start", "end")
# names(SWs) <- cols; rm(cols)
# #testSlide <- SWs[c(1:5),]




#Stuff to do each slide
sumAndMeansSizes <- function(slide, dataset) {
    
    Xchr   <- slide$chr
    Xstart <- slide$start
    Xend   <- slide$end
    
    wsize  <- Xend-Xstart+1
    slide$wsize <- wsize
    
    slideData<- subset(dataset, chr == Xchr & start >= Xstart & start <= Xend, select=c(size, newPI, in_or_del))
    
    slide$numINDELS <- nrow(slideData)
    slide$numIN <- nrow(subset(slideData, in_or_del == "IN"))
    slide$numDEL <- nrow(subset(slideData, in_or_del == "DEL"))
    
    slide$affectedSites <- sum(slideData$size)
    slide$gain <- sum(subset(slideData, in_or_del == "IN")$size)
    slide$loss <- sum(subset(slideData, in_or_del == "DEL")$size)
    
    ###PIs
    slide$PI <- mean(slideData$newPI)
    slide$sPI <- sum(slideData$newPI)
    
    slide$PIindel_adim <- slide$sPI/wsize
    
    ###Sizes
    slide$meanAbsSize <- mean(slideData$size)
    slide$meanInSize <- mean(subset(slideData, in_or_del == "IN")$size)
    slide$meanDelSize <- mean(subset(slideData, in_or_del == "DEL")$size)
    
    slide$medianAbsSize <- median(slideData$size)
    slide$medianInSize <- median(subset(slideData, in_or_del == "IN")$size)
    slide$medianDelSize <- median(subset(slideData, in_or_del == "DEL")$size)
    
    return(slide)
}


system.time( PI_sizes <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), sumAndMeansSizes, dataset=subset(new.fdataCl_dsim, size > 100 & in_or_del != "NA") ))) )







### Plot sizes by chr
ggplot(PI_sizes[c("chr", "start", "meanAbsSize", "meanInSize", "meanDelSize")], aes(x=start, colour=chr)) + 
    #geom_line()+
    #geom_point(aes(size=numINDELS)) +
    geom_line(aes(y=meanInSize), linetype="dashed") + 
    geom_line(aes(y=meanDelSize)) + 
    geom_line(size=1.2, aes(y=meanAbsSize), alpha=0.5) + 
    facet_grid(chr ~ ., scales="free_x") +
    #facet_wrap(~ chr, ncol=2, scales = "free_x") +
    #coord_cartesian(ylim = c(0, 0.001)) +
    ylab("Mean size (bp)") +
    theme(legend.position = "none",
          axis.text.x=element_blank(), 
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(), 
          plot.margin=unit(rep(.5, 4), "lines")) +
    labs(title = "INDEL size") #+ scale_y_log10()


### Plot counts by chr
ggplot(PI_sizes[c("chr", "start",  "numINDELS", "numIN", "numDEL")], aes(x=start, colour=chr)) + 
    #geom_line()+
    #geom_point(aes(size=numINDELS)) +
    geom_line(aes(y=numIN), linetype="dashed") + 
    geom_line(aes(y=numDEL)) + 
    geom_line(size=1.2, aes(y=numINDELS), alpha=0.5) + 
    facet_grid(chr ~ ., scales="free_x") +
    #facet_wrap(~ chr, ncol=2, scales = "free_x") +
    #coord_cartesian(ylim = c(0, 0.001)) +
    ylab("count") +
    theme(legend.position = "bottom",
          axis.text.x=element_blank(), 
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(), 
          plot.margin=unit(rep(.5, 4), "lines")) +
    labs(title = "INDEL count") #+ scale_y_log10()







######################
############# OLD

#Stuff to do each slide
sumAndMeans2 <- function(slide) {
    
    Xchr   <- slide$chr
    Xstart <- slide$start
    Xend   <- slide$end
    
    wsize  <- Xend-Xstart+1
    slide$wsize <- wsize
    
    drops <- c("chr", "start", "REF", "ALT")
    slideData<- subset(fdata, chr == Xchr & start >= Xstart & start <= Xend & abs(size) <= 100 & abs(size) > 1)[,!(names(fdata) %in% drops)]
    
    slide$numINDELS <- nrow(slideData)
    slide$numIN <- nrow(subset(slideData, size > 0))
    slide$numDEL <- nrow(subset(slideData, size < 0))
    
    slide$affectedSites <- sum(abs(slideData$size))
    slide$gain <- sum(subset(slideData, size > 0)$size)
    slide$loss <- sum(subset(slideData, size < 0)$size)
    
    ###PIs
    slide$PI <- mean(slideData$newPI)
    slide$sPI <- sum(slideData$newPI)
    
    slide$sPI_m_plus_gain <- slide$sPI/(wsize+slide$gain)
    slide$sPI_m_plus_abs <- slide$sPI/(wsize+slide$affectedSites)
    slide$sPI_analyzed <- slide$sPI/slide$affectedSites
    
    slide$PI_wsize <- (sum(slideData$wsizePI)*slide$numINDELS)/slide$affectedSites
    slide$PI_wsize2 <- (sum(slideData$wsizePI)*slide$numINDELS)/(wsize+slide$gain)
    slide$PI_wsize3 <- sum(slideData$wsizePI/(wsize+slide$gain))
    slide$PI_wsize4 <- sum(slideData$wsizePI)
    slide$PI_wsize5 <- mean(sum(slideData$wsizePI))
    
    ###Sizes
    slide$meanAbsSize <- mean(abs(slideData$size))
    slide$meanInSize <- mean(subset(slideData, size > 0)$size)
    slide$meanDelSize <- mean(subset(slideData, size < 0)$size)
    
    slide$medianAbsSize <- median(abs(slideData$size))
    slide$medianInSize <- median(subset(slideData, size > 0)$size)
    slide$medianDelSize <- median(subset(slideData, size < 0)$size)
    
    return(slide)
}




## Coordinate slides loop
system.time( more100 <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), sumAndMeans2) )) )
system.time( less100 <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), sumAndMeans2) )) )
system.time( less100no1 <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), sumAndMeans2) )) )

#sizes by window
#big
ggplot(more100[c("chr", "start","numINDELS")], aes(start, numINDELS, colour=chr)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x")

ggplot(more100[c("chr", "start","numIN")], aes(start, numIN, colour=chr)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") 

ggplot(more100[c("chr", "start","numDEL")], aes(start, numDEL, colour=chr)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") 

###small
ggplot(less100[c("chr", "start","numINDELS")], aes(start, numINDELS, colour=chr)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") 

ggplot(less100[c("chr", "start","numIN")], aes(start, numIN, colour=chr)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") 

ggplot(less100[c("chr", "start","numDEL")], aes(start, numDEL, colour=chr)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") 

###small no 1nuc indels
ggplot(less100no1[c("chr", "start","numINDELS")], aes(start, numINDELS, colour=chr)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") 

ggplot(less100no1[c("chr", "start","numIN")], aes(start, numIN, colour=chr)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") 

ggplot(less100no1[c("chr", "start","numDEL")], aes(start, numDEL, colour=chr)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") 


#### Pondered PI ??
ddply(less100, .(chr), summarize, totalChr = sum(numINDELS))

less100$totalChr <- ifelse(less100$chr == "2L", 119821 , 
       ifelse(less100$chr == "2R", 96874, 
       ifelse(less100$chr == "3L", 121565,
       ifelse(less100$chr == "3R", 121438, 
       ifelse(less100$chr == "4", 1465, 
       ifelse(less100$chr == "X", 82446,  NA ))))))

ggplot(less100[c("chr", "start","numINDELS", "sPI_m_plus_gain", "totalChr")], aes(start, (numINDELS*sPI_m_plus_gain)/totalChr, colour=chr)) +  # geom_bar(stat="identity") +
    geom_point(size=1.5) + facet_grid(. ~ chr, scale="free_x", space="free_x") + scale_y_continuous(limits=c(0, 0.00001)) + opts(legend.position = "none") 


