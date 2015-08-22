# LOAD LIBS
#custom libs
source("~/SV_analysis/R_analisis/lib_mike_SV.r")

#graphics
library(ggplot2)
library(reshape2)
library(gridExtra)

#ddply
library(plyr)    

### Check Working Directory! Rstudio changes it randomly ¬¬
##setwd('~/SV_analysis')

### GET DATA
# INDELS from other script


#get slidding windows
SWs <- read.delim("~/SV_analysis/slidingWindows_100k.bed", header=F, quote="")
cols <- c("chr", "start", "end")
names(SWs) <- cols; rm(cols)
#testSlide <- SWs[c(1:5),]


calculatePI <- function(slide, dataset) {
    Xchr   <- slide$chr
    Xstart <- slide$start
    Xend   <- slide$end
    
    wsize  <- Xend-Xstart+1
    slide$wsize <- wsize

    slideData<- subset(dataset, chr == Xchr & start >= Xstart & start <= Xend) #[,!(names(dataset) %in% drops)]
    
    slide$numINDELS <- nrow(slideData)
    
    ###PIs
    slide$PI  <- mean(slideData$newPI)
    slide$sPI <- sum(slideData$newPI)
    
    slide$PIindel_adim <- slide$sPI/wsize
    
    return(slide)
}

## Coordinate slides loop
system.time( 
    PI <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), calculatePI, dataset=subset(fdataCl_graphs, in_or_del != "NA") ))) 
)
system.time( 
    PI_DEL  <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), calculatePI, dataset=subset(fdataCl_graphs, in_or_del == "DEL") ))) 
)
system.time( 
    PI_IN  <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), calculatePI, dataset=subset(fdataCl_graphs, in_or_del == "IN") ))) 
)






### GRAPHS
##########

SW_pi_plot <- ggplot(PI [c("chr", "start","PIindel_adim", "numINDELS")], aes(start, PIindel_adim, colour=chr)) + 
    #geom_line()+
    #geom_point(aes(size=numINDELS)) +
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +
    #facet_wrap(~ chr, ncol=2, scales = "free_x") +
    coord_cartesian(ylim = c(0, 0.001)) +
    ylab("PI indel") +
    theme(legend.position = "none",axis.text.x=element_blank(), 
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(), 
          plot.margin=unit(rep(.5, 4), "lines")) +
    labs(title = "PI indel TOTAL") #+ scale_y_log10()

SW_pi_plot_DEL <- ggplot(PI_DEL [c("chr", "start","PIindel_adim", "numINDELS")], aes(start, PIindel_adim, colour=chr)) + 
    #geom_line()+
    #geom_point(aes(size=numINDELS)) +
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +
    #facet_wrap(~ chr, ncol=2, scales = "free_x") +
    coord_cartesian(ylim = c(0, 0.001)) +
    ylab("PI indel") +
    theme(legend.position = "none",axis.text.x=element_blank(), 
            axis.title.x=element_blank(), 
            axis.ticks=element_blank(), 
            plot.margin=unit(rep(.5, 4), "lines")) + 
    labs(title = "PI DELETIONS") #+ scale_y_log10()

SW_pi_plot_IN <- ggplot(PI_IN [c("chr", "start","PIindel_adim", "numINDELS")], aes(start, PIindel_adim, colour=chr)) + 
    #geom_line()+
    #geom_point(aes(size=numINDELS)) +
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +
    #facet_wrap(~ chr, ncol=2, scales = "free_x") +
    coord_cartesian(ylim = c(0, 0.001)) +
    ylab("PI indel") +
    theme(legend.position = "none",axis.text.x=element_blank(), 
            axis.title.x=element_blank(), 
            axis.ticks=element_blank(), 
            plot.margin=unit(rep(.5, 4), "lines")) +
    labs(title = "PI INSERTIONS") #+ scale_y_log10()

grid.arrange(SW_pi_plot, SW_pi_plot_DEL, SW_pi_plot_IN)


# VARIATION by size
ggplot(fdataCl_graphs, aes(factor(cut(size, c(0,1,2,3,4,5,10,20,50,100,1000,10000,Inf))), newPI, fill=in_or_del )) + 
    geom_boxplot(notch = TRUE, notchwidth = 0.8) + 
    facet_grid(. ~ in_or_del) +
    ylab("log PI") + xlab("INDEL size") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + # xlabels vertical, centered & aligned
    scale_y_log10() 







###################### 
#### OLD


## Coordinate slides loop
system.time( 
        PI <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), sumAndMeans, dataset=fdata ))) 
    )
system.time( 
    PI_DEL  <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), sumAndMeans, dataset=subset(fdata, size < 0) ))) 
)
system.time( 
    PI_IN  <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), sumAndMeans, dataset=subset(fdata, size > 0) ))) 
)


ggplot(PI [c("chr", "start","PIindel_adim", "numINDELS")], aes(start, PIindel_adim, colour=chr)) + 
    #geom_line()+
    #geom_point(aes(size=numINDELS)) +
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +
    #facet_wrap(~ chr, ncol=2, scales = "free_x") +
    #coord_cartesian(ylim = c(0, 0.01)) +
    theme(legend.position = "none",axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks=element_blank()) +
    labs(title = "PI INDEL adimensional") #+ scale_y_log10()

ggplot(PI_DEL [c("chr", "start","PIindel_adim", "numINDELS")], aes(start, PIindel_adim, colour=chr)) + 
    #geom_line()+
    #geom_point(aes(size=numINDELS)) +
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +
    #facet_wrap(~ chr, ncol=2, scales = "free_x") +
    #coord_cartesian(ylim = c(0, 0.01)) +
    theme(legend.position = "none",axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks=element_blank()) +
    labs(title = "PI DELETION adimensional") #+ scale_y_log10()

ggplot(PI_IN [c("chr", "start","PIindel_adim", "numINDELS")], aes(start, PIindel_adim, colour=chr)) + 
    #geom_line()+
    #geom_point(aes(size=numINDELS)) +
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +
    #facet_wrap(~ chr, ncol=2, scales = "free_x") +
    #coord_cartesian(ylim = c(0, 0.01)) +
    theme(legend.position = "none",axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks=element_blank()) +
    labs(title = "PI INSERTION adimensional") #+ scale_y_log10()



# VARIATION by size
ggplot(fdata, aes(factor(cut(abs(size), c(-Inf,1,2,3,4,5,10,20,50,100,1000,10000,Inf))), newPI )) + 
    geom_boxplot(notch = TRUE, notchwidth = 0.8) + 
    scale_y_log10() 



#SIZES
#mean
ggplot(PI [c("chr", "start","meanAbsSize")], aes(start, meanAbsSize, colour=chr)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") + scale_y_log10()

ggplot(PI [c("chr", "start","meanInSize")], aes(start, meanInSize)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") + scale_y_log10()

ggplot(PI [c("chr", "start","meanDelSize")], aes(start, abs(meanDelSize))) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") + scale_y_log10()


#median
ggplot(PI [c("chr", "start","medianAbsSize")], aes(start, medianAbsSize)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") 

ggplot(PI [c("chr", "start","medianInSize")], aes(start, medianInSize)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") 

ggplot(PI [c("chr", "start","medianDelSize")], aes(start, abs(medianDelSize))) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") 




#numbers
ggplot(PI [c("chr", "start","numINDELS")], aes(start, numINDELS, colour=chr)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") + scale_x_continuous(labels = "comma")

ggplot(PI [c("chr", "start","numIN")], aes(start, numIN, colour=chr)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") 

ggplot(PI [c("chr", "start","numDEL")], aes(start, numDEL, colour=chr)) +  # geom_bar(stat="identity") +
    geom_line() + facet_wrap(~ chr, ncol=2, scale="free_x") 



#### TEST multiple variables
df.m2 <- melt(PI [c("chr", "start","meanAbsSize", "meanInSize", "meanDelSize")], id=c("chr", "start"))
ggplot(df.m2, aes(start, abs(value), colour=variable)) + 
    geom_bar(stat="identity") + 
    facet_wrap(variable ~ chr, nrow = 3, scales = "free_x") # + coord_cartesian(ylim = c(0, 0.01)) +





