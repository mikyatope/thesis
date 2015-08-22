# Pi, Theta and Tajima's D by regions from SNP/INDEL frequencies

#IMPORTS

#ddply
library(plyr)    



# FUNCTIONS

##### Global measures/means of population statistics by Regions
# (Stuff to do each region --> to be used with a by loop)
# input: (region) coordinates dataset, (dataset) to be cut -> [must contain columns: chr, start, refc, altc (reference and alternative allele counts), PI]
# returns: dataframe
calcPopStatsRegions <- function(region, dataset, out) {
    
    regionSize  <- region$end-region$start+1
    region$regionSize <- regionSize
    
    # take a slide from of dataset with the variants in that region
    slideData<- subset(dataset, chr == region$chr & start %in% region$start:region$end)
    
    #n <- out.fixNum #number of sequences
    S <- nrow(slideData) #Segregating sites
    
    
    if (!S) {
        warning("no segregating sites")
        region$THETA.w       <- NA
        region$D             <- NA
        region$D.Pval.normal <- NA
        region$D.Pval.beta   <- NA
        region$regionPI      <- NA
        region$S             <- NA
   
        return(region)
    }
    
    region$S <- S # output S
        
    ###Summary PI
    sumPI <- sum(slideData$PI)
    PI <- sumPI/regionSize
    
    region$regionPI <- PI # output PI
    

    ###### Theta and D Tajima
    # adapted from PEGAS 0.5 package
    
    # intermediate values from tajima's and theta formulae
    a1  <- out$a1
    a2  <- out$a2
    b1  <- out$b1
    b2  <- out$b2
    c1  <- out$c1
    c2  <- out$c2
    e1  <- out$e1
    e2  <- out$e2
    
    
    THETA.w <- (S/a1)/regionSize  #Waterson's Theta
    D       <- (sumPI - (S/a1))/sqrt(e1 * S + e2 * S * (S - 1)) #corrected, estimates for region, NOT sites

    # save D, theta outputs
    region$THETA.w <- THETA.w
    region$D       <- D
    
    # P-val
    Dmin <- out$Dmin
    Dmax <- out$Dmax
    tmp1 <- out$tmp1
    tmp2 <- out$tmp2
    a    <- out$a
    b    <- out$b
    
#     n <- 183
#     
#     Dmin <- (2/n - 1/a1)/sqrt(e2)
#     Dmax <- ((n + 1)/(2 * n) - 1/a1)/sqrt(e2)
#     tmp1 <- 1 + Dmin * Dmax
#     tmp2 <- Dmax - Dmin
#     a    <- -tmp1 * Dmax/tmp2
#     b    <- tmp1 * Dmin/tmp2
 
    p    <- pbeta((D - Dmin)/tmp2, b, a) 
    p    <- ifelse(p < 0.5, 2 * p, 2 * (1 - p) )
    
    
    D.Pval.normal <- 2 * pnorm(-abs(D)) 
    D.Pval.beta   <- p
    
    # save P-val's outputs
    region$D.Pval.normal <- D.Pval.normal
    region$D.Pval.beta   <- D.Pval.beta
    
    return(region)
}



calcSumPIRegions <- function(region, dataset) {
    
    # take a slide from of dataset with the variants in that region
    slideData<- subset(dataset, chr == region$chr & start %in% region$start:region$end)
    
    ###Summary PI
    sumPI <- sum(slideData$PI)
    region$sumPI <- sumPI # output sumPI
    return(region)
}














####### MAIN #######


# Estimates by sliding windows

# get slidding windows
SWs <- read.delim("~/SV_analysis/slidingWindows_100k.bed", header=F, quote="")
cols <- c("chr", "start", "end")
names(SWs) <- cols; rm(cols)


#PREPARE DATASET AND REGIONS

# SNPS
fixNum_SNPs <- subset(fSNPdata, select=c(chr,start,refc,altc))
write.table(fixNum_SNPs, "fixedNum/SNPs.csv", sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE) # -> do fix nums outside with Python

#recover fix nums
fixNum_SNPs <- read.delim("fixedNum/SNPs.fixed183.csv", header=F, quote="")
names(fixNum_SNPs) <- c("chr", "start", "foo", "refc", "altc", "fixNum")

fixNum_SNPs$PI <- 2*(fixNum_SNPs$altc*fixNum_SNPs$refc)/(((fixNum_SNPs$altc+fixNum_SNPs$refc)*(fixNum_SNPs$altc+fixNum_SNPs$refc-1)))

fixNum_SNPs_foo <- merge(fixNum_SNPs, subset(fSNPdata, select=c(chr,start,class, class2, classSize)), all.y=FALSE)




# CHROMATIN REGIONS
pre_chromatin_regions<- read.delim("chromatin_regions/YELLOW.bed", header=F, quote="")
names(pre_chromatin_regions) <- c("chr", "start", "end")
pre_chromatin_regions$chromatinType <- "Yellow"
chromatin_regions <- pre_chromatin_regions

pre_chromatin_regions<- read.delim("chromatin_regions/BLUE.bed", header=F, quote="")
names(pre_chromatin_regions) <- c("chr", "start", "end")
pre_chromatin_regions$chromatinType <- "Blue"
chromatin_regions <- rbind(chromatin_regions, pre_chromatin_regions)

pre_chromatin_regions<- read.delim("chromatin_regions/RED.bed", header=F, quote="")
names(pre_chromatin_regions) <- c("chr", "start", "end")
pre_chromatin_regions$chromatinType <- "Red"
chromatin_regions <- rbind(chromatin_regions, pre_chromatin_regions)

pre_chromatin_regions<- read.delim("chromatin_regions/BLACK.bed", header=F, quote="")
names(pre_chromatin_regions) <- c("chr", "start", "end")
pre_chromatin_regions$chromatinType <- "Black"
chromatin_regions <- rbind(chromatin_regions, pre_chromatin_regions)

pre_chromatin_regions<- read.delim("chromatin_regions/GREEN.bed", header=F, quote="")
names(pre_chromatin_regions) <- c("chr", "start", "end")
pre_chromatin_regions$chromatinType <- "Green"
chromatin_regions <- rbind(chromatin_regions, pre_chromatin_regions)

pre_chromatin_regions<- read.delim("chromatin_regions/all.coding.bed", header=F, quote="")
names(pre_chromatin_regions) <- c("chr", "start", "end")
pre_chromatin_regions$chromatinType <- "All coding"
chromatin_regions <- rbind(chromatin_regions, pre_chromatin_regions)

pre_chromatin_regions<- read.delim("chromatin_regions/all.non.coding.bed", header=F, quote="")
names(pre_chromatin_regions) <- c("chr", "start", "end")
pre_chromatin_regions$chromatinType <- "All non coding"
chromatin_regions <- rbind(chromatin_regions, pre_chromatin_regions)

rm(pre_chromatin_regions)
chromatin_regions <- chromatin_regions[with(chromatin_regions, order(chr, start)), ]
row.names(chromatin_regions)<-NULL #reset row.names




chromatin_regions$size <- chromatin_regions$end - chromatin_regions$start

ggplot(chromatin_regions, aes(x=chromatinType, y=size, fill=chromatinType)) + 
    geom_boxplot(notch = TRUE, notchwidth = 0.8) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)#, 
                         #labels = trans_format("log10", math_format(10^.x))
                  ) +
    facet_grid(.~chromatinType, scales="free") + 
    theme_bw() +
    scale_fill_manual(values=c("azure3", "azure4", "grey33", "steelblue3", "olivedrab4", "red3", "gold" )) +
    theme(legend.position = "none",
          axis.title.x=element_blank()) +
    labs(y="Size (bp)")


# Table medians and size
ddply(chromatin_regions, .(chromatinType), summarize, meanSize = mean(size), medianSize = median(size), sumSize=sum(size))










# Precomp Tajima arguments
n <- 183

tmp <- 1:(n - 1)
a1  <- sum(1/tmp) 
a2  <- sum(1/tmp^2)
b1  <- (n + 1)/(3 * (n - 1))
b2  <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
c1  <- b1 - 1/a1
c2  <- b2 - (n + 2)/(a1 * n) + a2/a1^2
e1  <- c1/a1
e2  <- c2/(a1^2 + a2)

Dmin <- (2/n - 1/a1)/sqrt(e2)
Dmax <- ((n + 1)/(2 * n) - 1/a1)/sqrt(e2)
tmp1 <- 1 + Dmin * Dmax
tmp2 <- Dmax - Dmin
a    <- -tmp1 * Dmax/tmp2
b    <- tmp1 * Dmin/tmp2

tajimas.list <- list(tmp=tmp, a1=a1, a2=a2, b1=b1, b2=b2, c1=c1, c2=c2, e1=e1, e2=e2, Dmin=Dmin, Dmax=Dmax, tmp1=tmp1, tmp2=tmp2, a=a, b=b)


system.time( 
    popStats_SNP_SW  <- data.frame(
        do.call(
            "rbind",  
            by(SWs, 1:nrow(SWs), calcPopStatsRegions, dataset=subset(fixNum_SNPs, PI> 0), tajimas.list ) 
        )
    ) 
)


fixNum_SNPs_DT <- data.table(fixNum_SNPs)
system.time( 
    popStats_SNP_chromatin_TEST  <- data.frame(
        do.call(
            "rbind",  
            by(head(chromatin_regions, 10), 1:10, calcPopStatsRegions, dataset=subset(fixNum_SNPs_DT, PI> 0), tajimas.list  ) 
        )
    ) 
)
rm(fixNum_SNPs_DT)


# Table medians D chromatin types
ddply(subset(popStats_SNP_chromatin, D != "NA"), .(chromatinType), summarize, meanD = mean(D), medianD = median(D), sd_D = sd(D), meanDpbeta = mean(D.Pval.beta), medianDpbeta = median(D.Pval.beta), sd_Dpbeta = sd(D.Pval.beta), numRegions = length(D))

# Table medians PI chromatin types
ddply(subset(popStats_SNP_chromatin, regionPI != "NA"), .(chromatinType), summarize, meanPI = mean(regionPI), medianPI = median(regionPI), sd_PI = sd(regionPI), numRegions = length(regionPI))

# Table medians count S
ddply(subset(popStats_SNP_chromatin, regionPI != "NA"), .(chromatinType), summarize, totalS = sum(S))




# Theta/D correction

popStats_SNP_chromatin$regionSize <- chromatin_regions$size

preDF <- popStats_SNP_chromatin

preDF$sumPI <- popStats_SNP_sumPI$sumPI

preDF$THETA.w.sites <- popStats_SNP_chromatin$THETA.w/popStats_SNP_chromatin$regionSize

preDF$D <- (preDF$sumPI - preDF$THETA.w)/sqrt(e1 * preDF$S + e2 * preDF$S * (preDF$S - 1))


p    <- pbeta((preDF$D - Dmin)/tmp2, b, a) 
p    <- ifelse(p < 0.5, 2 * p, 2 * (1 - p) )

preDF$D.Pval.normal <- 2 * pnorm(-abs(preDF$D)) 
preDF$D.Pval.beta   <- p



# Table medians D chromatin types
ddply(subset(preDF, D != "NA"), .(chromatinType), summarize, meanD = mean(D), medianD = median(D), sd_D = sd(D), numRegions = length(D))






# (PI - THETA.w)/sqrt(e1 * S + e2 * S * (S - 1))

# n <- 183
# 
# tmp <- 1:(n - 1)
# a1  <- sum(1/tmp) 
# a2  <- sum(1/tmp^2)
# b1  <- (n + 1)/(3 * (n - 1))
# b2  <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
# c1  <- b1 - 1/a1
# c2  <- b2 - (n + 2)/(a1 * n) + a2/a1^2
# e1  <- c1/a1
# e2  <- c2/(a1^2 + a2)
# 
# Dmin <- (2/n - 1/a1)/sqrt(e2)
# Dmax <- ((n + 1)/(2 * n) - 1/a1)/sqrt(e2)
# tmp1 <- 1 + Dmin * Dmax
# tmp2 <- Dmax - Dmin
# a    <- -tmp1 * Dmax/tmp2
# b    <- tmp1 * Dmin/tmp2
# 
# 
# (0.8649493 - 0.8644586)/sqrt(e1 * 5 + e2 * 5 * (5 - 1))
# 
# 
# # 2*(altc*refc)/(((altc+refc)*(altc+refc-1)))
# 2*(2*3)/(((2+3)*(2+3-1)))

pbeta((-0.4779254 - Dmin)/tmp2, b, a) 

ifelse( 0.388157 < 0.5, 2 *  0.388157, 2 * (1 -  0.388157) )





# Global PIs

fixNum_SNPs_DT <- data.table(fixNum_SNPs)
system.time( 
    popStats_SNP_sumPI  <- data.frame(
        do.call(
            "rbind",  
            by(chromatin_regions, 1:nrow(chromatin_regions), calcSumPIRegions, dataset=subset(fixNum_SNPs_DT, PI> 0) ) 
        )
    ) 
)
rm(fixNum_SNPs_DT)


ddply(subset(popStats_SNP_sumPI, sumPI != "NA" & (chromatinType=="All coding" | chromatinType=="All non coding")), .(), summarize,  sumsumPI = sum(sumPI), numRegions = length(sumPI))
ddply(subset(popStats_SNP_sumPI, sumPI != "NA" & (chromatinType!="All coding" & chromatinType!="All non coding")), .(), summarize,  sumsumPI = sum(sumPI), numRegions = length(sumPI))




#popStats_SNP_chromatin_corrected <- preDF 
#rm(preDF)



# Boxplot D chromatin Types

ggplot(preDF, aes(x=chromatinType, y=D, fill=chromatinType)) + 
    geom_boxplot(notch = TRUE, notchwidth = 0.8) + #, outlier.shape = NA) + 
    facet_grid(.~chromatinType, scales="free") + 
#     coord_cartesian(ylim = c(0.001,-0.001)) +
#     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
#                   labels = trans_format("log10", math_format(10^.x))
#     ) +
    theme_bw() +
    scale_fill_manual(values=c("azure3", 
                               "azure4", 
                               "grey33", 
                               "steelblue3", 
                               "olivedrab4", 
                               "red3", 
                               "gold" )) +
    theme(legend.position = "none",
          axis.title.x=element_blank()) +
    labs(y="Tajima's D")



ggplot(preDF, aes(x=D, fill=chromatinType)) + 
    geom_histogram(binwidth=0.00005) + 
    facet_grid(.~chromatinType, scales="free") + 
    coord_cartesian(ylim = c(0,500)) +
    #     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)#, 
    #                   #labels = trans_format("log10", math_format(10^.x))
    #     ) +
    theme_bw() +
    scale_fill_manual(values=c("azure3", 
                               "azure4", 
                               "grey33", 
                               "steelblue3", 
                               "olivedrab4", 
                               "red3", 
                               "gold" )) +
    theme(legend.position = "none")






#correlations

for (chromatinTypeText in unique(chromatin_regions$chromatinType)){
    
#     sliceDF_size <- subset(chromatin_regions, chromatinType==chromatinTypeText)
    sliceDF_D <- subset(popStats_SNP_chromatin, chromatinType==chromatinTypeText)
    
    cat(chromatinTypeText, "correlation: Size/D ")
    print(cor.test(~ sliceDF_size$size + sliceDF_D$D, method = "spearman"))
    
    cat(chromatinTypeText, "correlation: Size/PI ")
    print(cor.test(~ sliceDF_size$size + sliceDF_D$regionPI, method = "spearman"))
}












######### TRASH / TESTS #########
#################################


# ### Individual Population Statistics 
# 
# # PI calculations simplified for bialelic alleles
# # input: dataframe (df). MUST have columns named 'refc' and 'altc' (reference/alternative alleles count)
# # returns: vector
# bialelicPI <- function(df) {
#  PI_vec <- 2*(df$altc*df$refc)/(((df$altc+df$refc)*(df$altc+df$refc-1)))
#  return(PI_vec)
# }
# 
# 


#temp_df <- subset(fSNPdata, select=c(chr,start,refc, altc))
#names(temp_df)[names(temp_df) == 'jgilPI'] <- 'PI'

# fix number segregating sites (95% of SNPs have 183 alleles or more)

# 
# fixNum <- function(df, fixnum) {
#     if (df$refc+df$altc >= fixnum){ 
#         countsDF    <- as.data.frame(table(sample(rep(c("r","a"), c(df$refc, df$altc)), fixnum)))
#         df$fix_refc <- countsDF[countsDF$Var1=="r",]$Freq
#         df$fix_altc <- countsDF[countsDF$Var1=="a",]$Freq
#         df$fix_PI   <- 2*(df$fix_altc*df$fix_refc)/(((df$fix_altc+df$fix_refc)*(df$fix_altc+df$fix_refc-1)))
#     } else {
#         df$fix_refc <- NA
#         df$fix_altc <- NA
#         df$fix_PI   <- NA
#     }
#     return(df)
# }
# 
# fixNum <- function(refc, altc) {
#     if((refc+altc) >= 183){ 
#         return( subset(as.data.frame(table(sample(rep(c("r","a"), c(refc, altc)), 183))), Var1=="r")$Freq )
#     } else {
#         return(NA)
#     }
# }
# 
# fixNum <- function(x) {
#     y <- sample([0:x],1)
#     y
# }
# v1<-pmin(temp_df1$refc,temp_df1$altc)
# aaply(v1, 1, function(x) sample(c(0:x),1))
# 
# 
# ifelse((temp_df1$refc+temp_df1$altc) >= 183, min(temp_df1$refc,temp_df1$altc), NA )
# ifelse((temp_df1$refc+temp_df1$altc) >= 183, subset(as.data.frame(table(sample(rep(c("r","a"), c(temp_df1$refc,temp_df1$altc)), 183))), Var1=="r")$Freq, NA )
# 
# subset(as.data.frame(table(sample(rep(c("r","a"), c(temp_df1$refc,temp_df1$altc)), 183))), Var1=="r")$Freq




# temp_df2 <- data.frame(do.call("rbind", by(temp_df1 , 1:nrow(temp_df1), fixNum, fixnum=183) ) ) 
# 
# temp_df2 <- adply(head(temp_df1), 1, summarise, foo = fixNum(refc, altc, 183)) 

# library(data.table)
# temp_DT1 <- data.table(head(temp_df1, 30))
# temp_DT1[, foo := floor(runif(1, 0, pmin(refc,altc)))]
# temp_DT1[, foo := sample(0:min(refc,altc), 1)]


p    <- pbeta((-1.301523e-04 - Dmin)/tmp2, b, a) 
ifelse(p < 0.5, 2 * p, 2 * (1 - p) )




##### Test Trudy Means/Complete

test_chromatin_regions_meanscomplete <- read.delim("testTajimasDmeanscomplete.csv", header=F, quote="")
names(test_chromatin_regions_meanscomplete) <- c("chr", "start", "end", "chromatinType", "size")
test_chromatin_regions_meanscomplete <- transform(test_chromatin_regions_meanscomplete, chr = as.character(chr))

fixNum_SNPs_DT <- data.table(fixNum_SNPs)
system.time( 
    popStats_SNP_chromatin_TEST_meanscomplete  <- data.frame(
        do.call(
            "rbind",  
            by(test_chromatin_regions_meanscomplete, 1:nrow(test_chromatin_regions_meanscomplete), calcPopStatsRegions, dataset=subset(fixNum_SNPs_DT, PI> 0), tajimas.list  ) 
        )
    ) 
)
rm(fixNum_SNPs_DT)