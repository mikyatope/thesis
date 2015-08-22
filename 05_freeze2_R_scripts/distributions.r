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
data <- getAndFilterFreeze2("~/SV_analysis/freeze2.Indels.UNIQUE2.pi.het.vcf")
#### rows 596640 

## remove duplicate/non-biallelic indels (multiple rows with same chr and start/coordinate)
fdata <- removeDuplicatedChrStart(data) 
#### rows 590674 

# QC filters
fdata <- subset(fdata, MISS < 102 & RC > 0 & FAC > 0 )  # get indels with less than half missing lines and at least 1 line with the reference allele in homozigosis
row.names(fdata)<-NULL #reset row.names
#### rows 559314


# Output to file simulating BED to do blasts
outdf <- fdata[c("chr","start")]
outdf$end <- outdf$start
write.table(outdf, "dmel_anotations/filteredINDELS.bed", sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
rm(outdf)




## Merge with polarizations
fdata <- subset(fdata, select=-c(REF, ALT)) # Remove sequences REF/ALT (not really used)
fdata$indelType <- "noClass"

#### with Dsim
fdata_dsim <- merge(fdata, blasts_dsim)
#merge with MS and TR
fdata_dsim <- rbind(fdata_dsim, fdata_MS_dsim, fdata_TR_dsim)
#sort
fdata_dsim <- fdata_dsim[with(fdata_dsim, order(chr, start)), ]
row.names(fdata_dsim)<-NULL #reset row.names


### ******Filter JGIL Indels******** (this step had to be done before :( )

# # old.fdata_dsim <- fdata_dsim

# jgil <- read.delim("../../share/dgrp.freeze2/freeze2.jgil.NO_SNPs_MNPs.csv", header=T, quote="", sep="") # separator is "white space"
# names(jgil) <- c("chr", "start", "id", "refc", "altc", "qual", "cov")
# 
# fdata_dsim <- merge(fdata_dsim, jgil, all.y=FALSE)
# fdata_dsim <- fdata_dsim[with(fdata_dsim, order(chr, start)), ]
# row.names(fdata_dsim)<-NULL #reset row.names




#OUTPUT to file (BED start is zero based for BEDtools!!)

# REF length (RLN) == 1
outdf <- subset(fdata_dsim, blast_stat=="OK" & RLN == 1, select=c(chr,start))
outdf$end <- outdf$start
outdf$start <- outdf$start-1
write.table(outdf, "dmel_anotations/polarized_INDELS_1cord_dsim.bed", sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
rm(outdf)


# REF length (RLN) > 1
outdf <- subset(fdata_dsim, blast_stat=="OK" & RLN > 1, select=c(chr,start, RLN))
outdf$end <- (outdf$start + outdf$RLN) - 1
outdf$start <- outdf$start-1
outdf <- subset(outdf, select=-c(RLN))
write.table(outdf, "dmel_anotations/polarized_INDELS_2cord_dsim.bed", sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
rm(outdf)




#### with Dyak
fdata_dyak <- merge(fdata, blasts_dyak)
#merge with MS and TR
fdata_dyak <- rbind(fdata_dyak, fdata_MS_dyak, fdata_TR_dyak)
#sort
fdata_dyak <- fdata_dyak[with(fdata_dyak, order(chr, start)), ]
row.names(fdata_dyak)<-NULL #reset row.names








###### size distribution
## fdata_dsim$size    <- abs(fdata_dsim$ALN-fdata_dsim$RLN)
## fdata_dyak$size    <- abs(fdata_dyak$ALN-fdata_dyak$RLN)

less50 <- subset(fdata_dsim, size <= 50 & blast_stat == "OK" )
ggplot(less50, aes(x=size, fill=in_or_del)) + geom_histogram(binwidth=1, origin=-0.5) + 
    facet_grid(chr ~ in_or_del, scales="free_y") + 
    scale_y_continuous(expand = c(0,0)) + # x axis starts at y zero
    scale_x_continuous(expand = c(0,0)) + # y axis starts at x zero
    theme(legend.position="bottom")


### Counts
prop.table(table(subset(fdata_dsim, blast_stat == "OK" )$in_or_del))*100
# DEL       IN 
# 223866    185328
# 54.70901  45.29099

# summary(as.factor(subset(fdata_dsim, blast_stat == "OK")$indelType))
# MS noClass      TR 
# 31415  348842   28937 




prop.table(table(subset(fdata_dyak, blast_stat == "OK" )$in_or_del))*100
# DEL       IN 
# 138103    108122 
# 56.08813  43.91187

# summary(as.factor(subset(fdata_dyak, blast_stat == "OK")$indelType))
# MS noClass      TR 
# 23802  202394   20029 



prop.table(
    table(
        subset(fdata_dsim, blast_stat == "OK" )$in_or_del, 
        cut( subset(fdata_dsim, blast_stat == "OK" )$size, c(0,1,2,3,4,5,6,7,8,9,10,20,50,100,1000,10000,Inf))
    )
)*100 ## table counts, cut create bins, prop.table count frequency

#     (0,1] (1,2] (2,3] (3,4] (4,5] (5,6] (6,7] (7,8] (8,9] (9,10] (10,20] (20,50] (50,100] (100,1e+03]
# DEL 37399 20103 19320 16671 13351 17909 13001 12297 10972   8859   42927   11020       25          12
# IN  37736 33336 21053 20026  9902 14636  8474  9072  5928   4757   15983    3382      309         686
# 
#     (1e+03,1e+04] (1e+04,Inf]
# DEL             0           0
# IN             48           0

#            (0,1]        (1,2]        (2,3]        (3,4]        (4,5]        (5,6]        (6,7]        (7,8]
# DEL  9.139674580  4.912828634  4.721476854  4.074106658  3.262755563  4.376652639  3.177221562  3.005176029
# IN   9.222031604  8.146747020  5.144992351  4.894011154  2.419879079  3.576787538  2.070900355  2.217041306
# 
#            (8,9]       (9,10]      (10,20]      (20,50]     (50,100]  (100,1e+03] (1e+03,1e+04]  (1e+04,Inf]
# DEL  2.681368740  2.164987756 10.490623030  2.693099117  0.006109571  0.002932594   0.000000000  0.000000000
# IN   1.448701594  1.162529265  3.905971251  0.826502832  0.075514304  0.167646642   0.011730377  0.000000000

prop.table(
    table(
        subset(fdata_dsim, blast_stat == "OK" )$in_or_del, 
        cut( subset(fdata_dsim, blast_stat == "OK" )$size, c(0,100,Inf))
    )
)*100 ## table counts, cut create bins, prop.table count frequency

#     (0,100] (100,Inf]
# DEL  223854        12
# IN   184594       734

#          (0,100]    (100,Inf]
# DEL 54.706080734  0.002932594
# IN  45.111609652  0.179377019


prop.table(
    table(
        subset(fdata_dyak, blast_stat == "OK" )$in_or_del, 
        cut( subset(fdata_dyak, blast_stat == "OK" )$size, c(0,1,2,3,4,5,6,7,8,9,10,20,50,100,1000,10000,Inf))
    )
)*100 ## table counts, cut create bins, prop.table count frequency

#           (0,1]       (1,2]       (2,3]       (3,4]       (4,5]       (5,6]       (6,7]       (7,8]       (8,9]
# DEL 9.466544827 5.994517210 5.264696924 4.506041222 3.294953802 5.064067418 3.314448167 3.182455072 2.718651640
# IN  8.719666971 7.745760991 5.038481064 4.519849731 2.270281247 3.536602701 2.009137984 2.153721190 1.456797644
# 
#          (9,10]     (10,20]     (20,50]    (50,100] (100,1e+03] (1e+03,1e+04] (1e+04,Inf]
# DEL 2.143161742 9.286628084 1.835313230 0.008934917 0.004061326   0.003655193 0.000000000
# IN  1.123768911 3.584526348 0.955630013 0.184790334 0.547466748   0.065387349 0.000000000

prop.table(
    table(
        subset(fdata_dyak, blast_stat == "OK" )$in_or_del, 
        cut( subset(fdata_dyak, blast_stat == "OK" )$size, c(0,100,Inf))
    )
)*100 ## table counts, cut create bins, prop.table count frequency

#          (0,100]    (100,Inf]
# DEL 56.080414255  0.007716519
# IN  43.299015128  0.612854097



# Ranges of sizes
# summary(subset(fdata_dsim, blast_stat == "OK" & in_or_del == "DEL" )$size)
# 
# summary(subset(fdata_dsim, blast_stat == "OK" & in_or_del == "IN" )$size)
# 
# summary(subset(fdata_dyak, blast_stat == "OK" & in_or_del == "DEL" )$size)
# 
# summary(subset(fdata_dyak, blast_stat == "OK" & in_or_del == "IN" )$size)




### Counts size by Chr
ddply(subset(fdata_dsim, indelType=="noClass"), .(chr), summarize, meanSize = mean(size), medianSize = median(size), maxSize = max(size) , sd = sd(size), numRow = length(size))

###### POST jgil
#   chr meanSize medianSize maxSize       sd numRow
# 1  2L 6.241390          4    2945 17.70288  51510
# 2  2R 6.527628          4    4698 40.56623  42185
# 3  3L 6.451567          4    4110 33.02855  51917
# 4  3R 6.580835          4    6389 36.35092  55768
# 5   4 7.206897          4      81 10.37105    232
# 6   X 6.061766          4     671 10.90065  35861


#    chr in_or_del meanSize medianSize maxSize        sd numRow
# 1   2L       DEL 6.840894          5      49  6.003952  31036
# 2   2L        IN 5.332617          3    2945 27.064050  20474
# 3   2R       DEL 6.897422          5     187  6.192216  24830
# 4   2R        IN 5.998559          3    4698 62.807851  17355
# 5   3L       DEL 6.917557          5      75  6.088600  31355
# 6   3L        IN 5.740979          3    4110 51.933514  20562
# 7   3R       DEL 7.022426          5     370  6.407751  34246
# 8   3R        IN 5.878171          3    6389 57.947868  21522
# 9    4       DEL 6.648148          5      32  6.355388    162
# 10   4        IN 8.500000          3      81 16.230004     70
# 11   X       DEL 6.849546          5     511  7.076109  18298
# 12   X        IN 5.241018          3     671 13.752851  17563


###### PRE jgil
#    chr in_or_del  meanSize medianSize maxSize numRow
# 1   2L       DEL  7.280124          6     175  49196
# 2   2L        IN  6.588971          3    5454  38153
# 3   2R       DEL  7.340574          6     373  39545
# 4   2R        IN  6.794129          4    4698  32671
# 5   3L       DEL  7.371663          6     156  49822
# 6   3L        IN  7.088439          3    5977  38637
# 7   3R       DEL  7.440418          6     370  56024
# 8   3R        IN  7.244478          4    8634  41337
# 9    4       DEL  6.426778          4      34    239
# 10   4        IN 52.207207          3    4598    111
# 11   X       DEL  7.224725          6     511  29040
# 12   X        IN  6.772103          4    2697  34419

chisq.test( 
    data.frame(
    numIndels = c(83353,66315,83434,89924,381,54188), 
    chrSize = c(23011544,21146708,24543557,27905053,1351857,22422827)
    )
)



ddply(subset(fdata_dyak, blast_stat == "OK"), .(chr,in_or_del), summarize, meanSize = mean(size), medianSize = median(size), maxSize = max(size) , numRow = length(size))
#    chr in_or_del   meanSize medianSize maxSize numRow
# 1   2L       DEL   6.952909        5.0    4587  30409
# 2   2L        IN  12.587687        3.0    8578  23715
# 3   2R       DEL   6.824140        5.0    2915  25725
# 4   2R        IN  14.164956        4.0    8555  20854
# 5   3L       DEL   6.648035        5.0     188  29878
# 6   3L        IN  13.944812        3.0    6133  21780
# 7   3R       DEL   6.909028        5.0    2915  34571
# 8   3R        IN  14.468929        4.0    8213  24492
# 9    4       DEL   5.614286        3.5      34     70
# 10   4        IN 274.500000        5.0    5390     44
# 11   X       DEL   7.233582        5.0    2915  17450
# 12   X        IN  13.293439        4.0    7549  17237








##################################
############# OLD!


#clr <- ifelse(subset(fdata, abs(size) < 100 )$size > 0, "in", "del")
less100 <- subset(fdata, abs(size) < 100 )
less100$clr <- clr
ggplot(less100, aes(x=size, fill=clr)) + geom_histogram(binwidth=1) + facet_grid(chr ~ .) + theme(legend.position="bottom")

ggplot(data, aes(x=HET)) + geom_histogram(binwidth=1) + theme(legend.position="bottom") + facet_grid(chr ~ .)
ggplot(fdata, aes(x=HET)) + geom_histogram(binwidth=1) + theme(legend.position="bottom") + facet_grid(chr ~ .)

ggplot(data, aes(x=HETQ)) + geom_histogram(binwidth=1) +  theme(legend.position="bottom") + facet_grid(chr ~ .)
ggplot(fdata, aes(x=HETQ)) + geom_histogram(binwidth=1) +  theme(legend.position="bottom") + facet_grid(chr ~ .)

ggplot(data, aes(x=MISS)) + geom_histogram(binwidth=1) +  theme(legend.position="bottom") + facet_grid(chr ~ .)
ggplot(fdata, aes(x=MISS)) + geom_histogram(binwidth=1) +  theme(legend.position="bottom") + facet_grid(chr ~ .)




# table(cut(data$size, c(-Inf, 0, Inf)))
prop.table(table(cut(fdata$size, c(-Inf, 0, Inf))))*100
# (-Inf,0] (0, Inf] 
#  483997   332488 
# 59.27813 40.72187 

#######
x <- as.data.frame(with(data, tapply(AC, size, mean)))
plot( x$AC , xlim=c(-100,100), ylim=(0,80000) )
#######




summary(fdata$size)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -26410.00    -10.00     -3.00    -13.04      1.00  96900.00 

summary(subset(fdata, size < 0)$size)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -26410.00    -14.00     -6.00    -29.58     -2.00     -1.00 

summary(subset(fdata, size > 0)$size)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1.00     2.00     4.00    24.63     8.00 96900.00 



##### count intervals
prop.table(table(cut(abs(fdata$size), c(-Inf,1,2,5,10,20,50,100,1000,10000,Inf))))*100 ## table counts, cut create bins, prop.table count frequency
# (-Inf,1]         (1,2]         (2,5]        (5,10]        (10,20]         (20,50]     (50,100]     (100,1e+03]  (1e+03,1e+04]   (1e+04, Inf]
# 19.054096081  10.555402637  22.367136979  19.717204720  15.265261329   8.486003224   1.721937636   2.508016286   0.319936513     0.005004594 
# 106605         59056        125141        110315          85407           47478         9634         14032          1790             28  

prop.table(table(cut(subset(fdata, size > 0)$size, c(-Inf,1,2,5,10,20,50,100,1000,10000,Inf))))*100 ## table counts, cut create bins, prop.table count frequency
# (-Inf,1]         (1,2]         (2,5]        (5,10]       (10,20]       (20,50]        (50,100]     (100,1e+03]   (1e+03,1e+04]   (1e+04, Inf]
# 24.192726485  13.908591979  26.441905007  17.785496313  11.306807865   4.011035032   0.622624159   1.366493097   0.355534209      0.008785854
# 41304         23746         45144         30365         19304          6848          1063          2333           607              15

prop.table(table(cut(subset(fdata, size < 0)$size, c(-Inf,-10001, -1001, -101, -51, -21, -11, -6, -3, -2, Inf))))*100 ## table counts, cut create bins, prop.table count frequency
# (-Inf,-1e+04] (-1e+04,-1e+03]   (-1e+03,-101]   (-101,-51]       (-51,-21]       (-21,-11]        (-11,-6]         (-6,-3]      (-3, -2]     (-2, Inf] 
# 0.003343991     0.304303202     3.009334880     2.204719143    10.451258755    17.003680963    20.565546087    20.577635901    9.082794651   16.797382427
# 13               1183           11699            8571           40630           66103           79950          79997              35310          65301




### Counts size by Chr
ddply(fdata, .(chr), summarize, meanSize = mean(abs(size)), medianSize = median(abs(size)), maxSize = max(abs(size)) , numRow = length(abs(size)))
ddply(subset(fdata, size > 0), .(chr), summarize, meanSize = mean(size), medianSize = median(size), maxSize = max(size) , numRow = length(size))
ddply(subset(fdata, size < 0), .(chr), summarize, meanSize = mean(size), medianSize = median(size), maxSize = min(size) , numRow = length(size))

## Counts MISS by chr
ddply(data, .(chr), summarize, meanMiss = mean(MISS), medianMiss = median(MISS), maxMiss = max(MISS) , numRow = length(MISS))
ddply(fdata, .(chr), summarize, meanMiss = mean(MISS), medianMiss = median(MISS), maxMiss = max(MISS) , numRow = length(MISS))

## Counts HET by chr
ddply(data, .(chr), summarize, meanh = mean(HET), medianh = median(HET), maxh = max(HET) , numRow = length(HET))
ddply(fdata, .(chr), summarize, meanh = mean(HET), medianh = median(HET), maxh = max(HET) , numRow = length(HET))

## Counts HETQ by chr
ddply(data, .(chr), summarize, meanhq = mean(HETQ), medianhq = median(HETQ), maxhq = max(HETQ) , numRow = length(HETQ))
ddply(fdata, .(chr), summarize, meanhq = mean(HETQ), medianhq = median(HETQ), maxhq = max(HETQ) , numRow = length(HETQ))






########INFO RSTUDIO SERVER update R in debian
# --> http://cran.r-project.org/bin/linux/debian/
