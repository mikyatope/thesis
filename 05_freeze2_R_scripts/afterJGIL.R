
####   AFTER JGIL!
########################
new.fdata <- rbind(fdata, fdata_MS, fdata_TR)
new.fdata <- new.fdata[with(new.fdata, order(chr, start)), ]
row.names(new.fdata)<-NULL #reset row.names

#### Merge with JGIL 
# new.fdata <- merge(new.fdata, jgil, all.y=FALSE)
# new.fdata <- new.fdata[with(new.fdata, order(chr, start)), ]
# row.names(new.fdata)<-NULL #reset row.names

new.fdata$size <- ifelse(new.fdata$ALN > new.fdata$RLN, new.fdata$ALN-1, new.fdata$RLN-1)
#########################

ddply(subset(fdata_dsim, indelType=="MS"), .(), summarize, meanSize = mean(size), medianSize = median(size), maxSize = max(size) , sd = sd(size), numRow = length(size))


# output semi.vcf to do blasts and polarizations
outdf <- merge(jgil, data, by=(c("chr","start")), all.y=FALSE)[,c("chr","start","id", "REF", "ALT", "qual")]
write.table(outdf, "blasts/method2/jgil.semi.vcf", sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
rm(outdf)



#################
#########  INDEL SIZE DISTRIBUTION
less50 <- subset(new.fdata, size <= 50)
ggplot(less50, aes(x=size, fill=chr)) + geom_histogram(binwidth=1, origin=-0.5) + 
    facet_grid(indelType ~ chr, scales="free_y") + 
    scale_y_continuous(expand = c(0,0)) + # x axis starts at y zero
    scale_x_continuous(expand = c(0,0)) + # y axis starts at x zero
    labs(title = "INDEL size distribution per chromosome") + 
    #theme(legend.position="bottom", plot.title=element_text(vjust=-60))
rm(less50)




########## BY CLASSES

library("GenomicRanges")

# get dmel filtered annotations
dmelAnnotations  <- read.delim("~/SV_analysis/dmel_anotations/dmel_FINAL_w_intergenic.gff", header=F, quote="")
# create GenomicRanges object for annotations
dmelAnnotations.gr <- with(dmelAnnotations, GRanges(V1, IRanges(V4, V5), functionalClass=V3))

# create GenomicRanges object for INDELS
gr_test <- with(new.fdata, GRanges(chr, IRanges(start-1, start)))

foo <- findOverlaps(gr_test, dmelAnnotations.gr, type="within") 


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
new.fdataCl <- merge(as.data.frame(foo.DT.ord)[,c("rn", "queryHits", "classCode")],foo_merged)[,c("queryHits", "classCode","V3","V4","V5")]
names(new.fdataCl) <- c("queryHits", "classCode","class", "classStart","classEnd")
row.names(new.fdataCl)<-NULL #reset row.names

# get INDEL info names merging by queryHit (queryHit == annotation row.name)
new.fdataCl <- merge(new.fdata, new.fdataCl, by.x="row.names", by.y="queryHits", all.x=FALSE)

#sort
new.fdataCl <- new.fdataCl[with(new.fdataCl, order(chr, start)), ]
row.names(new.fdataCl)<-NULL #reset row.names
new.fdataCl <- subset(new.fdataCl, select=-c(Row.names))

rm(foo, foo_merged, foo.DT, foo.DT.ord)


#########  INDEL SIZE DISTRIBUTION BY CLASS
less50 <- subset(new.fdataCl, size <= 50)
ggplot(less50, aes(x=size, fill=class)) + geom_histogram(binwidth=1, origin=-0.5) + 
    facet_grid(class ~ indelType, scales="free_y") + 
    scale_y_continuous(expand = c(0,0)) + # x axis starts at y zero
    scale_x_continuous(expand = c(0,0)) + # y axis starts at x zero
    labs(title = "INDEL size distribution by functional class") + 
    #theme(legend.position="bottom", plot.title=element_text(vjust=-60))
rm(less50)



##### add classes

#annotate Small/Long introns

new.fdataCl$classSize <- (new.fdataCl$classEnd-new.fdataCl$classStart)+1
new.fdataCl$class <- as.character(new.fdataCl$class) # make sure is not factor!



new.fdataCl$class2[new.fdataCl$class != "intron"] <- new.fdataCl$class[new.fdataCl$class != "intron"]
new.fdataCl$class2[new.fdataCl$class == "intron"] <- ifelse(new.fdataCl$classSize[new.fdataCl$class == "intron"] <= 100, "small_intron", "long_intron")

### New PI
new.fdataCl$newPI <-  2*(new.fdataCl$FAC*new.fdataCl$RC)/(((new.fdataCl$FAC+new.fdataCl$RC)*(new.fdataCl$FAC+new.fdataCl$RC-1)))


### Annotate frameshift indels
new.fdataCl$class3 <- new.fdataCl$class2
new.fdataCl$class3[new.fdataCl$size %% 3 == 0 & new.fdataCl$class2 == "CDS"] <- "CDS_non_frameshift"
new.fdataCl$class3[new.fdataCl$size %% 3 != 0 & new.fdataCl$class2 == "CDS"] <- "CDS_frameshift"
# Annotate X/Autosome in new column
new.fdataCl$aut_or_X <- ifelse(new.fdataCl$chr=="X", "X", "autosome")
row.names(new.fdataCl)<-NULL #reset row.names




############## MAF

# Get rellevant data from other sets, infer Minor allele

### INDELS
# search for minor allele
new.fdataCl$MAF <- ifelse(new.fdataCl$refc < new.fdataCl$altc , new.fdataCl$refc, new.fdataCl$altc)
MAF_indels <- subset(new.fdataCl, select=c("chr", "start", "MAF", "refc", "altc", "class", "class2", "classSize", "size", "class3", "aut_or_X"))


# SNPs
MAF_snps <- subset(fSNPdata, (class=="SYN" | class=="NOSYN") & class2=="CDS", select=c("chr", "refc", "altc", "class", "class2"))
# Add small introns
MAF_snps <- rbind(MAF_snps, subset(fSNPdata, class2=="intron" & classSize <= 100, select=c("chr", "refc", "altc", "class", "class2")))
# Select derived freq. in new column
MAF_snps$MAF <- ifelse(MAF_snps$refc < MAF_snps$altc, MAF_snps$refc, MAF_snps$altc)
# Annotate X/Autosome in new column
MAF_snps$aut_or_X <- ifelse(MAF_snps$chr=="X", "X", "autosome")
# create class3 column. Make sure SYN/NOSYN freeze2 annotations are from our CDS annotations. Anything else must be small_intron
MAF_snps$class3 <- "SNPs_small_intron"                                     
MAF_snps$class3[MAF_snps$class=="SYN" & MAF_snps$class2=="CDS"] <- "SNPs_synonymous"   
MAF_snps$class3[MAF_snps$class=="NOSYN" & MAF_snps$class2=="CDS"] <- "SNPs_non_synonymous" 
row.names(MAF_snps)<-NULL   #reset row.names


# Create dataframes with counts
MAF_indels_graph <- ddply(MAF_indels, .(class3,  aut_or_X), function(i) prop.table( table( cut(  ((i$MAF*100)/205), c(0,5,10,15,20,25,30,35,40,45,50) ) ) ) )
MAF_snps_graph <- ddply(MAF_snps, .(class3, aut_or_X), function(i) prop.table( table( cut(  ((i$MAF*100)/205), c(0,5,10,15,20,25,30,35,40,45,50) ) ) ) )


# merge data.frames filling empty columns
MAF_graph <- rbind.fill( MAF_snps_graph, MAF_indels_graph )


##### Counts
MAF_indels_counts<- ddply(MAF_indels, .(class3, aut_or_X), function(i)  table( cut(  ((i$MAF*100)/205), c(0,5,10,15,20,25,30,35,40,45,50) ) )  )
MAF_snps_counts <- ddply(MAF_snps, .(class3, aut_or_X), function(i)  table( cut(  ((i$MAF*100)/205), c(0,5,10,15,20,25,30,35,40,45,50) ) )  )

# merge data.frames filling empty columns
MAF_counts <- rbind.fill( MAF_snps_counts, MAF_indels_counts )



# reorder levels
# Convert class3 into factor
melted_MAF_graph <- melt(MAF_graph)
melted_MAF_graph <- transform(melted_MAF_graph, class3 = as.factor(class3))
# reorder levels using factor()
melted_MAF_graph$class3 <- factor(melted_MAF_graph$class3, levels=c("SNPs_non_synonymous",
                                                                    "SNPs_small_intron",
                                                                    "SNPs_synonymous",
                                                                    "CDS_frameshift",
                                                                    "CDS_non_frameshift",
                                                                    "five_prime_UTR",
                                                                    "three_prime_UTR",
                                                                    "small_intron",
                                                                    "long_intron",
                                                                    "intergenic" ))


ggplot(melted_MAF_graph, aes(x=variable, y=value*100,  fill=class3))+
    geom_bar(stat="identity", position="dodge", colour="black")+
    facet_grid(. ~ aut_or_X) +
    theme_bw() +
    labs(title="Minor Allele Frequency spectrum", x="Minor Allele Frequency (%)", y="Allele count (%)") +
    guides(fill = guide_legend(keywidth = 0.8, keyheight = 0.8, label.theme=element_text(size=14, angle=0))) +
    scale_x_discrete(labels=c("≤5%","≤10%","≤15%","≤20%","≤25%","≤30%","≤35%","≤40%","≤45%","≤50%"))+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.title = element_text(size=16),
          strip.text = element_text(size=16),
          axis.text.x = element_text(hjust = 0.5, vjust=0.5)) + # xlabels vertical, centered & aligned
    scale_fill_discrete(labels = c("non-synonymous SNPs ",
                                       "small intron SNPs",
                                       "synonymous SNPs",
                                       "CDS frameshift",
                                       "CDS non-frameshift",
                                       "5' prime UTR",
                                       "3'prime UTR",
                                       "small intron",
                                       "long intron",
                                       "intergenic" ))
rm(melted_MAF_graph)




########## POLARIZATION AND DERIVED ALLELES!


new.polarizations <- read.delim("blasts/method2/results/merged_results.csv", header=F, quote="", sep="") 
names(new.polarizations) <- c("chr", "start", "derived", "REF_num_hits", "REF_num_hsps", "REF_evalue", "REF_qstart", "REF_qend", "ALT_num_hits", "ALT_num_hsps", "ALT_evalue", "ALT_qstart", "ALT_qend", "aln_slice", "REF", "ALT")

new.fdataCl_dsim <- merge(new.polarizations, new.fdataCl, by=(c("chr", "start")), all.y=FALSE)
new.fdataCl_dsim <- new.fdataCl_dsim[with(new.fdataCl_dsim, order(chr, start)), ]
row.names(new.fdataCl_dsim)<-NULL #reset row.names

new.fdataCl_dsim$polarizationStatus <-ifelse((new.fdataCl_dsim$REF_qstart <=50 & new.fdataCl_dsim$REF_qend >= (150+new.fdataCl_dsim$RLN)) &
                                             (new.fdataCl_dsim$ALT_qstart <=50 & new.fdataCl_dsim$ALT_qend >= (150+new.fdataCl_dsim$ALN)) &
                                             (new.fdataCl_dsim$REF_num_hits == 1 & new.fdataCl_dsim$ALT_num_hits == 1) &
                                             (new.fdataCl_dsim$REF_num_hits >= 0.0000000001 & new.fdataCl_dsim$ALT_num_hits >= 0.0000000001) &
                                             (new.fdataCl_dsim$derived != "NA")       
                                              , "OK", "NA" ) 


new.fdataCl_dsim$polStat <-ifelse((new.fdataCl_dsim$REF_qstart <=50 & new.fdataCl_dsim$REF_qend >= (150+new.fdataCl_dsim$RLN)) &
                                                 (new.fdataCl_dsim$ALT_qstart <=50 & new.fdataCl_dsim$ALT_qend >= (150+new.fdataCl_dsim$ALN)) &
                                                 (new.fdataCl_dsim$REF_num_hits == 1 & new.fdataCl_dsim$ALT_num_hits == 1) &
                                                 (new.fdataCl_dsim$REF_num_hits >= 0.0000000001 & new.fdataCl_dsim$ALT_num_hits >= 0.0000000001) &
                                                 (new.fdataCl_dsim$derived != "NA")       
                                                , "OK", ifelse((new.fdataCl_dsim$size > 25) &
                                                               (new.fdataCl_dsim$REF_num_hits == 1 | new.fdataCl_dsim$ALT_num_hits == 1) &
                                                               (new.fdataCl_dsim$REF_num_hits >= 0.0000000001 | new.fdataCl_dsim$ALT_num_hits >= 0.0000000001) &
                                                               (new.fdataCl_dsim$derived != "NA")
                                                               , "longIndel", "NA"))
                                             


# check Insertion or Deletion
# (When REF/ALT is the derived, if smaller allele, always deletion)
new.fdataCl_dsim$in_or_del <- ifelse(new.fdataCl_dsim$polStat != "NA" & new.fdataCl_dsim$derived=="REF", 
                                     ifelse(new.fdataCl_dsim$RLN < new.fdataCl_dsim$ALN, "DEL", "IN"), 
                                        ifelse(new.fdataCl_dsim$polStat != "NA" & new.fdataCl_dsim$derived=="ALT", 
                                            ifelse(new.fdataCl_dsim$ALN < new.fdataCl_dsim$RLN, "DEL", "IN"), "NA"))

# Check dataset : 
# rstudio::viewData(new.fdataCl_dsim[,c("chr","start","polarizationStatus", "derived", "in_or_del", "RLN", "ALN", "size", "aln_slice", "REF", "ALT")])



prop.table(table(new.fdataCl_dsim$polarizationStatus))
# NA        OK 
# 117573    200777
# 0.3693199 0.6306801

table(subset(new.fdataCl_dsim, polarizationStatus=="OK")$in_or_del)
# DEL     IN 
# 136420  64357



#### Compare Reference IN/DEL with Polarized IN/DEL
new.fdataCl_dsim$in_or_del_byREF <- ifelse(new.fdataCl_dsim$RLN > new.fdataCl_dsim$ALN, "DEL", "IN") 
new.fdataCl_dsim$ref_pol_concordance <- ifelse(new.fdataCl_dsim$in_or_del  == "NA", "NA",
                                               ifelse(new.fdataCl_dsim$in_or_del == new.fdataCl_dsim$in_or_del_byREF, "SAME", "INVERSE"))

# rstudio::viewData(new.fdataCl_dsim[,c("chr","start","polarizationStatus", "derived", "in_or_del", "in_or_del_byREF", "ref_pol_concordance", "RLN", "ALN", "size", "aln_slice", "REF", "ALT")])

prop.table(table(subset(new.fdataCl_dsim, in_or_del  != "NA")$ref_pol_concordance))
#   INVERSE      SAME 
#     37166    163611 
# 0.1851108 0.8148892


#prop.table(table(subset(new.fdataCl_dsim, in_or_del  != "NA" & in_or_del == "DEL")$ref_pol_concordance))
#prop.table(table(subset(new.fdataCl_dsim, in_or_del  != "NA" & in_or_del == "IN")$ref_pol_concordance))



prop.table(table(subset(new.fdataCl_dsim, in_or_del  != "NA")$ref_pol_concordance, subset(new.fdataCl_dsim, in_or_del  != "NA")$derived))
#             ALT    REF
# INVERSE       0  37166
# SAME     163611      0

#prop.table(table(subset(new.fdataCl_dsim, polarizationStatus=="OK")$ref_pol_concordance, subset(new.fdataCl_dsim, polarizationStatus=="OK")$in_or_del))

#ddply(subset(new.fdataCl_dsim, polarizationStatus=="OK"), .(in_or_del,ref_pol_concordance), summarize, count=table(ref_pol_concordance),prop= prop.table(table(ref_pol_concordance)) )
# in_or_del ref_pol_concordance  count prop
# 1       DEL             INVERSE  20470    1
# 2       DEL                SAME 115950    1
# 3        IN             INVERSE  16696    1
# 4        IN                SAME  47661    1

#manually
#              INVERSE      SAME
# DELETIONS      20470    115950
# INSERTIONS     16696     47661

#              INVERSE      SAME
# DELETIONS  0.1500513 0.8499487
# INSERTIONS 0.2594279 0.7405721



ddply(subset(new.fdataCl_dsim, in_or_del !="NA" & size <= 150), .(in_or_del), summarize, meanSize = mean(size), medianSize = median(size), maxSize = max(size), sd = sd(size), numRow = length(size))





##### To graphs

fdataCl_graphs <- subset(new.fdataCl_dsim, size > 100)