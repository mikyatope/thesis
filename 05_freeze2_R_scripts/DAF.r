########## Derived Allele Frequency

# Get rellevant data from other sets, infer Derived frequency

### INDELS
DAF_indels <- subset(fdataCl_graphs, in_or_del != "NA", select=c("chr", "start", "refc", "altc", "derived", "in_or_del", "class3", "aut_or_X"))
# Select derived freq. in new column
DAF_indels$DAF <- ifelse(DAF_indels$derived=="REF", DAF_indels$refc, DAF_indels$altc)
# # Annotate frameshift/non-frameshift in new column
# DAF_indels$class3 <- DAF_indels$class2
# DAF_indels$class3[DAF_indels$size %% 3 == 0 & DAF_indels$class2 == "CDS"] <- "CDS_non_frameshift"
# DAF_indels$class3[DAF_indels$size %% 3 != 0 & DAF_indels$class2 == "CDS"] <- "CDS_frameshift"
# # Annotate X/Autosome in new column
# DAF_indels$aut_or_X <- ifelse(DAF_indels$chr=="X", "X", "autosome")
row.names(DAF_indels)<-NULL #reset row.names

# SNPs
DAF_snps <- subset(fSNPdata, (class=="SYN" | class=="NOSYN") & derived!="NA" & class2=="CDS", select=c("chr", "refc", "altc", "derived", "class", "class2"))
# Add small introns
DAF_snps <- rbind(DAF_snps, subset(fSNPdata, class2=="intron" & classSize <= 100, select=c("chr", "refc", "altc", "derived", "class", "class2")))
# Select derived freq. in new column
DAF_snps$DAF <- ifelse(DAF_snps$derived=="REF", DAF_snps$refc, DAF_snps$altc)
# Annotate X/Autosome in new column
DAF_snps$aut_or_X <- ifelse(DAF_snps$chr=="X", "X", "autosome")
# create class3 column. Make sure SYN/NOSYN freeze2 annotations are from our CDS annotations. Anything else must be small_intron
DAF_snps$class3 <- "SNPs_small_intron"                                     
DAF_snps$class3[DAF_snps$class=="SYN" & DAF_snps$class2=="CDS"] <- "SNPs_synonymous"   
DAF_snps$class3[DAF_snps$class=="NOSYN" & DAF_snps$class2=="CDS"] <- "SNPs_non_synonymous" 
row.names(DAF_snps)<-NULL   #reset row.names


# Create dataframes with counts
DAF_indels_graph <- ddply(DAF_indels, .(class3, in_or_del, aut_or_X), function(i) prop.table( table( cut(  ((i$DAF*100)/205), c(0,10,20,30,40,50,60,70,80,90,100) ) ) ) )
DAF_snps_graph <- ddply(DAF_snps, .(class3, aut_or_X), function(i) prop.table( table( cut(  ((i$DAF*100)/205), c(0,10,20,30,40,50,60,70,80,90,100) ) ) ) )

# ddply(DAF_indels, .(class3, in_or_del, aut_or_X, cut(  ((DAF*100)/205), c(0,10,20,30,40,50,60,70,80,90,100) )), function(i) nrow(i) ) 

# merge data.frames filling empty columns
DAF_graph <- rbind.fill( DAF_snps_graph, DAF_indels_graph )


##### Counts
DAF_indels_counts<- ddply(DAF_indels, .(class3, in_or_del, aut_or_X), function(i)  nrow(i)  )
DAF_snps_counts <- ddply(DAF_snps, .(class3, aut_or_X), function(i)  nrow(i)  )

# merge data.frames filling empty columns
DAF_counts <- rbind.fill( DAF_snps_counts, DAF_indels_counts )



## Temporary duplication of SNPs values to add to both Insertions & Deletions graph
distribute.na.in_or_del <- function(dat) {
    rbind(
        transform(subset(dat, in_or_del %in% c("IN", NA)), in_or_del="IN"),
        transform(subset(dat, in_or_del %in% c("DEL", NA)), in_or_del="DEL")
    )
}

ggplot(distribute.na.in_or_del(melt(DAF_graph)), aes(x=variable, y=value*100,  fill=class3))+
    geom_bar(stat="identity", position="dodge", colour="black")+
    #geom_errorbar(ymax=(value*100)+0.1, ymin=(value*100)-0.1)+
    facet_grid(in_or_del ~ aut_or_X) +
    xlab('Derived Allele Frequency (%)')+
    ylab("Allele count (%)")+
    guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, label.theme=element_text(size=8, angle=0))) +
    theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) # xlabels vertical, centered & aligned
    