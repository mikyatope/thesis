
### CALC SFS for DFE-alpha

#################
# SFS for INDELS
fixNum_INDELs <- subset(new.fdataCl_dsim, select=c(chr,start,refc,altc))
write.table(fixNum_INDELs, "fixedNum/INDELs.csv", sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE) 

########################## -> do fix nums outside with Python
# > OUTSIDE WORLD <
####### <- recover fix nums

fixNum_INDELs <- read.delim("fixedNum/INDELs.fixed149.csv", header=F, quote="")
names(fixNum_INDELs) <- c("chr", "start", "refc", "altc", "fixNum")

# recover annotations
fixNum_INDELs <- merge(fixNum_INDELs, subset(new.fdataCl_dsim, select=c(chr, start, class, class2, class3, classSize, size, derived, in_or_del, aut_or_X, newPI)), all.y=FALSE)

# remove NAs in frequencies or derived
fixNum_INDELs <- subset(fixNum_INDELs, refc != "NA" | altc != "NA" )
fixNum_INDELs <- subset(fixNum_INDELs, derived != "NA")

# MAF
fixNum_INDELs$SFS <- ifelse(fixNum_INDELs$refc <= fixNum_INDELs$altc, fixNum_INDELs$refc, fixNum_INDELs$altc)


# SFS
# intervals from 1 to 148, 149 does not exist since is fixation!!

# # UNFOLDED
# cat(
#     as.data.frame (
#         table( 
#             cut( subset(fixNum_INDELs, chr == "X" & class == "CDS" & in_or_del == "IN")$SFS, c(0:148) ) 
#         ) 
#     )$Freq
# , sep=",")


# FOLDED class3
write.table(
            
            ddply(subset(fixNum_INDELs, in_or_del != "NA"), .(class3, aut_or_X, in_or_del), function(i)  table( cut(  (i$SFS), c(1:148) ) )  )
            
            , "selection/DFEalpha_inputs/SFS.folded.indels.byclass3.txt",sep=" ", col.names=FALSE, quote=FALSE, row.names=FALSE) 

#### SUM by class
write.table(
    
    ddply(subset(fixNum_INDELs, in_or_del != "NA"), .(class2, aut_or_X, in_or_del), function(i)  nrow(i)  )
    
    , "selection/DFEalpha_inputs/SFS.TOTAL_COUNT.indels.byclass2.txt",sep=" ", col.names=FALSE, quote=FALSE, row.names=FALSE) 


# FOLDED class2
write.table(
    
    ddply(subset(fixNum_INDELs, in_or_del != "NA"), .(class2, aut_or_X, in_or_del), function(i)  table( cut(  (i$SFS), c(1:148) ) )  )
    
    , "selection/DFEalpha_inputs/SFS.folded.indels.byclass2.txt",sep=" ", col.names=FALSE, quote=FALSE, row.names=FALSE) 



# PI by class

write.table(
    
    ddply(subset(fixNum_INDELs, in_or_del != "NA"), .(class2, aut_or_X, in_or_del), function(i)  sum(i$newPI)  )
    
    , "selection/DFEalpha_inputs/SFS.sum_PI.indels.byclass.txt",sep=" ", col.names=FALSE, quote=FALSE, row.names=FALSE) 



################
# divergent INDELS counts

write.table(
    
    ddply(subset(fdivergenceCl, in_or_del != "NA" & sp=="polDmel"), .(class2, aut_or_X, in_or_del), function(i) nrow(i)  )
    
    , "selection/DFEalpha_inputs/DIVERGENCE_counts.INDELS.byclass2.txt",sep=" ", col.names=FALSE, quote=FALSE, row.names=FALSE) 



write.table(
    
    ddply(subset(fdivergenceCl, in_or_del != "NA"), .(class2, aut_or_X, in_or_del), function(i) nrow(i)  )
    
    , "selection/DFEalpha_inputs/DIVERGENCE_counts.INDELS.byclass2.txt",sep=" ", col.names=FALSE, quote=FALSE, row.names=FALSE) 




################
# Classes total nucleotide counts

dmelAnnotations$aut_or_X <- ifelse(dmelAnnotations$V1 != "X", "autosome", "X")

# V3 is class column
dmelAnnotations$V3 <- as.character(dmelAnnotations$V3) # make sure is not factor!
dmelAnnotations$class <- dmelAnnotations$V3 
dmelAnnotations$class2[dmelAnnotations$class != "intron"] <- dmelAnnotations$class[dmelAnnotations$class != "intron"] # the same if annotation is not 'intron'
dmelAnnotations$class2 <- as.character(dmelAnnotations$class2) # make sure is not factor!
dmelAnnotations$class2[dmelAnnotations$class == "intron"] <- ifelse(dmelAnnotations$size[dmelAnnotations$class == "intron"] <= 100, "small_intron", "long_intron")


write.table(
    
    ddply(dmelAnnotations, .(class2, aut_or_X), function(i) sum(i$size)  )
    
    , "selection/DFEalpha_inputs/NUCLEOTIDE_counts.byclass.txt",sep=" ", col.names=FALSE, quote=FALSE, row.names=FALSE) 










#############################################################
#############################################################
## SNPS SFS, not used finally

# SFS SNPs!
# intervals from 1 to 148, 149 does not exist since is fixation!!

fixNum_SNPs <- read.delim("fixedNum/SNPs.fixed149.csv", header=F, quote="")
names(fixNum_SNPs) <- c("chr", "start", "refc", "altc", "fixNum")

# recover annotations for SNPs
fixNum_SNPs <- merge(fixNum_SNPs, subset(fSNPdata, select=c(chr,start, class, class2, classSize, derived)), all.y=FALSE)

#annotate Small/Long introns
fixNum_SNPs$class3 <- NA
fixNum_SNPs$class3[fixNum_SNPs$class2 != "intron"] <- fixNum_SNPs$class2[fixNum_SNPs$class2 != "intron"] # the same if annotation is not 'intron'
fixNum_SNPs$class3[fixNum_SNPs$class2 == "intron"] <- ifelse(fixNum_SNPs$classSize[fixNum_SNPs$class2 == "intron"] <= 100, "small_intron", "long_intron")


# Autosome or X
fixNum_SNPs$aut_or_X <- NA
fixNum_SNPs$aut_or_X[fixNum_SNPs$chr != "X"] <- "autosome"
fixNum_SNPs$aut_or_X[fixNum_SNPs$chr == "X"] <- "X"

# remove NAs in frequencies or derived
fixNum_SNPs <- subset(fixNum_SNPs, refc != "NA" | altc != "NA" )
fixNum_SNPs <- subset(fixNum_SNPs, derived != "NA")

# MAF
fixNum_SNPs$SFS <- ifelse(fixNum_SNPs$refc <= fixNum_SNPs$altc, fixNum_SNPs$refc, fixNum_SNPs$altc)



# # UNFOLDED
# cat(
#     as.data.frame (
#         table( 
#             cut( subset(fixNum_SNPs, chr != "X" & class3 == "small_intron" )$SFS, c(0:148) ) 
#         ) 
#     )$Freq
#     , sep=",")

# FOLDED

# FOLDED
write.table(
    
    ddply(subset(fixNum_SNPs, class3 == "small_intron"), .(aut_or_X), function(i)  table( cut(  (i$SFS), c(1:148) ) )  )
    
    , "selection/DFEalpha_inputs/SFS.folded.SNPs.smallIntrons.txt",sep=" ", col.names=FALSE, quote=FALSE, row.names=FALSE) 










#################################
#################################
#### DFE alpha results

DFEalpha.results <- read.delim("selection/DFEalpha_inputs/DFEalpha.results.csv", header=T, quote="")


####  HEATMAP 

DFEalpha.results$chrindel <- interaction(DFEalpha.results$chr, DFEalpha.results$indel) #mix categories into single column

ggplot(DFEalpha.results, aes(x=chrindel, y=class)) + 
    geom_raster(aes(fill=omega.small_introns, width=.1, height=.1)) + 
    geom_text(aes(label=omega.small_introns)) +
    scale_fill_gradient(low="navajowhite", high="orange3")





