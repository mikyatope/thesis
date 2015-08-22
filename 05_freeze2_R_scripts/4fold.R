
##### 4fold

all4fold  <- read.delim("~/SV_analysis/4fold_annotations/4fold_positions.csv", header=F, quote="")
names(all4fold) <- c("chr", "start")

all4fold <- merge(all4fold, subset(fSNPdata, select=c("chr","start","altc","refc", "class2")), all.x=TRUE, all.y=FALSE)

polymorphic4fold <- subset(all4fold, class2 == "CDS")
fixed4fold <- subset(all4fold, is.na(class2) )



##### FixNum

adply(fdataCl_graphs[sample(nrow(fdataCl_graphs), 15),c("chr","in_or_del","altc","refc")], 1, transform, 
      fixnum = ifelse(sum(altc,refc) > 175, min(table(sample(c(rep("REF", refc), rep("ALT",altc)), 175))), NA)  
)


fixnum_fdataCl_graphs <- fdataCl_graphs