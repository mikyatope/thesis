


##### r2
pre_r2_dels_data <- getAndFilterFreeze2r2("~/SV_analysis/r2_dels.ld")
r2_dels_data<- subset(pre_r2_dels_data, r2 != "NaN"); rm(pre_r2_dels_data)

pre_r2_ins_data <- getAndFilterFreeze2r2("~/SV_analysis/r2_ins.ld")
r2_ins_data<- subset(pre_r2_ins_data, r2 != "NaN"); rm(pre_r2_ins_data)

pre_r2_total_data <- getAndFilterFreeze2r2("~/SV_analysis/r2_total.ld")
r2_total_data<- subset(pre_r2_total_data, r2 != "NaN"); rm(pre_r2_total_data)


#slideData<- subset(r2_total_data, chr == "2L" & pos1 >= 1 & pos2 <= 100000);rm(slideData)
#mean(combn(slideData$pos1,2, diff))

system.time( 
    R2_TOTAL  <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), sumAndMeansR2, dataset=r2_total_data) )) 
)
system.time(     
    R2_INS  <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), sumAndMeansR2, dataset=r2_ins_data) )) 
)
system.time( 
    R2_DELS  <- data.frame(do.call("rbind",  by(SWs, 1:nrow(SWs), sumAndMeansR2, dataset=r2_dels_data) )) 
)


ggplot(R2_TOTAL [c("chr", "start","r2")], aes(start, r2, colour=chr)) + 
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +  #coord_cartesian(ylim = c(0, 0.01)) +
    theme(legend.position = "none",axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks=element_blank()) +
    labs(title = "r² TOTAL") 

ggplot(R2_INS [c("chr", "start","r2")], aes(start, r2, colour=chr)) + 
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +  #coord_cartesian(ylim = c(0, 0.01)) +
    theme(legend.position = "none",axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks=element_blank()) +
    labs(title = "r² INSERTIONS") 

ggplot(R2_DELS [c("chr", "start","r2")], aes(start, r2, colour=chr)) + 
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +  #coord_cartesian(ylim = c(0, 0.01)) +
    theme(legend.position = "none",axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks=element_blank()) +
    labs(title = "r² DELETIONS") 



#CORRELATIONS
#PI

#Original (contains freeze1 pi and comeron rec!)
#correls_PI_adim  <- read.delim("pis_recombinations_original.csv", header=T, quote="")
correls_PI_adim  <- read.delim("pis_recombinations.csv", header=T, quote="")


#correls_PI_adim          <- PI[c('chr', 'start', 'PIindel_adim')]
#### correls_PI_adim$PI_indel    <- PI$PIindel_adim

correls_PI_adim$PI_in    <- PI_IN$PIindel_adim
correls_PI_adim$PI_del   <- PI_DEL$PIindel_adim
correls_PI_adim$PI_snp   <- PI_SNP$PI_SNP_adim
correls_PI_adim$r2_total <- R2_TOTAL$r2
correls_PI_adim$r2_in    <- R2_INS$r2
correls_PI_adim$r2_del   <- R2_DELS$r2

correls_PI_adim$DIVdmel   <- DIV$kdmel
correls_PI_adim$pi_div   <- (correls_PI_adim$PI_indel/DIV$kdmel)

#colnames(correls_PI_adim)[3] <- "PI_indel"

#write.table(correls_PI_adim, file="pis_recombinations.csv", sep="\t", row.names=F, quote=F)





####function
loopCorrels <- function(chrom, method, var1, var2, lab1, lab2){
    
    cat(chrom, " MAIN ", method, lab1, " + ", lab2)
    print(cor.test(~ var1 + var2, method = method))
    
    #     data <- cbind(var1, var2)
    #     #print(data)
    #     rowSlice <- nrow(data)
    #     sliceDF1 <- data[c(1:(rowSlice/3)),]
    #     sliceDF2 <- data[((rowSlice/3)+1):((rowSlice/3)*2),]
    #     sliceDF3 <- data[(((rowSlice/3)*2)+1):rowSlice,]
    #     
    #     cat(chrom, " 0 SEGMENT ", method, lab1, " + ", lab2)
    #     print(cor.test(~ var1 + var2, data = sliceDF1, method = method))
    #     cat(chrom, " 1 SEGMENT ", method, lab1, " + ", lab2)
    #     print(cor.test(~ var1 + var2, data = sliceDF2, method = method))
    #     cat(chrom, " 2 SEGMENT ", method, lab1, " + ", lab2)
    #     print(cor.test(~ var1 + var2, data = sliceDF3, method = method))
}


sink("correlations.out")
#####
for(method in c("spearman")){ #, "kendall", "pearson"
    
    
    loopCorrels (chrom, method, correls_PI_adim$DIVdmel, correls_PI_adim$rec_comeron, "DIVdmel", "rec_comeron") }
    for (chrom in c("2L","2R","3L","3R","X")){
        sliceDF <- subset(correls_PI_adim, chr==chrom)
        loopCorrels (chrom, method, sliceDF$DIVdmel, sliceDF$rec_comeron, "DIVdmel", "rec_comeron")
    }
}    
    
#     for (chrom in c("2L","2R","3L","3R","4", "X")){
#         
#         sliceDF <- subset(correls_PI_adim, chr==chrom)
#         
#         if (chrom != "4"){
#             loopCorrels (chrom, method, sliceDF$Pi_freeze1, sliceDF$PI_snp, "Pi_freeze1", "PI_snp")
#             
#             loopCorrels (chrom, method, sliceDF$PI_indel, sliceDF$rec_comeron, "PI_indel", "rec_comeron")
#             loopCorrels (chrom, method, sliceDF$PI_in, sliceDF$rec_comeron, "PI_in", "rec_comeron")
#             loopCorrels (chrom, method, sliceDF$PI_del, sliceDF$rec_comeron, "PI_del", "rec_comeron")
#             loopCorrels (chrom, method, sliceDF$rec_comeron, sliceDF$PI_snp, "PI_snp", "rec_comeron")
#             loopCorrels (chrom, method, sliceDF$rec_comeron, sliceDF$pi_div, "Pi / divergence", "rec_comeron")
#         }
#         
#         loopCorrels (chrom, method, sliceDF$PI_indel, sliceDF$pi_div, "Pi / divergence", "PI_indel")
#         loopCorrels (chrom, method, sliceDF$PI_in, sliceDF$pi_div, "Pi / divergence", "PI_in")
#         loopCorrels (chrom, method, sliceDF$PI_del, sliceDF$pi_div, "Pi / divergence", "PI_del")
#         loopCorrels (chrom, method, sliceDF$PI_snp, sliceDF$pi_div, "Pi / divergence", "PI_snp")
#         
#         loopCorrels (chrom, method, sliceDF$PI_indel, sliceDF$PI_snp, "PI_indel", "PI_snp")
#         loopCorrels (chrom, method, sliceDF$PI_in, sliceDF$PI_snp, "PI_in", "PI_snp")    
#         loopCorrels (chrom, method, sliceDF$PI_del, sliceDF$PI_snp, "PI_del", "PI_snp")
#         loopCorrels (chrom, method, sliceDF$PI_in, sliceDF$PI_del, "PI_in", "PI_del") 
#         
#         loopCorrels (chrom, method, sliceDF$PI_indel, sliceDF$DIVdmel, "PI_indel", "DIVdmel")
#         loopCorrels (chrom, method, sliceDF$PI_in, sliceDF$DIVdmel, "PI_in", "DIVdmel")    
#         loopCorrels (chrom, method, sliceDF$PI_del, sliceDF$DIVdmel, "PI_del", "DIVdmel")
#         loopCorrels (chrom, method, sliceDF$PI_snp, sliceDF$DIVdmel, "PI_snp", "DIVdmel")
#         
#         
#         
#         
#         #loopCorrels (chrom, method, sliceDF$r2_total, sliceDF$PI_snp, "r2_total", "PI_snp")
#         #         loopCorrels (chrom, method, sliceDF$r2_in, sliceDF$PI_snp, "r2_in", "PI_snp")    
#         #         loopCorrels (chrom, method, sliceDF$r2_del, sliceDF$PI_snp, "r2_del", "PI_snp")
#         
#         # loopCorrels (chrom, method, sliceDF$r2_total, sliceDF$PI_indel, "r2_total", "PI_indel")
#         #         loopCorrels (chrom, method, sliceDF$r2_in, sliceDF$PI_indel, "r2_in", "PI_indel")    
#         #         loopCorrels (chrom, method, sliceDF$r2_del, sliceDF$PI_indel, "r2_del", "PI_indel")
#         
#         #         loopCorrels (chrom, method, sliceDF$r2_total, sliceDF$PI_in, "r2_total", "PI_in")
#         #         loopCorrels (chrom, method, sliceDF$r2_in, sliceDF$PI_in, "r2_in", "PI_in")    
#         #         loopCorrels (chrom, method, sliceDF$r2_del, sliceDF$PI_in, "r2_del", "PI_in")
#         #         
#         #         loopCorrels (chrom, method, sliceDF$r2_total, sliceDF$PI_del, "r2_total", "PI_del")
#         #         loopCorrels (chrom, method, sliceDF$r2_in, sliceDF$PI_del, "r2_in", "PI_del")    
#         #         loopCorrels (chrom, method, sliceDF$r2_del, sliceDF$PI_del, "r2_del", "PI_del")
#         
#         
#         rm(sliceDF)
#     }
    
}
sink()



##### Extra Graphs

ggplot(correls_PI_adim, aes(start, rec_comeron, colour=chr)) + 
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +  #coord_cartesian(ylim = c(0, 0.01)) +
    theme(legend.position = "none",axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks=element_blank()) +
    labs(title = "Recombination (Comeron)") 

ggplot(correls_PI_adim, aes(start, Pi_freeze1, colour=chr)) + 
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +  #coord_cartesian(ylim = c(0, 0.01)) +
    theme(legend.position = "none",axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks=element_blank()) +
    labs(title = "Pi SNPs Freeze1") 