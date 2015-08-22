

##### Figure1 
##### INDEL count per chromosome in 100kb sliding windows (dsim)
## needs PI_sizes from variation_estimates_sizes_distribution.r (depends on fdataCl_graphs)  !!!!

ggplot(PI_sizes[c("chr", "start",  "numINDELS", "numIN", "numDEL")], aes(x=start, colour=chr)) + 
    #geom_line()+
    #geom_point(aes(size=numINDELS)) +
    geom_line(aes(y=numIN), linetype="dashed") + 
    geom_line(aes(y=numDEL)) + 
    geom_line(size=1.2, aes(y=numINDELS), alpha=0.5) + 
    facet_grid(chr ~ ., scales="free_x") +
    #facet_wrap(~ chr, ncol=2, scales = "free_x") +
    #coord_cartesian(ylim = c(0, 0.001)) +
    theme_bw() +
    theme(legend.position = "none",
          #axis.text.x=element_blank(), 
          #axis.title.x=element_blank(), 
          #axis.ticks=element_blank(), 
          plot.margin=unit(rep(.5, 4), "lines")) +
    scale_x_continuous(breaks=c(0, 5000000,10000000,15000000,20000000,25000000), labels = comma) +
    labs(x="size (bp)", y="count") 
    #labs(title = "INDEL count per chromosome in 100kb sliding windows") +
    #+ scale_y_log10()
ggsave("plots/over100bp/PG3_count_x_chr_100kb.png", dpi=150)

################



##### Figure2 
##### INDEL average size per chromosome in 100kb sliding windows (dsim)
## needs PI_sizes from variation_estimates_sizes_distribution.r (depends on fdataCl_graphs) -> sizes <= 100 !!!!

ggplot(PI_sizes[c("chr", "start", "meanAbsSize", "meanInSize", "meanDelSize")], aes(x=start, colour=chr)) + 
    #geom_line()+
    #geom_point(aes(size=numINDELS)) +
    geom_line(aes(y=meanInSize), linetype="dashed") + 
    geom_line(aes(y=meanDelSize)) + 
    geom_line(size=1.2, aes(y=meanAbsSize), alpha=0.5) + 
    facet_grid(chr ~ ., scales="free_x") +
    #facet_wrap(~ chr, ncol=2, scales = "free_x") +
    #coord_cartesian(ylim = c(0, 0.001)) +
    theme_bw() +
    theme(legend.position = "none",
          #axis.text.x=element_blank(), 
          #axis.title.x=element_blank(), 
          #axis.ticks=element_blank(), 
          plot.margin=unit(rep(.5, 4), "lines"))+
    scale_x_continuous(breaks=c(0, 5000000,10000000,15000000,20000000,25000000), labels = comma) +
    labs(x="size (bp)", y="Mean size (bp)") 
    #labs(title = "INDEL average size per chromosome in 100kb sliding windows") +
    #+ scale_y_log10()
ggsave("plots/over100bp/PG4_averageSize_x_chr_100kb.png", dpi=150)

###################



##### Figure3 
##### INDEL size distribution per chromosome (dsim)
## needs new.fdataCl_dsim from afterJGIL.r

less50 <- subset(new.fdataCl_dsim, size > 100 & in_or_del != "NA" )
ggplot(less50, aes(x=size, fill=in_or_del)) + geom_histogram(binwidth=1, origin=-0.5) + 
    facet_grid(chr ~ in_or_del, scales="free_y") + 
    scale_y_continuous(expand = c(0,0)) + # x axis starts at y zero
    scale_x_continuous(expand = c(0,0)) + # y axis starts at x zero
    #labs(title = "INDEL size distribution per chromosome") + 
    theme_bw() + 
    theme(legend.position="none") +
    labs(x="size (bp)")
rm(less50)
ggsave("plots/over100bp/PG1_size_distribution_x_chr.png", dpi=150)

###################



##### Figure4 
##### INDEL size distribution by functional class (dsim)
## needs fdataCl_graphs from variation_estimates_class.r

ggplot(subset(fdataCl_graphs, size > 100 & in_or_del != "NA" ), aes(x=size, fill=in_or_del )) + 
    geom_histogram(binwidth=1, origin=-0.5) + # origin centers x label tick
    facet_grid(class2 ~ in_or_del , scales="free_y") +
    scale_y_continuous(expand = c(0,0)) + # x axis starts at y zero
    scale_x_continuous(expand = c(0,0)) + # y axis starts at x zero
    #labs(title = "INDEL size distribution by functional class") + 
    theme_bw() +
    theme(#panel.grid.major = element_line(colour = "grey", size = 0.7),
        #panel.grid.minor = element_line(colour = "grey"),
        strip.text.y = element_text(size = 8),
        legend.position="none") +
    labs(x="size (bp)")
ggsave("plots/over100bp/PG2_size_distribution_by_class.png", dpi=150)

###################



##### Figure5 
##### INDEL size by functional class (dsim)
## needs fdataCl_graphs from variation_estimates_class.r

ggplot(subset(fdataCl_graphs, in_or_del != "NA"), aes(class3, size, fill=in_or_del)) + #stat_boxplot(geom ='errorbar', width=0.2) +
    geom_boxplot(notch = TRUE, notchwidth = 0.8, outlier.shape = NA) +
    #coord_cartesian(ylim = ylim1*1.05) +
    coord_cartesian(ylim = c(0,25)) +
    facet_grid(. ~ in_or_del ) +
    theme_bw()+
    labs( y="size (bp)") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),                    # xlabels vertical, centered & aligned
          legend.position = "none",
          axis.title.x = element_blank()) +
    
ggsave("plots/size_by_class.png", dpi=150)

##################




##### Figure6 
##### INDEL diversity (PI indel) in 100kb sliding windows (dsim)
## needs PI, PI_IN, PI_DEL df's from variation_estimates_windows.r (depends on fdataCl_graphs)

SW_pi_plot <- ggplot(PI [c("chr", "start","PIindel_adim")], aes(start, PIindel_adim, colour=chr)) + 
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +
    coord_cartesian(ylim = c(0, 0.000015)) +
    theme_bw() +
    ylab("π indel") +
    theme(legend.position = "none",axis.text.x=element_blank(), 
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(), 
          plot.margin=unit(rep(.5, 4), "lines")) +
    scale_y_continuous(labels = comma) +
    labs(title = "π indel TOTAL") 

SW_pi_plot_DEL <- ggplot(PI_DEL [c("chr", "start","PIindel_adim")], aes(start, PIindel_adim, colour=chr)) + 
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +
    coord_cartesian(ylim = c(0, 0.000015)) +
    theme_bw() +
    ylab("π indel") +
    theme(legend.position = "none",axis.text.x=element_blank(), 
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(), 
          plot.margin=unit(rep(.5, 4), "lines")) + 
    scale_y_continuous(labels = comma) +
    labs(title = "π DELETIONS") 

SW_pi_plot_IN <- ggplot(PI_IN [c("chr", "start","PIindel_adim")], aes(start, PIindel_adim, colour=chr)) + 
    geom_point(size=2.2) + 
    facet_grid(. ~ chr, scales="free_x", space="free_x") +
    coord_cartesian(ylim = c(0, 0.000015)) +
    theme_bw() +
    ylab("π indel") +
    theme(legend.position = "none",
          axis.text.x=element_blank(), 
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(), 
          plot.margin=unit(rep(.5, 4), "lines")) +
    scale_y_continuous(labels = comma) +
    labs(title = "π INSERTIONS") 

grid.arrange(SW_pi_plot, SW_pi_plot_DEL, SW_pi_plot_IN)
ggsave("plots/PG_5_indel_diversity_chr_100kb.png", dpi=150)

########### 



##### Figure7 
##### INDEL diversity (PI indel) in 100kb sliding windows (dsim)
## needs PI, PI_IN, PI_DEL df's from variation_estimates_windows.r





##### Figure12
##### INDEL variation (Pi INDEL) by functional class (dsim)
## needs fdataCl_graphs from variation_estimates_class.r

ggplot(subset(fdataCl_graphs, in_or_del != "NA"), aes(class3, newPI, fill=in_or_del)) + 
    geom_boxplot(notch = TRUE, notchwidth = 0.8, outlier.shape = NA) +
    #coord_cartesian(ylim = c(0,25)) +
    facet_grid(. ~ in_or_del ) +
    theme_bw() +
    labs( y="π indel") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),                                 # xlabels vertical, centered & aligned
          legend.position = "none",
          axis.title.x = element_blank()) +
    
ggsave("plots/indel_diversity_by_class.png", dpi=150)

###############




##### Figure13
##### Derived Allele Frequency spectrum (dsim)
## needs DAF_graph from DAF.r (depends on fdataCl_graphs)

## Temporary duplication of SNPs values to add to both Insertions & Deletions graph
fdistribute.na.in_or_del <- function(dat) {
    rbind(
        transform(subset(dat, in_or_del %in% c("IN", NA)), in_or_del="IN"),
        transform(subset(dat, in_or_del %in% c("DEL", NA)), in_or_del="DEL")
    )
}

# reorder levels
# Convert class3 into factor
melted_DAF_graph <- melt(DAF_graph)
melted_DAF_graph <- transform(melted_DAF_graph, class3 = as.factor(class3))
# reorder levels using factor()
melted_DAF_graph$class3 <- factor(melted_DAF_graph$class3, levels=c("SNPs_non_synonymous",
                                                                    "SNPs_small_intron",
                                                                    "SNPs_synonymous",
                                                                    "CDS_frameshift",
                                                                    "CDS_non_frameshift",
                                                                    "five_prime_UTR",
                                                                    "three_prime_UTR",
                                                                    "small_intron",
                                                                    "long_intron",
                                                                    "intergenic" ))


ggplot(distribute.na.in_or_del(melted_DAF_graph), aes(x=variable, y=value*100,  fill=class3))+
    geom_bar(stat="identity", position="dodge", colour="black")+
    #geom_errorbar(ymax=(value*100)+0.1, ymin=(value*100)-0.1)+
    facet_grid(in_or_del ~ aut_or_X) +
    theme_bw() +
    labs(title="Derived Allele Frequency spectrum (dsim)", x="Derived Allele Frequency (%)", y="Allele count (%)") +
    guides(fill = guide_legend(keywidth = 0.8, keyheight = 0.8, label.theme=element_text(size=12, angle=0))) +
    scale_x_discrete(labels=c("<10%","<20%","<30%","<40%","<50%","<60%","<70%","<80%","<90%","<100%"))+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.title = element_text(size=14),
          strip.text = element_text(size=12)) +
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

rm(melted_DAF_graph)

ggsave("plots/over100bp/PG6_DAF.png", dpi=250, scale=1.5)

#################



