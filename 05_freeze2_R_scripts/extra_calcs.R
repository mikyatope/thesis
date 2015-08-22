########### EXTRA CALCS 


################### Size distribution by chr of FixedNum indels


ggplot(subset(fixNum_INDELs, in_or_del!= "NA" & size < 50), aes(x=size, fill=in_or_del)) + geom_histogram(binwidth=1, origin=-0.5) + 
    facet_grid(chr ~ in_or_del, scales="free_y") + 
    scale_y_continuous(expand = c(0,0)) + # x axis starts at y zero
    scale_x_continuous(expand = c(0,0)) + # y axis starts at x zero
    #labs(title = "INDEL size distribution per chromosome") + 
    theme_bw() + 
    theme(legend.position="none") +
    labs(x="size (bp)")





################## PI indel by chr of FixNUM indels


ddply(subset(fixNum_INDELs, in_or_del != "NA"), .(chr, in_or_del), function(i)  sum(i$newPI)  )