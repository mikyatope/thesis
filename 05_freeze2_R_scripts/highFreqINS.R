#### Select high frequency INSERTIONS

high_freq_INS <- subset(
                        merge(DAF_indels[, c("chr","start","derived","in_or_del","DAF", "type")], data[,c("chr","start","REF","ALT")], all.y=FALSE), 
                        DAF >= 184 & in_or_del == "IN"
                    )


#output
sampled.df <- subset(high_freq_INS, type=="very_good" | type=="longIndel")
sampled.df <- sampled.df[sample(nrow(sampled.df), 50),c("chr","start","start","REF","ALT", "derived", "in_or_del", "DAF")]
sampled.df <- sampled.df[with(sampled.df, order(chr, start)), ]

write.table(
    sampled.df,
    file="blasts/high_freq_INS_HQ.csv",
    sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE
)

rm(sampled.df)