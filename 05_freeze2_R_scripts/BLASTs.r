
######## Functions

####### Main function. Read CSV files with REF/ALT BLAST results

readBLASTs <- function(path_REF, path_ALT) {
    
    raw_REF_blasts  <- read.delim(path_REF, header=F, quote="")
    raw_ALT_blasts  <- read.delim(path_ALT, header=F, quote="")
        
    cols_ALT            <- c("chr", "start", "len_indel_ALT", "num_hits_ALT", "num_hsps_ALT", "evalue_ALT", "qstart_ALT", "qend_ALT", "gaps_ALT", "ident_ALT", "posit_ALT", "total_aln_ALT" )
    cols_REF            <- c("chr", "start", "len_indel_REF", "num_hits_REF", "num_hsps_REF", "evalue_REF", "qstart_REF", "qend_REF", "gaps_REF", "ident_REF", "posit_REF", "total_aln_REF" )
    names(raw_ALT_blasts) <- cols_ALT; rm(cols_ALT)
    names(raw_REF_blasts) <- cols_REF; rm(cols_REF)
    
    
    #### REMOVE multiallelic residuals (duplicated). Make sure we only work with biallelic.
    ## get duplicate pairs chr-start
    dup_REF_coordinates <- raw_REF_blasts[duplicated(raw_REF_blasts[c("chr","start")]),c("chr","start")]
    dup_ALT_coordinates <- raw_ALT_blasts[duplicated(raw_ALT_blasts[c("chr","start")]),c("chr","start")]
    # dup_ORIG_coordinates <- ORIGINALdataset[duplicated(ORIGINALdataset[c("V1","V2")]),c("V1","V2")]
    
    ### Do the same for REF and ALT
    if(nrow(dup_REF_coordinates)>0 | nrow(dup_ALT_coordinates)>0){
        
        if(nrow(dup_REF_coordinates)>0){
            ## ALL duplicated rows (find/merge duplicate coordinates with whole dataset)
            merged_REF_coord <- merge(dup_REF_coordinates, raw_REF_blasts, by=c("chr","start"), all.x = FALSE)        
            
            ## ALL UNIQUE rows (exclude ALL duplicate rows, even the first one, from whole dataset)
            library(sqldf)
            uniq_REF_blasts <- sqldf('SELECT * FROM raw_REF_blasts EXCEPT SELECT * FROM merged_REF_coord')
        } else {
            uniq_REF_blasts <- raw_REF_blasts
        }
    
        if(nrow(dup_ALT_coordinates)>0){
            ## ALL duplicated rows (find/merge duplicate coordinates with whole dataset)
            merged_ALT_coord <- merge(dup_ALT_coordinates, raw_ALT_blasts, by=c("chr","start"), all.x = FALSE)
                        
            ## ALL UNIQUE rows (exclude ALL duplicate rows, even the first one, from whole dataset)
            library(sqldf)
            uniq_ALT_blasts <- sqldf('SELECT * FROM raw_ALT_blasts EXCEPT SELECT * FROM merged_ALT_coord')
        } else {
            uniq_ALT_blasts <- raw_ALT_blasts
        }
        
    } else {
        uniq_REF_blasts <- raw_REF_blasts
        uniq_ALT_blasts <- raw_ALT_blasts
    }
    
    ## MERGE REF/ALT
    merged_BLASTS <- merge(uniq_REF_blasts, uniq_ALT_blasts, by=c("chr","start"), all = FALSE)
    
    ## remove intermediate data frames
    #rm(merged_ALT_coord, merged_REF_coord, raw_REF_blasts, raw_ALT_blasts, dup_ALT_coordinates, merged_REF_coord, uniq_ALT_blasts, uniq_REF_blasts)
    
    
    ## Check each blast
    
    # blastSAMPLE <- merged_BLASTS[sample(1:nrow(merged_BLASTS), 1000, replace=FALSE),]
    # blastSAMPLEr3 <- data.frame(do.call("rbind",  by(blastSAMPLE, 1:nrow(blastSAMPLE), checkBLAST )))
    
    result <- data.frame(do.call(  "rbind",  by(merged_BLASTS, 1:nrow(merged_BLASTS), checkBLAST )  ))
    #rm(merged_BLASTS)
    return(result)
}



####### Check BLASTs and decide derived allele

checkBLAST <- function(row) {
    row$hits_stat  <- NA        
    row$hsps_stat  <- NA
    row$type       <- NA
    row$blast_stat <- NA
    row$derived    <- NA
    row$in_or_del  <- NA
    
    # Normal cases
    if(row$num_hits_REF == 1 & row$num_hits_ALT == 1 ) {
        
        row$hits_stat <- "unique"
        if (row$num_hsps_REF == 1 & row$num_hsps_ALT == 1 ) { 
            row$hsps_stat <- "unique" 
        } else if (row$num_hsps_REF > 1 | row$num_hsps_ALT > 1 ) {
            row$hsps_stat <- "multi"
        }
        
        # Check if indel between coordinates, and diference of gaps equal indel length qstart_REF qend_REF gaps_REF
        if ( (row$qstart_REF < 100 & row$qend_REF > (100+row$len_indel_REF) ) &
             (row$qstart_ALT < 100 & row$qend_ALT > (100+row$len_indel_ALT) ) & 
             (row$gaps_REF != row$gaps_ALT) ) {
            
            # Maybe we have good matches!
            # Check small evals
            # Check 50& or more of flanking regions (remember qstart/end are relative coordinates)
            if (  row$evalue_REF < 1e-10 & 
                  row$evalue_ALT < 1e-10 &
                 (row$qstart_REF <=50 & row$qend_REF >= (150+row$len_indel_REF) ) &
                 (row$qstart_ALT <=50 & row$qend_ALT >= (150+row$len_indel_ALT) ) 
               ) {
                row$blast_stat <- "OK"
                row$type       <- ifelse(row$hsps_stat == "unique", "good", "normal" )
                
                # check if almost perfect BLAST align.
                if ((row$qstart_REF == 1 & row$qend_REF == (200+row$len_indel_REF) ) &
                    (row$qstart_ALT == 1 & row$qend_ALT == (200+row$len_indel_ALT) ) &
                    (abs(row$gaps_REF - row$gaps_ALT) == max(row$len_indel_REF, row$len_indel_ALT)-1 ) &
                    (row$hsps_stat == "unique") ) {
                    
                    row$type       <- "very_good"
                }
                
                # Check if derived allele is REF or ALT
                if (row$gaps_REF > row$gaps_ALT) {
                    row$derived    <- "REF"
                    row$in_or_del  <- ifelse(row$len_indel_REF == 1, "DEL", "IN")
                } else if (row$gaps_REF < row$gaps_ALT) {
                    row$derived    <- "ALT"
                    row$in_or_del  <- ifelse(row$len_indel_ALT == 1, "DEL", "IN")
                }
            }
        }
        
        
    # Long indels create multiple hits/hsps    
    } else if( (row$len_indel_REF >= 25 | row$len_indel_ALT >= 25 ) 
               & (row$num_hits_REF >= 1 | row$num_hits_ALT >= 1) ) { 
        
        row$hits_stat <- "multi"
        row$type      <- "longIndel"

        
        # Check for only one alignment REF or ALT  
        # Check if derived allele is REF or ALT
          
        if (is.na(row$evalue_ALT)) {
            if( (row$qstart_REF < (100+0.2*row$len_indel_REF) 
                 & row$qend_REF > (100+0.8*row$len_indel_REF) )  
                 & ( row$evalue_REF < 1e-10 ) ) {
                
                row$blast_stat <- "OK"
                row$derived     <- "REF"
                row$in_or_del  <- ifelse(row$len_indel_REF == 1, "DEL", "IN") 
            }

        } else if (is.na(row$evalue_REF)) {
            if ( (row$qstart_ALT < (100+0.2*row$len_indel_ALT) 
                  & row$qend_ALT > (100+0.8*row$len_indel_ALT) )  
                  & ( row$evalue_ALT < 1e-10 ) ) {
            
                row$blast_stat <- "OK"
                row$derived     <- "ALT"
                row$in_or_del  <- ifelse(row$len_indel_ALT == 1, "DEL", "IN")
            }
        }
          
    # others    
    } else if(row$num_hits_REF > 1 | row$num_hits_ALT > 1) {  
        row$hits_stat <- "multi"       
    } 
    
    return(row[c("chr","start","qstart_REF", "qend_REF", "len_indel_REF", "gaps_REF", "evalue_REF", "qstart_ALT", "qend_ALT", "len_indel_ALT", "gaps_ALT", "evalue_ALT", "hits_stat", "hsps_stat", "type", "blast_stat","derived", "in_or_del")])
}


    
##########
## MAIN ##
##########


##### NonClass Indels
system.time (
    blasts_dyak <- readBLASTs("blasts/results/noclass_indels/REF_dyak_blasts", "blasts/results/noclass_indels/ALT_dyak_blasts")    
)

system.time (
    blasts_dsim <- readBLASTs("blasts/results/noclass_indels/REF_dsim_blasts", "blasts/results/noclass_indels/ALT_dsim_blasts")    
)


##### MS Indels
system.time (
    blasts_MS_dyak <- readBLASTs("blasts/results/MS_indels/REF_dyak_blasts", "blasts/results/MS_indels/ALT_dyak_blasts")    
)

system.time (
    blasts_MS_dsim <- readBLASTs("blasts/results/MS_indels/REF_dsim_blasts", "blasts/results/MS_indels/ALT_dsim_blasts")    
)


##### TR Indels
system.time (
    blasts_TR_dyak <- readBLASTs("blasts/results/TR_indels/REF_dyak_blasts", "blasts/results/TR_indels/ALT_dyak_blasts")    
)

system.time (
    blasts_TR_dsim <- readBLASTs("blasts/results/TR_indels/REF_dsim_blasts", "blasts/results/TR_indels/ALT_dsim_blasts")    
)









# remove 'false' indels? (REF size = ALT size)
blasts_dyak$size <- abs(blasts_dyak$len_indel_REF - blasts_dyak$len_indel_ALT)
blasts_dsim$size <- abs(blasts_dsim$len_indel_REF - blasts_dsim$len_indel_ALT)

blasts_dyak <- subset(blasts_dyak, size > 0)
blasts_dsim <- subset(blasts_dsim, size > 0)

# Sort
blasts_dyak <- blasts_dyak[with(blasts_dyak, order(chr, start)), ]
row.names(blasts_dyak)<-NULL #reset row.names
blasts_dsim <- blasts_dsim[with(blasts_dsim, order(chr, start)), ]
row.names(blasts_dsim)<-NULL #reset row.names


## general Counts
summary.factor( blasts_dsim$type ) 
# good longIndel    normal very_good      NA's 
#    200701      8911     13008    184276    183833 
summary.factor( blasts_dsim$blast_stat )
# OK   NA's 
# 399822 190907 
summary.factor( blasts_dsim$derived ) 
# ALT    REF   NA's 
# 290703 109119 190907 
summary.factor( blasts_dsim$in_or_del )
# DEL     IN   NA's 
# 225372 174450 190907 

summary.factor( blasts_dyak$type ) 
# good longIndel    normal very_good      NA's 
#    160617      9792     12082     69539    338699 
summary.factor( blasts_dyak$blast_stat ) 
# OK   NA's 
# 245291 345438 
summary.factor( blasts_dyak$derived )
# ALT    REF   NA's 
# 171838  73453 345438 
summary.factor( blasts_dyak$in_or_del )
# DEL     IN   NA's 
# 137347 107944 345438


## Count "OK" by chr
table(subset(blasts_dsim, blast_stat == "OK")$chr)
# 2L    2R    3L    3R     4     X 
# 87853 69808 87425 95544   460 58732 

table(subset(blasts_dyak, blast_stat == "OK")$chr)
# 2L    2R    3L    3R     4     X 
# 55309 45708 51595 59158   162 33359 


# Size distributions
ggplot(subset(blasts_dsim[c("size", "blast_stat", "in_or_del")], blast_stat == "OK"), aes(size))+
    geom_histogram(binwidth=1, origin=0.5)+
    coord_cartesian(xlim = c(0,50) )+
    facet_grid(in_or_del ~.)+
    labs(title="Dsim INDEL size distribution") 


ggplot(subset(blasts_dyak[c("size", "blast_stat", "in_or_del")], blast_stat == "OK"), aes(size))+
    geom_histogram(binwidth=1)+
    coord_cartesian(xlim = c(0,50) )+
    facet_grid(in_or_del ~.)+
    labs(title="Dyak INDEL size distribution") 


# Size distributions by type
ggplot(subset(blasts_dsim[c("size", "blast_stat", "type", "in_or_del")], blast_stat == "OK"), aes(size))+
    geom_histogram(binwidth=1)+
    coord_cartesian(xlim = c(0,50) )+
    facet_grid(in_or_del ~ type)+
    labs(title="Dsim INDEL size distribution") 


ggplot(subset(blasts_dyak[c("size", "blast_stat", "type", "in_or_del")], blast_stat == "OK"), aes(size))+
    geom_histogram(binwidth=1)+
    coord_cartesian(xlim = c(0,50) )+
    facet_grid(in_or_del ~ type)+
    labs(title="Dyak INDEL size distribution") 


