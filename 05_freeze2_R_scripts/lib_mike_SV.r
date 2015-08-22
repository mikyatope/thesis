   
getAndFilterFreeze2 <- function(path) {
  
    dataset  <- read.delim(path, header=F, quote="")
    
    ### Filter raw dataset
    filtered <- data.frame(do.call(rbind, strsplit(as.character(dataset$V8), ";")))   # Split INFO column using semicolons ";" (strsplit returns vectors, rbind/do.cal/data.frame converts it to dataframe)
    filtered <- as.data.frame(sapply(filtered, gsub, pattern=".*=", replacement=""))  # Replace keys from VCF, get only values (use sapply to do it in all columns)
    filtered <- cbind(dataset[c(1:2,4:5)], filtered)                                  # Merge dataframes
    
    rm(dataset)
    
    ##### ADD column headers
    ## INFO column legend:
    ## AC (#valid/homozigous major ALT allele lines), 
    ## RC (#homozigous REF allele lines)
    ## FAC (AC + rand(HETQ) fixed value)
    ## RLN (REF length), 
    ## ALN (ALT length), 
    ## MISS (missing genotypes), 
    ## HET (heterozigous genotypes)
    ## HETQ (high quality het. (>= 99))
    ## PI (k, or nucl. diversity for a single indel)
    cols            <- c("chr", "start", "REF", "ALT", "AC", "RC", "FAC", "RLN", "ALN","MISS", "HET", "HETQ", "PI" )
    names(filtered) <- cols; rm(cols)
    
    #### Refactor column data type
    filtered$chr  <- as.character(filtered$chr)
    #start gets the correct INT type
    filtered$REF  <- as.character(filtered$REF)
    filtered$ALT  <- as.character(filtered$ALT)
    filtered$AC   <- as.numeric(as.character(filtered$AC))
    filtered$RC   <- as.numeric(as.character(filtered$RC))
    filtered$FAC  <- as.numeric(as.character(filtered$FAC))
    filtered$RLN  <- as.numeric(as.character(filtered$RLN))
    filtered$ALN  <- as.numeric(as.character(filtered$ALN))
    filtered$MISS <- as.numeric(as.character(filtered$MISS))
    filtered$HET  <- as.numeric(as.character(filtered$HET))
    filtered$HETQ <- as.numeric(as.character(filtered$HETQ))
    filtered$PI   <- as.numeric(as.character(filtered$PI))
    
    
    
    ##### lenghts 
    #create empty 'size' column and populate with lenght data (ALT - REF allele lenghts)
    filtered["size"] <- NA
    filtered$size    <- filtered$ALN-filtered$RLN
    
    # remove 'zero' lenght indels (previous filtering errors?)
    filtered         <- subset(filtered, size != 0)
    
    return(filtered)   
    
}

getAndFilterFreeze2SNP <- function(path) {
    
    dataset  <- read.delim(path, header=F, quote="")
    
    ### Filter raw dataset
    filtered <- data.frame(do.call(rbind, strsplit(as.character(dataset$V8), ";")))   # Split INFO column using semicolons ";" (strsplit returns vectors, rbind/do.cal/data.frame converts it to dataframe)
    filtered <- as.data.frame(sapply(filtered, gsub, pattern=".*=", replacement=""))  # Replace keys from VCF, get only values (use sapply to do it in all columns)
    filtered <- cbind(dataset[c(1:2,4:5)], filtered)                                  # Merge dataframes
    
    #rm(dataset)
    
    ##### ADD column headers
    ## INFO column legend:
    ## AC (#valid/homozigous major ALT allele lines), 
    ## RC (#homozigous REF allele lines)
    ## FAC (AC + rand(HETQ) fixed value)
    ## RLN (REF length), 
    ## ALN (ALT length), 
    ## MISS (missing genotypes), 
    ## HET (heterozigous genotypes)
    ## HETQ (high quality het. (>= 99))
    ## PI (k, or nucl. diversity for a single indel)
    cols            <- c("chr", "start", "REF", "ALT", "AC", "RC", "FAC", "RLN", "ALN","MISS", "HET", "HETQ", "PI" )
    names(filtered) <- cols; rm(cols)
    
    #### Refactor column data type
    filtered$chr  <- as.character(filtered$chr)
    #start gets the correct INT type
    filtered$REF  <- as.character(filtered$REF)
    filtered$ALT  <- as.character(filtered$ALT)
    filtered$AC   <- as.numeric(as.character(filtered$AC))
    filtered$RC   <- as.numeric(as.character(filtered$RC))
    filtered$FAC  <- as.numeric(as.character(filtered$FAC))
    filtered$RLN  <- as.numeric(as.character(filtered$RLN))
    filtered$ALN  <- as.numeric(as.character(filtered$ALN))
    filtered$MISS <- as.numeric(as.character(filtered$MISS))
    filtered$HET  <- as.numeric(as.character(filtered$HET))
    filtered$HETQ <- as.numeric(as.character(filtered$HETQ))
    filtered$PI   <- as.numeric(as.character(filtered$PI))
    
    
    
    ##### lenghts 
    #create empty 'size' column and populate with lenght data (ALT - REF allele lenghts)
    filtered["size"] <- NA
    filtered$size    <- filtered$ALN-filtered$RLN
    
    # remove 'zero' lenght indels (previous filtering errors?)
    # filtered         <- subset(filtered, size != 0)
    
    return(filtered)   
    
}




getAndFilterFreeze2Class <- function(path) {
    
    dataset  <- read.delim(path, header=F, quote="")
    
    ##### ADD column headers
    ## class/classStart/classEnd/classInfo -> data from BEDtools intersect                  
    cols            <- c("chr", "start", "class", "classStart", "classEnd", "classInfo" )
    names(dataset) <- cols; rm(cols)
    
    #### Refactor column data type
    dataset$chr  <- as.character(dataset$chr)
    #start gets the correct INT type
    dataset$class       <- as.character(dataset$class)
    dataset$classStart  <- as.numeric(as.character(dataset$classStart))                      
    dataset$classEnd    <- as.numeric(as.character(dataset$classEnd)) 
    dataset$classInfo   <- as.character(dataset$classInfo)                      
                      
    ##### lenghts 
    #create empty 'size' column and populate with lenght data for the annotation
    dataset["classSize"] <- NA
    dataset$classSize    <- dataset$classEnd-dataset$classStart+1  # +1!!!!
    
    return(dataset)   
    
}


#### Get r2
getAndFilterFreeze2r2 <- function(path) {
    
    dataset  <- read.delim(path, header=T, quote="")
    
    cols            <- c("chr", "pos1", "pos2", "n_indv", "r2" )
    names(dataset) <- cols; rm(cols)
    
    dataset$chr  <- as.character(dataset$chr)
    dataset$pos1  <- as.numeric(as.character(dataset$pos1))
    dataset$pos2  <- as.numeric(as.character(dataset$pos2))
    dataset$n_indv  <- as.numeric(as.character(dataset$n_indv))
    dataset$r2  <- as.numeric(as.character(dataset$r2))
    
    return(dataset)   
    
}



##### Sliding windows INDELs loop
#Stuff to do each slide --> to be used with a by loop
#input: (slide) of a window coordinates dataset, (dataset) to be cut
sumAndMeans <- function(slide, dataset) {
    
    Xchr   <- slide$chr
    Xstart <- slide$start
    Xend   <- slide$end
    
    wsize  <- Xend-Xstart+1
    slide$wsize <- wsize
    
    #drops <- c("chr", "start", "REF", "ALT")
    slideData<- subset(dataset, chr == Xchr & start >= Xstart & start <= Xend) #[,!(names(dataset) %in% drops)]
    
    slide$numINDELS <- nrow(slideData)
    #slide$numIN <- nrow(subset(slideData, size > 0))
    #slide$numDEL <- nrow(subset(slideData, size < 0))
    
    #slide$affectedSites <- sum(abs(slideData$size))
    #slide$gain <- sum(subset(slideData, size > 0)$size)
    #slide$loss <- sum(subset(slideData, size < 0)$size)
    
    ###PIs
    slide$PI  <- mean(slideData$newPI)
    slide$sPI <- sum(slideData$newPI)
    
    slide$PIindel_adim <- slide$sPI/wsize
    
    return(slide)
}


##### Sliding windows SNP loop
#Stuff to do each slide --> to be used with a by loop
#input: (slide) of a window coordinates dataset, (dataset) to be cut
sumAndMeansSNP <- function(slide, dataset) {
    Xchr   <- slide$chr
    Xstart <- slide$start
    Xend   <- slide$end
    
    wsize  <- Xend-Xstart+1
    slide$wsize <- wsize
    
    drops <- c("chr", "start", "REF", "ALT")
    slideData<- subset(dataset, chr == Xchr & start >= Xstart & start <= Xend)[,!(names(dataset) %in% drops)]
    
    
    ###PIs
    slide$PI <- mean(slideData$jgilPI)
    slide$sPI <- sum(slideData$jgilPI)
    
    slide$PI_SNP_adim <- slide$sPI/wsize
    
    return(slide)
}


##### Sliding windows r2 loop
#Stuff to do each slide --> to be used with a by loop
#input: (slide) of a window coordinates dataset, (dataset) to be cut, 
sumAndMeansR2 <- function(slide, dataset) {
    Xchr   <- slide$chr
    Xstart <- slide$start
    Xend   <- slide$end
    
    wsize  <- Xend-Xstart+1
    slide$wsize <- wsize
    
    slideData<- subset(dataset, chr == Xchr & pos1 >= Xstart & pos2 <= Xend)
    
    ###PIs
    slide$r2<- mean(slideData$r2)
    combins <- combn(slideData$pos1, 2, diff)
    slide$pairws_mean_distance<- mean(combins)
    slide$pairws_median_distance<- median(combins)
    slide$r2__mean_distance <-  slide$r2/slide$pairws_mean_distance
    slide$r2__median_distance <-  slide$r2/slide$pairws_median_distance
    slide$n <- length(slideData$chr)

    return(slide)
}


##### Sliding windows r2 loop
removeDuplicatedChrStart <- function(df) {
    
    # unique duplicated coordinates
    df_unique_dup_coords <- unique( df[duplicated(df[c("chr", "start")]), c("chr", "start")] )
    
    # only duplicated entries (ALL cases, also 1st one)
    df_all_dup_rows <- merge(df_unique_dup_coords, df, by=c("chr","start"), all = FALSE)
    
    # Remove duplicate dataset from whole dataset
    require(sqldf)
    df_real_unique <- sqldf('SELECT * FROM df EXCEPT SELECT * FROM df_all_dup_rows')
    
    return(df_real_unique)
}
