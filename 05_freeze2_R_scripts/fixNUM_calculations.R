
###### FIX NUM calcs!!!

ggplot((fdataCl_graphs), aes(x=(RC+FAC))) + 
    geom_histogram(binwidth=1, origin=-0.5) +
    geom_vline(xintercept = 174) +
    geom_vline(xintercept = 160) +
    geom_vline(xintercept = 123) 

# 90% !!
prop.table(c(nrow(subset(fdataCl_graphs, RC+FAC > 173)), 
             nrow(fdataCl_graphs)-nrow(subset(fdataCl_graphs, RC+FAC > 165))
))

# 95% !! (94.9876)
prop.table(c(nrow(subset(fdataCl_graphs, RC+FAC > 160)), 
             nrow(fdataCl_graphs)-nrow(subset(fdataCl_graphs, RC+FAC > 149))
))

# 99% !!
prop.table(c(nrow(subset(fdataCl_graphs, RC+FAC > 127)), 
             nrow(fdataCl_graphs)-nrow(subset(fdataCl_graphs, RC+FAC > 105))
))



###### FIX NUM SNPs!

ggplot((fSNPdata), aes(x=(refc+altc))) + 
    geom_histogram(binwidth=1, origin=-0.5) 


# 90% (.89572)
prop.table(c(nrow(subset(fSNPdata, refc+altc > 188)), 
             nrow(fSNPdata)-nrow(subset(fSNPdata, refc+altc  > 188))
))

# 95 %
prop.table(c(nrow(subset(fSNPdata, refc+altc > 183)), 
             nrow(fSNPdata)-nrow(subset(fSNPdata, refc+altc  > 183))
))

# 99 %
prop.table(c(nrow(subset(fSNPdata, refc+altc > 152)), 
             nrow(fSNPdata)-nrow(subset(fSNPdata, refc+altc  > 152))
))
