###################################################
### Sequential correlations of species richness ###
###################################################

### Hannah White 05.05.2021

### Code looks at correlations between spatial patterns of species richness
### in the Machair system from 2008
### Correlation method used is Pearson

### Code can be adjusted depending on spatial scale under investigation

set.seed(32)

library(abind)

## area of occupancy function
aoo.func <- function(x){length(which(x > 0))} # calculates number of squares cover is > 0

# read in raw data and relabel plot 824Z as 842 Z
species <- read.csv('Sitexspeciescover.csv', header = TRUE)
library(plyr)
species$Plot.ID <- revalue(species$Plot.ID, c('824Z'='842Z'))

# remove composition columns which aren't species
species <- species[,-c(150:158)]

plot.label <- paste(species$Plot.ID, species$Year, species$Scale, species$Refno., species$Refletter., sep = '_')

species <- data.frame(species[,1:6], plot.label, species[,7:149])

## Extract data from 2008 at scale of 1 x 1 m
species08.4cm <- species[grep('2008_0.04', species$plot.label),]
species08.20cm <- species[grep('2008_0.2', species$plot.label),]
species08.1m <- species[grep('2008_1', species$plot.label),]

## Delete as appropriate for scale under investigation
species.scale <- species08.4cm
species.scale <- species08.20cm
species.scale <- species08.1m

#re-order species names alphabetically
sp <- species.scale[,13:150]
sp.order <- sp[,order(colnames(sp))]

#remove species not present
pp <- which(colSums(sp.order)!=0)
present.sp <- sp.order[,pp]
site.sp <- data.frame(species.scale[,c(2, 5:7)], present.sp) # creates data frame of site info and species % cover

xdim <- dim(site.sp)[2] # number of columns in data frame
n <- xdim - 4 # number of species

## Change to presence-absence (didn't do this in loop as might want to use % cover later on)
site.sp[,5:xdim] <- ifelse(site.sp[,5:xdim] == 0, 0, 1) # change to presence absence 

# calculate geographical rarity of each species at this scale
occ <- apply(site.sp[,5:xdim], 2, FUN = aoo.func)

# have one ordered names just for reference later on
ord.c2r <- occ[rev(order(occ))]
common2rare <- names(ord.c2r)
ord.r2c <- occ[order(occ)]
rare2common <- names(ord.r2c)

## Calculate full species richness at each site 
full.sr <- rowSums(site.sp[,5:xdim])

##### Calculate species richnesses for subsets
### Common to rare

c2r <- common2rare

c2r.sr <- array(data = 0, dim = c(475, n))
row.names(c2r.sr) <- site.sp$plot.label

system.time(
  for (sp in 2:n) {
    sub.sp <- c2r[1:sp]
    sub.comm <- site.sp[ ,sub.sp] #create sub-assemblage
    
    # calculate species richness
    out <- rowSums(sub.comm)
    
    c2r.sr[, sp] <- out
    
  }
) 

c2r.sr[,1] <- ifelse(site.sp[,c2r[1]] > 0, 1, 0) # sets species richness of first column

### Rare to common

r2c <- rare2common

r2c.sr <- array(data = 0, dim = c(475, n))
row.names(r2c.sr) <- site.sp$plot.label

system.time(
  for (sp in 2:n) {
    sub.sp <- r2c[1:sp]
    sub.comm <- site.sp[ ,sub.sp] #create sub-assemblage
    
    # calculate species richness
    out <- rowSums(sub.comm)
    
    r2c.sr[, sp] <- out
    
  }
) 

r2c.sr[,1] <- ifelse(site.sp[,r2c[1]] > 0, 1, 0) # sets species richness of first column

############################
### Pearson correlations ###
############################

c2r.pear.sr <- apply(c2r.sr, 2, FUN = function(x2) cor(cbind(x2, c2r.sr[,n]), method = 'pearson', use = 'complete.obs')[[2]])
r2c.pear.sr <- apply(r2c.sr, 2, FUN = function(x2) cor(cbind(x2, r2c.sr[,n]), method = 'pearson', use = 'complete.obs')[[2]])

#### Subassemblage
subassem <- 1:n

pearcorrsSR <- data.frame(subassem, c2r.pear.sr, r2c.pear.sr)

### Calculate information (binomial variance) in subassemblages
#####Calculate pq
sp.p  <- apply(site.sp[,5:xdim], 2, FUN = function(x) (aoo.func(x))/475) # range size as proportion of total study area
sp.q <- 1-sp.p
pq <- sp.p*sp.q # expected binomial variance for each species
rm(sp.p, sp.q)

# common to rare accumulated information
sum.pq.c2r <- rep(NA, n)

for (i in 1:n){
  # order of equally ranked species does not matter as they will have same information
  sub.sp <- common2rare[1:i] # names of species ranked common to rare
  sub.pq <- pq[sub.sp] #create sub-assemblage
  sum.pq.c2r[i] <- sum(sub.pq)
}

# rare to common accumulated information
sum.pq.r2c <- rep(NA, n)

for (i in 1:n){
  # order of equally ranked species does not matter as they will have same information
  sub.sp <- rare2common[1:i] # names of species ranked rare to common
  sub.pq <- pq[sub.sp]
  sum.pq.r2c[i] <- sum(sub.pq)
}

# add to dataframe
SR.correlations <- data.frame(subassem, sum.pq.c2r, sum.pq.r2c, c2r.pear.sr, r2c.pear.sr)

###################################################################################
### Add randomisations where order species are added to subassemblage is random ###
###################################################################################

rand.name <- replicate(1000, sample(common2rare, length(common2rare), replace=FALSE)) # 1000 random orders of species names

##initiate empty list
rand.list <- list()

###create list of random species vectors
for(i in 1:1000){
  tmp <- rand.name[,i]
  rand.list[[i]] <- tmp
}

## function
sr.func <- function(x, y){
  sr.tmp <- array(data = 0, dim = c(475, n))
  for (sp in 2:n){
    sub.sp <- x[1:sp]
    sub.comm <- y[ ,sub.sp] #create sub-assemblage
    
    out <- apply(sub.comm, 1, FUN = aoo.func)
    sr.tmp[,sp] <- out
  }
  return(sr.tmp)
}

system.time(
  sr.list <- lapply(rand.list, function(x) sr.func(x, y = site.sp))
) 

sr.rand <- abind(sr.list, along = 3)

srrand.corr <- array(data = NA, dim = c(n, 1000))

for(rr in 1:1000){ # about 3 seconds per iteration
  sr.tmp <- sr.rand[,,rr]
  srrand.corr[,rr] <- apply(sr.tmp, 2, FUN = function(x2) cor(cbind(x2, sr.tmp[, n]), method = 'pearson', use = 'complete.obs')[[2]])
  print(rr)
}


#### Get binomial variance for random assemblages
rand.pq <- array(data = NA, dim = c(n, 1000))

for (rr in 1:1000){ # fills empty array with binomial variances of subassemblages built from randomisations
  rand.rarity <- rand.name[,rr] # takes a random order of species
  
  for (sp in 1:n){
    sub.sp <- rand.rarity[1:sp]
    sub.pq <- pq[sub.sp] #create sub-assemblage info
    sum.pq.rand <- sum(sub.pq)
    rand.pq[sp, rr] <- sum.pq.rand
  }
}


#################
##### Plots #####
#################

## Subassemblage 

rand.sub <- data.frame(subassem = SR.correlations$subassem, srrand.corr)
rand.sub <- melt(rand.sub, id.vars="subassem")

p.sub <- ggplot() + geom_line(data = rand.sub, col = 'lightsteelblue4', size = 0.2, alpha = 0.1,
                              aes(x = subassem, y = value, group = variable))
p.sub <- p.sub + geom_line(data = SR.correlations, col = 'gold1', size = 1.5,
                           aes(x = subassem, y = c2r.pear.sr))
p.sub <- p.sub + geom_line(data = SR.correlations, col = 'forestgreen', size = 1.5,
                           aes(x = subassem, y = r2c.pear.sr))
p.sub <- p.sub + xlab('Subassemblage size') + ylab('Correlation') + 
  theme_classic() + theme(axis.text = element_text(size = 11))



## Cumulative information
info <- data.frame(subassem = SR.correlations$subassem, rand.pq)
sumpq.melt <- melt(info, id.vars = 'subassem')
names(sumpq.melt) <- c('subassem', 'randomisation', 'sum.pq')

corrs <- data.frame(subassem = SR.correlations$subassem, srrand.corr)
corrs.melt <- melt(corrs, id.vars = 'subassem')
names(corrs.melt) <- c('subassem', 'randomisation', 'pear')

random.df <- merge(sumpq.melt, corrs.melt, by = c('subassem', 'randomisation'), all = TRUE)

p.info <- ggplot() + geom_line(data = random.df, col = 'lightsteelblue4', size = 0.2, alpha = 0.1,
                               aes(x = sum.pq, y = pear, group = randomisation))
p.info <- p.info + geom_line(data = SR.correlations, size = 1.5,
                             aes(x = sum.pq.c2r, y = c2r.pear.sr, col = 'gold1')) +
  geom_rug(data = SR.correlations, col = 'gold1', aes(x = sum.pq.c2r))
p.info <- p.info + geom_line(data = SR.correlations, size = 1.5,
                             aes(x = sum.pq.r2c, y = r2c.pear.sr, col = 'forestgreen')) +
  geom_rug(data = SR.correlations, col = 'forestgreen', alpha = 0.7, aes(x = sum.pq.r2c))
p.info <- p.info + scale_color_manual(name='', values = c('gold1' = 'gold1', 'forestgreen' = 'forestgreen'), labels = c('rare to common', 'common to rare'))
p.info <- p.info + xlab('Cumulative information') + ylab('Correlation') + 
  theme_classic() + theme(axis.text = element_text(size = 11))

