#########################################################################
### Sequential functional dispersion and RaoQ correlations of Machair ###
#########################################################################

### Hannah White 05.05.2021

### Calculates functional dispersion and RaoQ of sub-assemblages with species added sequentially 
### from common to rare and rare to common. Pearson correlations are then calculated 
### for the diversity of each sub-assemblage with the full assemblage.

### This is for use with presence-absence data

### This code should be adjusted depending on the scale under investigation.

set.seed(32)

library(FD)

# get traits data
traits <- read.csv('D:\\IRC_Trinity\\Data\\Machair\\FullTraits.csv', header = TRUE)
traits <- traits[, c('species', 'Canopy.height', 'Seed.mass', 'Leaf.size',  'SLA', 
                     'Insects', 'Selfing', 'Wind',
                     'Seed', 'Vegetative')]

traits$species <- gsub(' ', '.', traits$species)
traits$species <- gsub('-', '.', traits$species)

# change categorical traits to binary factor
traits[,6:10] <- ifelse(traits[,6:10] > 0, 1, 0)


# standardise traits
#traits[,2:5] <- scale(traits[,2:5]) # ONLY DO THIS IF USING EUCLIDEAN DISTANCE

## area of occupancy function
aoo.func <- function(x){length(which(x > 0))} # calculates number of squares cover is > 0

# read in raw data and relabel plot 824Z as 842 Z
species <- read.csv('D:/IRC_Trinity/Data/Machair/Sitexspeciescover.csv', header = TRUE)
names(species) <- gsub('Calastegia', 'Calystegia', names(species))

library(plyr)
species$Plot.ID <- revalue(species$Plot.ID, c('824Z'='842Z'))

# remove composition columns which aren't species
species <- species[,-c(150:158)]

plot.label <- paste(species$Plot.ID, species$Year, species$Scale, species$Refno., species$Refletter., sep = '_')

species <- data.frame(species[,1:6], plot.label, species[,7:149])

## Extract data from 2008 at each scale of 1 x 1 m
species08.4cm <- species[grep('2008_0.04', species$plot.label),]
species08.20cm <- species[grep('2008_0.2', species$plot.label),]
species08.1m <- species[grep('2008_1', species$plot.label),]

## Delete as appropriate for scale under investigation
species.scale <- species08.4cm
species.scale <- species08.20cm
species.scale <- species08.1m

#re-order species names alphabetically
specs <- species.scale[,13:150]
sp.order <- specs[,order(colnames(specs))]

#remove species not present
pp <- which(colSums(sp.order)!=0)
present.sp <- sp.order[,pp]
site.sp <- data.frame(species.scale[,c(2, 5:7)], present.sp) # creates data frame of site info and species % cover

xdim <- dim(site.sp)[2] # number of columns in data frame
n <- xdim - 4 # number of species

site.sp[,5:xdim] <- ifelse(site.sp[,5:xdim] == 0, 0, 1) # change to presence absence 

# calculate geographical rarity of each species at this scale
occ <- apply(site.sp[,5:xdim], 2, FUN = aoo.func)

# have one ordered names just for reference later on

ord.c2r <- occ[rev(order(occ,sample(n, n)))] # sample argument randomises ites
common2rare <- names(ord.c2r)
#ord.r2c <- occ[order(occ, runif(length(occ)))]
#rare2common <- names(ord.r2c)
rare2common <- rev(common2rare)

### Calculate FD for each site
rownames(traits) <- names(sp.order)
traits <- traits[,-1]

# only use traits in site.sp
traits <- traits[rownames(traits) %in% names(site.sp[,5:xdim]),]

## set weights so that binary pollen vector fuzzy coding is equal to continuous traits
w <- c(1, 1, 1, 1, 1/3, 1/3, 1/3, 1/2, 1/2)


#### Common to Rare FD

c2r <- common2rare

c2r.fdis <- array(data = NA, dim = c(475, n))
rownames(c2r.fdis) <- site.sp$plot.label

c2r.raoq <- array(data = NA, dim = c(475, n))
rownames(c2r.raoq) <- site.sp$plot.label
system.time(
  for (sp in 2:n) {
    sub.sp <- c2r[1:sp]
    sub.comm <- site.sp[ ,sub.sp] #create sub-assemblage
    rownames(sub.comm) <- site.sp$plot.label
    sub.traits <- traits[rownames(traits) %in% sub.sp, ] #create smaller trait matrix with species just in sub assemblage
    
    # can only calculate FD for sites with species present
    pop.sites <- sub.comm[rowSums(sub.comm)>0,]
    # change sub.comm columns to alphabetical order
    comm.alpha <- pop.sites[, order(colnames(pop.sites))] 
    
    # calculate functional diversity
    out <- FD::dbFD(sub.traits, comm.alpha, calc.CWM = FALSE, calc.FRic = FALSE, w = w)
    
    out.fdis <- out$FDis
    out.raoq <- out$RaoQ
    
    
    c2r.fdis[names(out.fdis), sp] <- out.fdis
    c2r.raoq[names(out.raoq), sp] <- out.raoq
    
  }
) 

c2r.fdis[,1] <- ifelse(site.sp[,c2r[1]] > 0, 0, NA)
c2r.raoq[,1] <- ifelse(site.sp[,c2r[1]] > 0, 0, NA)


## Pearsons

## Pearson correlations
c2r.pear.fdis <- apply(c2r.fdis[,1:n], 2, FUN = function(x2) cor(cbind(x2, c2r.fdis[,n]), method = 'pearson', use = 'complete.obs')[[2]])
c2r.pear.raoq <- apply(c2r.raoq[,1:n], 2, FUN = function(x2) cor(cbind(x2, c2r.raoq[,n]), method = 'pearson', use = 'complete.obs')[[2]])

c2r.pearcorrs <- data.frame(c2r.pear.fdis, c2r.pear.raoq)


### Rare to common 

r2c <- rare2common

r2c.fdis <- array(data = NA, dim = c(475, n))
rownames(r2c.fdis) <- site.sp$plot.label

r2c.raoq <- array(data = NA, dim = c(475, n))
rownames(r2c.raoq) <- site.sp$plot.label
system.time(
  for (sp in 2:n) {
    sub.sp <- r2c[1:sp]
    sub.comm <- site.sp[ ,sub.sp] #create sub-assemblage
    rownames(sub.comm) <- site.sp$plot.label
    sub.traits <- traits[rownames(traits) %in% sub.sp, ] #create smaller trait matrix with species just in sub assemblage
    
    # can only calculate FD for sites with species present
    pop.sites <- sub.comm[rowSums(sub.comm)>0,]
    # change sub.comm columns to alphabetical order
    comm.alpha <- pop.sites[, order(colnames(pop.sites))] 
    
    # calculate functional diversity
    out <- FD::dbFD(sub.traits, comm.alpha, calc.CWM = FALSE, calc.FRic = FALSE, w = w)
    out.fdis <- out$FDis
    out.raoq <- out$RaoQ
    
    
    r2c.fdis[names(out.fdis), sp] <- out.fdis
    r2c.raoq[names(out.raoq), sp] <- out.raoq
    
  }
) 

r2c.fdis[,1] <- ifelse(site.sp[,r2c[1]] > 0, 0, NA)
r2c.raoq[,1] <- ifelse(site.sp[,r2c[1]] > 0, 0, NA)

### Pearsons
## Pearson correlations
r2c.pear.fdis <- apply(r2c.fdis[,1:n], 2, FUN = function(x2) cor(cbind(x2, r2c.fdis[,n]), method = 'pearson', use = 'complete.obs')[[2]])
r2c.pear.raoq <- apply(r2c.raoq[,1:n], 2, FUN = function(x2) cor(cbind(x2, r2c.raoq[,n]), method = 'pearson', use = 'complete.obs')[[2]])

r2c.pearcorrs <- data.frame(r2c.pear.fdis, r2c.pear.raoq)


###### Calculate information (binomial variance) in subassemblages
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

#### Combine data into single plottable dataframe
subassem <- 1:n

pearFD.df <- data.frame(subassem, 
                        c2r.fdis = c2r.pearcorrs[,1], c2r.raoq = c2r.pearcorrs[,2],
                        r2c.fdis = r2c.pearcorrs[,1], r2c.raoq = r2c.pearcorrs[,2], 
                        sum.pq.c2r, sum.pq.r2c)

############################
###### Make some plots #####
############################

##### Plots in ggplot2 #####

library(ggplot2)
library(reshape2)

## load randomisations for correct scale - DELETE AS APPROPRIATE
load('randomcorrs.c2r4cm.RData')
load('randomcorrs.r2c4cm.RData')

load('randomcorrs.c2r20cm.RData')
load('randomcorrs.r2c20cm.RData')

load('randomcorrs.c2r.RData')
load('randomcorrs.r2c.RData')


subassem <- pearFD.df$subassem
sum.pq.c2r <- pearFD.df$sum.pq.c2r
sum.pq.r2c <- pearFD.df$sum.pq.r2c


### FDis
# Subassemblage 

c2r.fdis.rand <- data.frame(subassem, fdis.corr)
c2r.fdis.rand <- melt(c2r.fdis.rand, id.vars="subassem")

r2c.fdis.rand <- data.frame(subassem, r2c.fdis.corr)
r2c.fdis.rand <- melt(r2c.fdis.rand, id.vars='subassem')


p.subfdis <- ggplot() + geom_line(data = c2r.fdis.rand, col = 'goldenrod', size = 0.2, alpha = 0.1,
                                  aes(x = subassem, y = value, group = variable))
p.subfdis <- p.subfdis + geom_line(data = r2c.fdis.rand, col = 'darkolivegreen2', size = 0.2, alpha = 0.1,
                                   aes(x = subassem, y =value, group = variable))
p.subfdis <- p.subfdis + geom_line(data = pearFD.df, col = 'gold1', size = 1.5,
                                   aes(x = subassem, y = c2r.fdis)) 
p.subfdis <- p.subfdis + geom_line(data = pearFD.df, col = 'forestgreen', size = 1.5,
                                   aes(x = subassem, y = r2c.fdis))
p.subfdis <- p.subfdis + xlab('Subassemblage size') + ylab('Correlation') + 
  theme_classic() + theme(axis.text = element_text(size = 11))


## Cumulative information

c2r.fdis.rand.pq <- data.frame(sum.pq.c2r = pearFD.df$sum.pq.c2r, fdis.corr)
c2r.fdis.rand.pq <- melt(c2r.fdis.rand.pq, id.vars="sum.pq.c2r")

r2c.fdis.rand.pq <- data.frame(sum.pq.r2c = pearFD.df$sum.pq.r2c, r2c.fdis.corr)
r2c.fdis.rand.pq <- melt(r2c.fdis.rand.pq, id.vars='sum.pq.r2c')

p.infofdis <- ggplot() + geom_line(data = c2r.fdis.rand.pq, col = 'goldenrod', size = 0.2, alpha = 0.1,
                                   aes(x = sum.pq.c2r, y = value, group = variable))
p.infofdis <- p.infofdis + geom_line(data = r2c.fdis.rand.pq, col = 'darkolivegreen2', size = 0.2, alpha = 0.1,
                                     aes(x = sum.pq.r2c, y =value, group = variable))
p.infofdis <- p.infofdis + geom_line(data = pearFD.df, size = 1.5,
                                     aes(x = sum.pq.c2r, y = c2r.fdis, col = 'gold1')) + 
  geom_rug(data = pearFD.df, col = 'gold1', aes(x = sum.pq.c2r))
p.infofdis <- p.infofdis + geom_line(data = pearFD.df, size = 1.5,
                                     aes(x = sum.pq.r2c, y = r2c.fdis, col = 'forestgreen')) + 
  geom_rug(data = pearFD.df, col = 'forestgreen', alpha = 0.7, aes(x = sum.pq.r2c))
p.infofdis <- p.infofdis + scale_color_manual(name = '', values = c('gold1' = 'gold1', 'forestgreen' = 'forestgreen'),
                                              labels = c('rare to common', 'common to rare'))
p.infofdis <- p.infofdis + xlab('Cumulative information') + ylab('Correlation') + 
  theme_classic() + theme(axis.text = element_text(size = 11))


### RaoQ
# Subassemblage 

c2r.raoq.rand <- data.frame(subassem, raoq.corr)
c2r.raoq.rand <- melt(c2r.raoq.rand, id.vars="subassem")

r2c.raoq.rand <- data.frame(subassem, r2c.raoq.corr)
r2c.raoq.rand <- melt(r2c.raoq.rand, id.vars='subassem')


p.subraoq <- ggplot() + geom_line(data = c2r.raoq.rand, col = 'goldenrod', size = 0.2, alpha = 0.1,
                                  aes(x = subassem, y = value, group = variable))
p.subraoq <- p.subraoq + geom_line(data = r2c.raoq.rand, col = 'darkolivegreen2', size = 0.2, alpha = 0.1,
                                   aes(x = subassem, y =value, group = variable))
p.subraoq <- p.subraoq + geom_line(data = pearFD.df, col = 'gold1', size = 1.5,
                                   aes(x = subassem, y = c2r.raoq)) 
p.subraoq <- p.subraoq + geom_line(data = pearFD.df, col = 'forestgreen', size = 1.5,
                                   aes(x = subassem, y = r2c.raoq))
p.subraoq <- p.subraoq + xlab('Subassemblage size') + ylab('Correlation') + 
  theme_classic() + theme(axis.text = element_text(size = 11))


## Cumulative information

c2r.raoq.rand.pq <- data.frame(sum.pq.c2r, raoq.corr)
c2r.raoq.rand.pq <- melt(c2r.raoq.rand.pq, id.vars="sum.pq.c2r")

r2c.raoq.rand.pq <- data.frame(sum.pq.r2c, r2c.raoq.corr)
r2c.raoq.rand.pq <- melt(r2c.raoq.rand.pq, id.vars='sum.pq.r2c')

p.inforaoq <- ggplot() + geom_line(data = c2r.raoq.rand.pq, col = 'goldenrod', size = 0.2, alpha = 0.1,
                                   aes(x = sum.pq.c2r, y = value, group = variable))
p.inforaoq <- p.inforaoq + geom_line(data = r2c.raoq.rand.pq, col = 'darkolivegreen2', size = 0.2, alpha = 0.1,
                                     aes(x = sum.pq.r2c, y =value, group = variable))
p.inforaoq <- p.inforaoq + geom_line(data = pearFD.df, size = 1.5,
                                     aes(x = sum.pq.c2r, y = c2r.raoq, col = 'gold1')) + 
  geom_rug(data = pearFD.df, col = 'gold1', aes(x = sum.pq.c2r))
p.inforaoq <- p.inforaoq + geom_line(data = pearFD.df, size = 1.5,
                                     aes(x = sum.pq.r2c, y = r2c.raoq, col = 'forestgreen')) + 
  geom_rug(data = pearFD.df, col = 'forestgreen', alpha = 0.7, aes(x = sum.pq.r2c))
p.inforaoq <- p.inforaoq + scale_color_manual(name = '', values = c('gold1' = 'gold1', 'forestgreen' = 'forestgreen'),
                                              labels = c('rare to common', 'common to rare'))
p.inforaoq <- p.inforaoq + xlab('Cumulative information') + ylab('Correlation') + 
  theme_classic() + theme(axis.text = element_text(size = 11))

