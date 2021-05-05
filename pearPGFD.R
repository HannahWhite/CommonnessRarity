###########################################################
### Common rare correlations of Petchey and Gaston's FD ###
###########################################################

### Hannah White 05.05.2021

### Calculates Petcehy and Gaston's FD of sub-assemblages with species added sequentially 
### from common to rare and rare to common. Pearson correlations are then calculated 
### for the diversity of each sub-assemblage with the full assemblage.

### This code should be adjusted depending on the scale under investigation.


library(vegan)
library(cluster)
library(reshape2)

set.seed(32)

# get traits data
traits <- read.csv('FullTraits.csv', header = TRUE)
traits <- traits[, c('species', 'Canopy.height', 'Seed.mass', 'Leaf.size',  'SLA', 'Insects', 'Selfing', 'Wind', 'Seed', 'Vegetative')]

traits$species <- gsub(' ', '.', traits$species)
traits$species <- gsub('-', '.', traits$species)

# change categorical traits to binary factor
traits[,6:10] <- ifelse(traits[,6:10] > 0, 1, 0)

rownames(traits) <- traits[,1]
traits <- traits[,-1]

## area of occupancy function
aoo.func <- function(x){length(which(x > 0))} # calculates number of squares cover is > 0

# read in raw data and relabel plot 824Z as 842 Z
species <- read.csv('Sitexspeciescover.csv', header = TRUE)
names(species) <- gsub('Calastegia', 'Calystegia', names(species))

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
specs <- species08.1[,13:150]
sp.order <- specs[,order(colnames(specs))]

#remove species not present
pp <- which(colSums(sp.order)!=0)
present.sp <- sp.order[,pp]
site.sp <- data.frame(species08.1[,c(2, 5:7)], present.sp) # creates data frame of site info and species % cover

xdim <- dim(site.sp)[2] # number of columns in data frame
n <- xdim - 4 # number of species

site.sp[,5:xdim] <- ifelse(site.sp[,5:xdim] == 0, 0, 1) # change to presence absence 

# calculate geographical rarity of each species at this scale
occ <- apply(site.sp[,5:xdim], 2, FUN = aoo.func)

# have one ordered names just for reference later on

ord.c2r <- occ[rev(order(occ,sample(n, n)))] # sample argument randomises species with equal aoo
common2rare <- names(ord.c2r)
rare2common <- rev(common2rare)

# only use traits in site.sp
traits <- traits[rownames(traits) %in% names(site.sp[,5:xdim]),]

#### Create functional dendrogram
w = c(1, 1, 1, 1, 1/3, 1/3, 1/3, 1/2, 1/2)

trait.dist <- daisy(traits, metric='gower', weights = w) 
ftree <- hclust(trait.dist, method='average') #UPGMA clustering
plot(ftree)

# full functional diversity
full.fd <- treedive(site.sp[,5:xdim], ftree) 

### Common to rare

c2r <- common2rare

c2r.fd <- array(data = NA, dim = c(475, n))
row.names(c2r.fd) <- site.sp$plot.label

system.time(
  for (sp in 2:n) {
    sub.sp <- c2r[1:sp]
    sub.comm <- site.sp[ ,sub.sp] #create sub-assemblage
    rownames(sub.comm) <- site.sp$plot.label
    sub.traits <- traits[rownames(traits) %in% sub.sp, ] #create smaller trait matrix with species just in sub assemblage
    
    # calculate functional diversity
    out <- treedive(sub.comm, ftree)
    
    c2r.fd[names(out), sp] <- out
    
  }
) 

c2r.fd[,1] <- ifelse(site.sp[,c2r[1]] > 0, 0, NA)

### Rare to common
r2c <- rare2common

r2c.fd <- array(data = NA, dim = c(475, n))
row.names(r2c.fd) <- site.sp$plot.label

system.time(
  for (sp in 2:n) {
    sub.sp <- r2c[1:sp]
    sub.comm <- site.sp[ ,sub.sp] #create sub-assemblage
    rownames(sub.comm) <- site.sp$plot.label
    sub.traits <- traits[rownames(traits) %in% sub.sp, ] #create smaller trait matrix with species just in sub assemblage
    
    # calculate functional diversity
    out <- treedive(sub.comm, ftree)
    
    r2c.fd[names(out), sp] <- out
    
  }
) 

r2c.fd[,1] <- ifelse(site.sp[,r2c[1]] > 0, 0, NA)

########################
##### Correlations #####
########################

## Pearson
c2r.pearfd <- apply(c2r.fd[,1:n], 2, 
                    FUN = function(x, y) cor(x, full.fd, method ='pearson', use = 'complete.obs'))


r2c.pearfd <- apply(r2c.fd[,1:n], 2,
                    FUN = function(x, y) cor(x, full.fd, method = 'pearson', use = 'complete.obs'))


FD.correlations <- data.frame(c2r.pearfd, r2c.pearfd)

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

## Combine into dataframe
subassem <- 1:115

petcheyFD.df <- data.frame(subassem, FD.correlations, sum.pq.c2r, sum.pq.r2c)

### Graphics

# load randomisations for correct spatial scale - DELETE AS APPROPRIATE
load('randomcorrsPetchey.c2r4cm.RData')
load('randomcorrsPetchey.r2c4cm.RData')

load('randomcorrsPetchey.c2r20cm.RData')
load('randomcorrsPetchey.r2c20cm.RData')

load('randomcorrsPetchey.c2r.RData')
load('randomcorrsPetchey.r2c.RData')

subassem <- petcheyFD1m$subassem
sum.pq.c2r <- petcheyFD1m$sum.pq.c2r
sum.pq.r2c <- petcheyFD1m$sum.pq.r2c

### Plots in ggplot
library(ggplot2)
library(reshape2)

## Subassemblage 

c2r.rand <- data.frame(subassem, fdc2r.corr)
c2r.rand <- melt(c2r.rand, id.vars="subassem")

r2c.rand <- data.frame(subassem, fdr2c.corr)
r2c.rand <- melt(r2c.rand, id.vars='subassem')


p.sub <- ggplot() + geom_line(data = c2r.rand, col = 'goldenrod', size = 0.2, alpha = 0.1,
                              aes(x = subassem, y = value, group = variable))
p.sub <- p.sub + geom_line(data = r2c.rand, col = 'darkolivegreen2', size = 0.2, alpha = 0.1,
                           aes(x = subassem, y =value, group = variable))
p.sub <- p.sub + geom_line(data = petcheyFD.df, col = 'gold1', size = 1.5,
                           aes(x = subassem, y = c2r.pearfd)) 
p.sub <- p.sub + geom_line(data = petcheyFD.df, col = 'forestgreen', size = 1.5,
                           aes(x = subassem, y = r2c.pearfd))
p.sub <- p.sub + xlab('Subassemblage size') + ylab('Correlation') + 
  theme_classic() + theme(axis.text = element_text(size = 11))


## Cumulative information

c2r.rand.pq <- data.frame(sum.pq.c2r, fdc2r.corr)
c2r.rand.pq <- melt(c2r.rand.pq, id.vars="sum.pq.c2r")

r2c.rand.pq <- data.frame(sum.pq.r2c, fdr2c.corr)
r2c.rand.pq <- melt(r2c.rand.pq, id.vars='sum.pq.r2c')

p.info <- ggplot() + geom_line(data = c2r.rand.pq, col = 'goldenrod', size = 0.2, alpha = 0.1,
                               aes(x = sum.pq.c2r, y = value, group = variable))
p.info <- p.info + geom_line(data = r2c.rand.pq, col = 'darkolivegreen2', size = 0.2, alpha = 0.1,
                             aes(x = sum.pq.r2c, y =value, group = variable))
p.info <- p.info + geom_line(data = petcheyFD.df, size = 1.5,
                             aes(x = sum.pq.c2r, y = c2r.pearfd, col = 'gold1')) + 
  geom_rug(data = petcheyFD.df, col = 'gold1', aes(x = sum.pq.c2r))
p.info <- p.info + geom_line(data = petcheyFD.df, size = 1.5,
                             aes(x = sum.pq.r2c, y = r2c.pearfd, col = 'forestgreen')) + 
  geom_rug(data = petcheyFD.df, col = 'forestgreen', alpha = 0.7, aes(x = sum.pq.r2c))
p.info <- p.info + scale_color_manual(name = '', values = c('gold1' = 'gold1', 'forestgreen' = 'forestgreen'),
                                      labels = c('rare to common', 'common to rare'))
p.info <- p.info + xlab('Cumulative information') + ylab('Correlation') + 
  theme_classic() + theme(axis.text = element_text(size = 11))



