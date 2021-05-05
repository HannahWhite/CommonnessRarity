#################################################################
### Aggregate to 5 x 5m and calculate Petchey and Gaston's FD ###
#################################################################

### Hannah White 05.05.2021

### pearPGFD.R adjusted so that data is aggregated to 5 m x 5 m and correlations
### applied across 19 plots at this spatial scale

library(ggplot2)
library(cluster)
library(vegan)
library(reshape2)

set.seed(32)

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
species08.1 <- species[grep('2008_1', species$plot.label),]
nplot <- length(unique(species08.1$Plot.ID))

species08.5 <- array(data = NA, dim = c(nplot, 138)) # 138 species (cols 13:150)

for (i in 1:nplot){
  plot <- unique(species08.1$Plot.ID)[i] # gets unique plot label
  sp.temp <- species08.1[species08.1$Plot.ID == plot, 13:150] # temp dataframe with records from single plot
  
  species08.5[i,] <- colSums(sp.temp) # gets total relative abundance measure (not % cover though)
}

species08.5 <- data.frame(unique(species08.1$Plot.ID), species08.5)
names(species08.5) <- c('Plot.ID', names(species08.1)[13:150])

#re-order species names alphabetically
specs <- species08.5[,2:139]
sp.order <- specs[,order(colnames(specs))]


#remove species not present
present.sp <- sp.order[,which(colSums(sp.order)!=0)]
site.sp <- data.frame(species08.5[1], present.sp) # creates data frame of site info and species % cover

## Change to presence-absence (didn't do this in loop as might want to use % cover later on)
site.sp[,2:116] <- ifelse(site.sp[,2:116] == 0, 0, 1) # change to presence absence 

rm(present.sp, specs, species, i)


### Load in and sort out traits
traits <- read.csv('D:\\IRC_Trinity\\Data\\Machair\\FullTraits.csv', header = TRUE)
traits <- traits[, c('species', 'Canopy.height', 'Seed.mass', 'Leaf.size',  'SLA', 
                     'Insects', 'Selfing', 'Wind',
                     'Seed', 'Vegetative')]

traits$species <- gsub(' ', '.', traits$species)
traits$species <- gsub('-', '.', traits$species)

# change categorical traits to binary factor
traits[,6:10] <- ifelse(traits[,6:10] > 0, 1, 0)

### Common and rare calculations
# calculate geographical rarity of each species at this scale
occ <- apply(site.sp[,2:116], 2, FUN = aoo.func)

# have one ordered names just for reference later on

ord.c2r <- occ[rev(order(occ,sample(115, 115)))] # sample argument randomises ites
common2rare <- names(ord.c2r)
rare2common <- rev(common2rare)

### Calculate FD for each site
rownames(traits) <- names(sp.order)
traits <- traits[,-1]

# only use traits in site.sp
traits <- traits[rownames(traits) %in% names(site.sp[,2:116]),]

## set weights so that binary pollen vector fuzzy coding is equal to continuous traits
w <- c(1, 1, 1, 1, 1/3, 1/3, 1/3, 1/2, 1/2)


#### Create functional dendrogram

trait.dist <- daisy(traits, metric='gower', weights = w) 
ftree <- hclust(trait.dist, method='average') #UPGMA clustering
plot(ftree)

# full functional diversity
full.fd <- treedive(site.sp[,2:116], ftree) 

### Common to rare

c2r <- common2rare

c2r.fd <- array(data = NA, dim = c(19, 115))
row.names(c2r.fd) <- site.sp$Plot.ID

system.time(
  for (sp in 2:115) {
    sub.sp <- c2r[1:sp]
    sub.comm <- site.sp[ ,sub.sp] #create sub-assemblage
    rownames(sub.comm) <- site.sp$Plot.ID
    sub.traits <- traits[rownames(traits) %in% sub.sp, ] #create smaller trait matrix with species just in sub assemblage
    
    # calculate functional diversity
    out <- treedive(sub.comm, ftree)
    
    c2r.fd[names(out), sp] <- out
    
  }
) 

c2r.fd[,1] <- ifelse(site.sp[,c2r[1]] > 0, 0, NA)

### Rare to common
r2c <- rare2common

r2c.fd <- array(data = NA, dim = c(19, 115))
row.names(r2c.fd) <- site.sp$Plot.ID

system.time(
  for (sp in 2:115) {
    sub.sp <- r2c[1:sp]
    sub.comm <- site.sp[ ,sub.sp] #create sub-assemblage
    rownames(sub.comm) <- site.sp$Plot.ID
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

## pearson
c2r.pearfd <- apply(c2r.fd[,1:115], 2, 
                    FUN = function(x, y) cor(x, full.fd, method ='pearson', use = 'complete.obs'))
#c2r.pearfd <- c(NA, c2r.pearfd)

r2c.pearfd <- apply(r2c.fd[,1:115], 2,
                    FUN = function(x, y) cor(x, full.fd, method = 'pearson', use = 'complete.obs'))
#r2c.pearfd <- c(NA, r2c.pearfd)

FD.correlations <- data.frame(c2r.pearfd, r2c.pearfd)


###### Calculate information (binomial variance) in subassemblages
#####Calculate pq
sp.p  <- apply(site.sp[,2:116], 2, FUN = function(x) (aoo.func(x))/19) # range size as proportion of total study area
sp.q <- 1-sp.p
pq <- sp.p*sp.q # expected binomial variance for each species
rm(sp.p, sp.q)

# common to rare accumulated information
sum.pq.c2r <- rep(NA, 115)

for (i in 1:115){
  # order of equally ranked species does not matter as they will have same information
  sub.sp <- common2rare[1:i] # names of species ranked common to rare
  sub.pq <- pq[sub.sp] #create sub-assemblage
  sum.pq.c2r[i] <- sum(sub.pq)
}

# rare to common accumulated information
sum.pq.r2c <- rep(NA, 115)

for (i in 1:115){
  # order of equally ranked species does not matter as they will have same information
  sub.sp <- rare2common[1:i] # names of species ranked rare to common
  sub.pq <- pq[sub.sp]
  sum.pq.r2c[i] <- sum(sub.pq)
}

## Combine into dataframe

subassem <- 1:115

petcheyFD.df <- data.frame(subassem, FD.correlations, sum.pq.c2r, sum.pq.r2c)
petcheyFD5m <- petcheyFD.df


### Graphics

# load randomisations

load('randomcorrsPetchey.c2r5m.RData')
load('randomcorrsPetchey.r2c5m.RData')

subassem <- petcheyFD5m$subassem
sum.pq.c2r <- petcheyFD5m$sum.pq.c2r
sum.pq.r2c <- petcheyFD5m$sum.pq.r2c

### Plots in ggplot

## Subassemblage 

c2r.rand <- data.frame(subassem, fdc2r.corr)
c2r.rand <- melt(c2r.rand, id.vars="subassem")

r2c.rand <- data.frame(subassem, fdr2c.corr)
r2c.rand <- melt(r2c.rand, id.vars='subassem')


p.sub <- ggplot() + geom_line(data = c2r.rand, col = 'goldenrod', size = 0.2, alpha = 0.1,
                              aes(x = subassem, y = value, group = variable))
p.sub <- p.sub + geom_line(data = r2c.rand, col = 'darkolivegreen2', size = 0.2, alpha = 0.1,
                           aes(x = subassem, y =value, group = variable))
p.sub <- p.sub + geom_line(data = petcheyFD5m, col = 'gold1', size = 1.5,
                           aes(x = subassem, y = c2r.pearfd)) 
p.sub <- p.sub + geom_line(data = petcheyFD5m, col = 'forestgreen', size = 1.5,
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
p.info <- p.info + geom_line(data = petcheyFD5m, size = 1.5,
                             aes(x = sum.pq.c2r, y = c2r.pearfd, col = 'gold1')) + 
  geom_rug(data = petcheyFD5m, col = 'gold1', aes(x = sum.pq.c2r))
p.info <- p.info + geom_line(data = petcheyFD5m, size = 1.5,
                             aes(x = sum.pq.r2c, y = r2c.pearfd, col = 'forestgreen')) + 
  geom_rug(data = petcheyFD5m, col = 'forestgreen', alpha = 0.7, aes(x = sum.pq.r2c))
p.info <- p.info + scale_color_manual(name = '', values = c('gold1' = 'gold1', 'forestgreen' = 'forestgreen'),
                                      labels = c('rare to common', 'common to rare'))
p.info <- p.info + xlab('Cumulative information') + ylab('Correlation') + 
  theme_classic() + theme(axis.text = element_text(size = 11))

