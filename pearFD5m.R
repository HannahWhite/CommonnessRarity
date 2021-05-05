##############################################################################
#### Aggregate data to 5 x 5 m  and calculate FD and pearson correlations ####
##############################################################################

## Hannah White 05.05.2021

## Pearson correlations of sequential sub-assemblages (common to rare, rare to common)
## a the 5 m x 5 m scale. This code is slightly different to pearFD.R as it 
## aggregates the data to 5 m x 5 m first and it is run only over 19 plots

library(FD)
library(ggplot2)

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

ord.c2r <- occ[rev(order(occ,sample(115, 115)))] # sample argument randomises sp with equal aoo
common2rare <- names(ord.c2r)
rare2common <- rev(common2rare)

### Calculate FD for each site
rownames(traits) <- names(sp.order)
traits <- traits[,-1]

# only use traits in site.sp
traits <- traits[rownames(traits) %in% names(site.sp[,2:116]),]

## set weights so that binary pollen vector fuzzy coding is equal to continuous traits
w <- c(1, 1, 1, 1, 1/3, 1/3, 1/3, 1/2, 1/2)

#func.div <- dbFD(x=traits, site.sp[,2:116], w = w, calc.CWM = FALSE, calc.FRic = FALSE)
#full.fric <- func.div$FRic 
#full.fdis <- func.div$FDis
#full.raoq <- func.div$RaoQ

#### Common to Rare FD

c2r <- common2rare

c2r.fdis <- array(data = NA, dim = c(19, 115))
rownames(c2r.fdis) <- site.sp$Plot.ID

c2r.raoq <- array(data = NA, dim = c(19, 115))
rownames(c2r.raoq) <- site.sp$Plot.ID

system.time(
  for (sp in 2:115) {
    sub.sp <- c2r[1:sp]
    sub.comm <- site.sp[ ,sub.sp] #create sub-assemblage
    rownames(sub.comm) <- site.sp$Plot.ID
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


## Pearson correlations
c2r.pear.fdis <- apply(c2r.fdis[,1:115], 2, FUN = function(x2) cor(cbind(x2, c2r.fdis[,115]), method = 'pearson', use = 'complete.obs')[[2]])
c2r.pear.raoq <- apply(c2r.raoq[,1:115], 2, FUN = function(x2) cor(cbind(x2, c2r.raoq[,115]), method = 'pearson', use = 'complete.obs')[[2]])


c2r.pearcorrs5m <- data.frame(c2r.pear.fdis, c2r.pear.raoq)

#### Rare to Common

r2c <- rare2common

r2c.fdis <- array(data = NA, dim = c(19, 115))
rownames(r2c.fdis) <- site.sp$Plot.ID

r2c.raoq <- array(data = NA, dim = c(19, 115))
rownames(r2c.raoq) <- site.sp$Plot.ID

system.time(
  for (sp in 2:115) {
    sub.sp <- r2c[1:sp]
    sub.comm <- site.sp[ ,sub.sp] #create sub-assemblage
    rownames(sub.comm) <- site.sp$Plot.ID
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


# Pearson correlations for rare to common
r2c.pear.fdis <- apply(r2c.fdis[,1:115], 2, FUN = function(x2) cor(cbind(x2, r2c.fdis[,115]), method = 'pearson', use = 'complete.obs')[[2]])
r2c.pear.raoq <- apply(r2c.raoq[,1:115], 2, FUN = function(x2) cor(cbind(x2, r2c.raoq[,115]), method = 'pearson', use = 'complete.obs')[[2]])

r2c.pearcorrs5m <- data.frame(r2c.pear.fdis, r2c.pear.raoq)

#### Subassemblage
subassem <- 1:115

pearcorrsFD.5m <- data.frame(subassem, c2r.pearcorrs5m, r2c.pearcorrs5m)


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


pearFD.df <- data.frame(pearcorrsFD.5m, sum.pq.c2r, sum.pq.r2c)

pearFD5m <- pearFD.df

#############
### Plots ###
#############


## load randomisations
load('randomcorrs.c2r5m.RData')
load('randomcorrs.r2c5m.RData')

subassem <- pearFD5m$subassem
sum.pq.c2r <- pearFD5m$sum.pq.c2r
sum.pq.r2c <- pearFD5m$sum.pq.r2c

### Functional dispersion

# Subassemblage 

c2r.fdis.rand <- data.frame(subassem, fdis.corr)
c2r.fdis.rand <- melt(c2r.fdis.rand, id.vars="subassem")

r2c.fdis.rand <- data.frame(subassem, r2c.fdis.corr)
r2c.fdis.rand <- melt(r2c.fdis.rand, id.vars='subassem')


p.subfdis <- ggplot() + geom_line(data = c2r.fdis.rand, col = 'goldenrod', size = 0.2, alpha = 0.1,
                                  aes(x = subassem, y = value, group = variable))
p.subfdis <- p.subfdis + geom_line(data = r2c.fdis.rand, col = 'darkolivegreen2', size = 0.2, alpha = 0.1,
                                   aes(x = subassem, y =value, group = variable))
p.subfdis <- p.subfdis + geom_line(data = pearFD5m, col = 'gold1', size = 1.5,
                                   aes(x = subassem, y = c2r.pear.fdis)) 
p.subfdis <- p.subfdis + geom_line(data = pearFD5m, col = 'forestgreen', size = 1.5,
                                   aes(x = subassem, y = r2c.pear.fdis))
p.subfdis <- p.subfdis + xlab('Subassemblage size') + ylab('Correlation') + 
  theme_classic() + theme(axis.text = element_text(size = 11))



## Cumulative information

c2r.fdis.rand.pq <- data.frame(sum.pq.c2r = pearFD5m$sum.pq.c2r, fdis.corr)
c2r.fdis.rand.pq <- melt(c2r.fdis.rand.pq, id.vars="sum.pq.c2r")

r2c.fdis.rand.pq <- data.frame(sum.pq.r2c = pearFD5m$sum.pq.r2c, r2c.fdis.corr)
r2c.fdis.rand.pq <- melt(r2c.fdis.rand.pq, id.vars='sum.pq.r2c')


p.infofdis <- ggplot() + geom_line(data = c2r.fdis.rand.pq, col = 'goldenrod', size = 0.2, alpha = 0.1,
                                   aes(x = sum.pq.c2r, y = value, group = variable))
p.infofdis <- p.infofdis + geom_line(data = r2c.fdis.rand.pq, col = 'darkolivegreen2', size = 0.2, alpha = 0.1,
                                     aes(x = sum.pq.r2c, y =value, group = variable))
p.infofdis <- p.infofdis + geom_line(data = pearFD5m, size = 1.5,
                                     aes(x = sum.pq.c2r, y = c2r.pear.fdis, col = 'gold1')) + 
  geom_rug(data = pearFD5m, col = 'gold1', aes(x = sum.pq.c2r))
p.infofdis <- p.infofdis + geom_line(data = pearFD5m, size = 1.5,
                                     aes(x = sum.pq.r2c, y = r2c.pear.fdis, col = 'forestgreen')) + 
  geom_rug(data = pearFD5m, col = 'forestgreen', alpha = 0.7, aes(x = sum.pq.r2c))
p.infofdis <- p.infofdis + scale_color_manual(name = '', values = c('gold1' = 'gold1', 'forestgreen' = 'forestgreen'),
                                              labels = c('rare to common', 'common to rare'))
p.infofdis <- p.infofdis + xlab('Cumulative information') + ylab('Correlation') + 
  theme_classic() + theme(axis.text = element_text(size = 11))



### Raoq
# Subassemblage 

c2r.raoq.rand <- data.frame(subassem, raoq.corr)
c2r.raoq.rand <- melt(c2r.raoq.rand, id.vars="subassem")

r2c.raoq.rand <- data.frame(subassem, r2c.raoq.corr)
r2c.raoq.rand <- melt(r2c.raoq.rand, id.vars='subassem')


p.subraoq <- ggplot() + geom_line(data = c2r.raoq.rand, col = 'goldenrod', size = 0.2, alpha = 0.1,
                                  aes(x = subassem, y = value, group = variable))
p.subraoq <- p.subraoq + geom_line(data = r2c.raoq.rand, col = 'darkolivegreen2', size = 0.2, alpha = 0.1,
                                   aes(x = subassem, y =value, group = variable))
p.subraoq <- p.subraoq + geom_line(data = pearFD5m, col = 'gold1', size = 1.5,
                                   aes(x = subassem, y = c2r.pear.raoq)) 
p.subraoq <- p.subraoq + geom_line(data = pearFD5m, col = 'forestgreen', size = 1.5,
                                   aes(x = subassem, y = r2c.pear.raoq))
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
p.inforaoq <- p.inforaoq + geom_line(data = pearFD5m, size = 1.5,
                                     aes(x = sum.pq.c2r, y = c2r.pear.raoq, col = 'gold1')) + 
  geom_rug(data = pearFD5m, col = 'gold1', aes(x = sum.pq.c2r))
p.inforaoq <- p.inforaoq + geom_line(data = pearFD5m, size = 1.5,
                                     aes(x = sum.pq.r2c, y = r2c.pear.raoq, col = 'forestgreen')) + 
  geom_rug(data = pearFD5m, col = 'forestgreen', alpha = 0.7, aes(x = sum.pq.r2c))
p.inforaoq <- p.inforaoq + scale_color_manual(name = '', values = c('gold1' = 'gold1', 'forestgreen' = 'forestgreen'),
                                              labels = c('rare to common', 'common to rare'))
p.inforaoq <- p.inforaoq + xlab('Cumulative information') + ylab('Correlation') + 
  theme_classic() + theme(axis.text = element_text(size = 11))
