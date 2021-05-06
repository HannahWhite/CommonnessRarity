###########################################
### Functional diversity randomisations ###
###########################################

### Hannah White 06.05.2021

### Uses name swap algorithm on species x trait matrix and then calculates FDis and RaoQ based
### on these randomised matrices. Sequential correlations are then calculated.
### This code needs to be adjusted depending on scale under investigation.

### Name swap algorithm

set.seed(32)

library(FD)

# load and sort out data
# get traits data
traits <- read.csv('FullTraits.csv', header = TRUE)
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
specs <- species.scale[,13:150]
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
ord.c2r <- occ[rev(order(occ,sample(n, n)))] # sample argument randomises sites
common2rare <- names(ord.c2r)
rare2common <- rev(common2rare)

### Calculate FD for each site
rownames(traits) <- names(sp.order)
traits <- traits[,-1]

# only use traits in site.sp
traits <- traits[rownames(traits) %in% names(site.sp[,5:xdim]),]

### Full FD
## set weights so that binary pollen vector fuzzy coding is equal to continuous traits
w <- c(1, 1, 1, 1, 1/3, 1/3, 1/3, 1/2, 1/2)

############## Randomisation ##########################
# randomise species names in trait matrix
rand.name <- replicate(500, sample(rownames(traits), dim(traits)[1], replace=FALSE))

##initiate empty list
rand.traits <- list()

###create dataframes of random trait species

for(i in 1:500){
  name <- paste('mat', i, sep='_')
  tmp <- data.frame(traits, row.names=rand.name[,i])
  rand.traits[[name]] <- tmp
}

#### Common to Rare FD
c2r <- common2rare

c2r.fdis <- array(data = 0, dim = c(475, n, 500))
rownames(c2r.fdis) <- site.sp$plot.label

c2r.raoq <- array(data = 0, dim = c(475, n, 500))
rownames(c2r.raoq) <- site.sp$plot.label
system.time(
  for (j in 1:500){ # 
    traits.tmp <- rand.traits[[j]]
    for (sp in 2:n) {
      sub.sp <- c2r[1:sp]
      sub.comm <- site.sp[ ,sub.sp] #create sub-assemblage
      rownames(sub.comm) <- site.sp$plot.label
      
      sub.traits <- traits.tmp[rownames(traits.tmp) %in% sub.sp, ] #create smaller trait matrix with species just in sub assemblage
      # change sub.traits to alphabetical order
      traits.alpha <- sub.traits[order(rownames(sub.traits)),]
      
      if(is.na(gowdis(traits.alpha, w=w))){
        out.fric <- NA
        #out.fdis <- NA 
        #out.raoq <- NA
      }
      else {
        # can only calculate FD for sites with species present
        pop.sites <- sub.comm[rowSums(sub.comm)>0,]
        # change sub.comm columns to alphabetical order
        comm.alpha <- pop.sites[, order(colnames(pop.sites))] 
        
        # calculate functional diversity
        out <- FD::dbFD(traits.alpha, comm.alpha, calc.CWM = FALSE, calc.FDiv = FALSE, w = w)
        out.fric <- out$FRic
        out.fdis <- out$FDis
        out.raoq <- out$RaoQ }
      
      c2r.fric[names(out.fric), sp, j] <- out.fric
      c2r.fdis[names(out.fdis), sp, j] <- out.fdis
      c2r.raoq[names(out.raoq), sp, j] <- out.raoq
      
    }
  }
) 

## Change 0s to NAs where no species present

c2r.fdis[,1,] <- ifelse(site.sp[,c2r[1]] > 0, 0, NA)
for (r in 1:500){
  for (s in 2:n){
    abs <- which(rowSums(site.sp[,c2r[1:s]]) == 0)
    c2r.fdis[abs, s, r] <- NA
    
  }
}

c2r.raoq[,1,] <- ifelse(site.sp[,c2r[1]] > 0, 0, NA)
for (r in 1:500){
  for (s in 2:n){
    abs <- which(rowSums(site.sp[,c2r[1:s]]) == 0)
    c2r.raoq[abs, s, r] <- NA
    
  }
}


#### Correlations

fdis.corr <- array(data = NA, dim = c(n, 500))

for(rr in 1:500){
  fdis.tmp <- c2r.fdis[,,rr]
  fdis.corr[,rr] <- apply(fdis.tmp[,1:n], 2, FUN = function(x2) cor(cbind(x2, fdis.tmp[,n]), method = 'pearson', use = 'complete.obs')[[2]])
  print(rr)
}

raoq.corr <- array(data = NA, dim = c(n, 500))

for(rr in 1:500){
  raoq.tmp <- c2r.raoq[,,rr]
  raoq.corr[,rr] <- apply(raoq.tmp[,1:n], 2, FUN = function(x2) cor(cbind(x2, raoq.tmp[,n]), method = 'pearson', use = 'complete.obs')[[2]])
  print(rr)
}


###### Rare to common FD
r2c <- rare2common

r2c.fdis <- array(data = 0, dim = c(475, n, 500))
rownames(r2c.fdis) <- site.sp$plot.label

r2c.raoq <- array(data = 0, dim = c(475, n, 500))
rownames(r2c.raoq) <- site.sp$plot.label
system.time( ### 
  for (j in 1:500){
    traits.tmp <- rand.traits[[j]]
    for (sp in 2:n) {
      sub.sp <- r2c[1:sp]
      sub.comm <- site.sp[ ,sub.sp] #create sub-assemblage
      rownames(sub.comm) <- site.sp$plot.label
      
      sub.traits <- traits.tmp[rownames(traits.tmp) %in% sub.sp, ] #create smaller trait matrix with species just in sub assemblage
      # change sub.traits to alphabetical order
      traits.alpha <- sub.traits[order(rownames(sub.traits)),]
      
      if(is.na(gowdis(traits.alpha, w=w))){
        
        out.fdis <- NA 
        out.raoq <- NA
      }
      else{  
        # can only calculate FD for sites with species present
        pop.sites <- sub.comm[rowSums(sub.comm)>0,]
        # change sub.comm columns to alphabetical order
        comm.alpha <- pop.sites[, order(colnames(pop.sites))] 
        
        # calculate functional diversity
        out <- FD::dbFD(traits.alpha, comm.alpha, calc.CWM = FALSE, calc.FDiv = FALSE, calc.FRic = FALSE, w = w)
        
        out.fdis <- out$FDis
        out.raoq <- out$RaoQ }
      
      
      
      r2c.fdis[names(out.fdis), sp, j] <- out.fdis
      r2c.raoq[names(out.raoq), sp, j] <- out.raoq
    } 
    
  }
) 

## Replace 0s with NAs where there are no species present

r2c.fdis[,1,] <- ifelse(site.sp[,r2c[1]] > 0, 0, NA)
for (r in 1:500){
  for (s in 2:n){
    abs <- which(rowSums(site.sp[,r2c[1:s]]) == 0)
    r2c.fdis[abs, s, r] <- NA
    
  }
}

r2c.raoq[,1,] <- ifelse(site.sp[,r2c[1]] > 0, 0, NA)
for (r in 1:500){
  for (s in 2:n){
    abs <- which(rowSums(site.sp[,r2c[1:s]]) == 0)
    r2c.raoq[abs, s, r] <- NA
    
  }
}


### Correlations

r2c.fdis.corr <- array(data = NA, dim = c(n, 500))

for(rr in 1:500){
  fdis.tmp <- r2c.fdis[,,rr]
  r2c.fdis.corr[,rr] <- apply(fdis.tmp[,1:n], 2, FUN = function(x2) cor(cbind(x2, fdis.tmp[,n]), method = 'pearson', use = 'complete.obs')[[2]])
}


r2c.raoq.corr <- array(data = NA, dim = c(n, 500))

for(rr in 1:500){
  raoq.tmp <- r2c.raoq[,,rr]
  r2c.raoq.corr[,rr] <- apply(raoq.tmp[,1:n], 2, FUN = function(x2) cor(cbind(x2, raoq.tmp[,115]), method = 'pearson', use = 'complete.obs')[[2]])
}

### save()

