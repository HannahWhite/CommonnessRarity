######################################################################################
### Modelled relationship between proportion common species in sites and diversity ###
######################################################################################

### Hannah White 05.05.2021

### Conditional autoregressive models of diversity measures as a function of the 
### proportion of the community made up of the 25% most common species

### NEED TO SET SEED FOR EACH MODEL

set.seed(32)

library(vegan)
library(FD)
library(cluster)

library(CARBayes)
library(MCMCglmm)

#### Read in site x species occurences ####

load('SpeciesOccurrences.RData')
load('SpeciesPercentCover.RData')

## area of occupancy function
aoo.func <- function(x){length(which(x > 0))} # calculates number of squares cover is > 0

## Calculate occurrences
occ4cm <- apply(site.sp4cm[, 5:62], 2, FUN = aoo.func)

occ20cm <- apply(site.sp20cm[,5:90], 2, FUN = aoo.func)

occ1m <- apply(site.sp1m[,5:119], 2, FUN = aoo.func)

occ5m <- apply(site.sp5m[,2:116], 2, FUN = aoo.func)

## Get orders
ord.c2r.4cm <- occ4cm[rev(order(occ4cm, sample(58, 58)))] # sample argument randomises sites
c2r4cm <- names(ord.c2r.4cm)
r2c4cm <- rev(c2r4cm)

ord.c2r.20cm <- occ20cm[rev(order(occ20cm, sample(86, 86)))] # sample argument randomises sites
c2r20cm <- names(ord.c2r.20cm)
r2c20cm <- rev(c2r20cm)

ord.c2r.1m <- occ1m[rev(order(occ1m, sample(115, 115)))] # sample argument randomises sites
c2r1m <- names(ord.c2r.1m)
r2c1m <- rev(c2r1m)

ord.c2r.5m <- occ5m[rev(order(occ5m, sample(115, 115)))] # sample argument randomises sites
c2r5m <- names(ord.c2r.5m)
r2c5m <- rev(c2r5m)

#### Top 25 % of each
top4cm <- c2r4cm[1:(ceiling(length(c2r4cm)/4))]
top20cm <- c2r20cm[1:(ceiling(length(c2r20cm)/4))]
top1m <- c2r1m[1:(ceiling(length(c2r1m)/4))]
top5m <- c2r5m[1:(ceiling(length(c2r5m)/4))]

common.prop <- function(x, y){
  sp.pres <- x[which(x==1)]
  common <- length(which(names(sp.pres) %in% y))
  rich <- sum(x)
  prop <- common/rich
  return(prop)
}

### Calculate proportion of each site that are 'common' species
common4cm <- apply(site.sp4cm[,5:62], 1, FUN = function(x) common.prop(x, top4cm))
common20cm <- apply(site.sp20cm[,5:90], 1, FUN = function(x) common.prop(x, top20cm))
common1m <- apply(site.sp1m[,5:119], 1, FUN = function(x) common.prop(x, top1m))
common5m <- apply(site.sp5m[,2:116], 1, FUN = function(x) common.prop(x, top5m))

## load traits
## load in traits
traits <- read.csv('FullTraits.csv', header = TRUE)
traits <- traits[, c('species', 'Canopy.height', 'Seed.mass', 'Leaf.size',  'SLA', 
                     'Insects', 'Selfing', 'Wind',
                     'Seed', 'Vegetative')]

traits$species <- gsub(' ', '.', traits$species)
traits$species <- gsub('-', '.', traits$species)

# change categorical traits to binary factor
traits[,6:10] <- ifelse(traits[,6:10] > 0, 1, 0)

### Calculate FD for each site
rownames(traits) <- traits[,1]
traits <- traits[,-1]

## set weights so that binary pollen vector fuzzy coding is equal to continuous traits
w <- c(1, 1, 1, 1, 1/3, 1/3, 1/3, 1/2, 1/2)

## create trait matrix for each spatial scale
traits4cm <- traits[rownames(traits) %in% names(site.sp4cm[,5:62]),]
traits20cm <- traits[rownames(traits) %in% names(site.sp20cm[,5:90]),]
traits1m <- traits[rownames(traits) %in% names(site.sp1m[,5:119]),]
traits5m <- traits[rownames(traits) %in% names(site.sp5m[,2:116]),]

##### Species Richness
sr4cm <- rowSums(site.sp4cm[,5:62])
sr20cm <- rowSums(site.sp20cm[,5:90])
sr1m <- rowSums(site.sp1m[,5:119])
sr5m <- rowSums(site.sp5m[,2:116])

##### P+G FD
trait.dist4cm <- daisy(traits4cm, metric='gower', weights = w) 
ftree4cm <- hclust(trait.dist4cm, method='average') #UPGMA clustering

full.fd4cm <- treedive(site.sp4cm[,5:62], ftree4cm) 

trait.dist20cm <- daisy(traits20cm, metric='gower', weights = w) 
ftree20cm <- hclust(trait.dist20cm, method='average') #UPGMA clustering

full.fd20cm <- treedive(site.sp20cm[,5:90], ftree20cm) 

trait.dist1m <- daisy(traits1m, metric='gower', weights = w) 
ftree1m <- hclust(trait.dist1m, method='average') #UPGMA clustering

full.fd1m <- treedive(site.sp1m[,5:119], ftree1m) 

trait.dist5m <- daisy(traits5m, metric='gower', weights = w) 
ftree5m <- hclust(trait.dist5m, method='average') #UPGMA clustering

full.fd5m <- treedive(site.sp5m[,2:116], ftree5m) 

##### Functional Dispersion and RaoQ
fd4cm <- dbFD(x=traits4cm, site.sp4cm[,5:62], w = w, calc.CWM = FALSE, calc.FRic = FALSE)
fdis4cm <- fd4cm$FDis
raoq4cm <- fd4cm$RaoQ

fd20cm <- dbFD(x=traits20cm, site.sp20cm[,5:90], w = w, calc.CWM = FALSE, calc.FRic = FALSE)
fdis20cm <- fd20cm$FDis
raoq20cm <- fd20cm$RaoQ

fd1m <- dbFD(x=traits1m, site.sp1m[,5:119], w = w, calc.CWM = FALSE, calc.FRic = FALSE)
fdis1m <- fd1m$FDis
raoq1m <- fd1m$RaoQ

fd5m <- dbFD(x=traits5m, site.sp5m[,2:116], w = w, calc.CWM = FALSE, calc.FRic = FALSE)
fdis5m <- fd5m$FDis
raoq5m <- fd5m$RaoQ

### Abundance weighted FDis and RaoQ
fdabund4cm <- dbFD(x=traits4cm, site.pc4cm[,5:62], w = w, calc.CWM = FALSE, calc.FRic = FALSE)
fdisabund4cm <- fdabund4cm$FDis
raoqabund4cm <- fdabund4cm$RaoQ

fdabund20cm <- dbFD(x=traits20cm, site.pc20cm[,5:90], w = w, calc.CWM = FALSE, calc.FRic = FALSE)
fdisabund20cm <- fdabund20cm$FDis
raoqabund20cm <- fdabund20cm$RaoQ

fdabund1m <- dbFD(x=traits1m, site.pc1m[,5:119], w = w, calc.CWM = FALSE, calc.FRic = FALSE)
fdisabund1m <- fdabund1m$FDis
raoqabund1m <- fdabund1m$RaoQ

fdabund5m <- dbFD(x=traits5m, site.pc5m[,2:116], w = w, calc.CWM = FALSE, calc.FRic = FALSE)
fdisabund5m <- fdabund5m$FDis
raoqabund5m <- fdabund5m$RaoQ

####### Models and Plots of Diversity~Proportion of common
commondiv4cm <- data.frame(plot.label = site.sp4cm$plot.label, common4cm, sr4cm, full.fd4cm, fdis4cm, raoq4cm, fdisabund4cm, raoqabund4cm)
commondiv20cm <- data.frame(plot.label = site.sp20cm$plot.label, common20cm, sr20cm, full.fd20cm, fdis20cm, raoq20cm, fdisabund20cm, raoqabund20cm)
commondiv1m <- data.frame(plot.label = site.sp1m$plot.label, common1m, sr1m, full.fd1m, fdis1m, raoq1m, fdisabund1m, raoqabund1m)
commondiv5m <- data.frame(plot.label = site.sp5m$Plot.ID, common5m, sr5m, full.fd5m, fdis5m, raoq5m, fdisabund5m, raoqabund5m)

### Get coordinates
## read in coordinates and create plot.label to match that in species df
sites4cm <- read.csv('coords0.04m.csv', header = TRUE)
sites4cm$plot.label <- paste(sites4cm$Name, '2008_0.04', sites4cm$Refno., sites4cm$Refletter., sep = '_')
sites4cm <- sites4cm[,c(7, 5:6)] # plot.label, east and north

sites20cm <- read.csv('coords0.2m.csv', header = TRUE)
sites20cm$plot.label <- paste(sites20cm$Name, '2008_0.2', sites20cm$Refno., sites20cm$Refletter., sep = '_')
sites20cm <- sites20cm[,c(7, 5:6)] # plot.label, east and north

sites1m <- read.csv('coords1m.csv', header = TRUE)
sites1m$plot.label <- paste(sites1m$Name, '2008_1', sites1m$Refno., sites1m$Refletter., sep = '_')
sites1m <- sites1m[,c(7, 5:6)] # plot.label, east and north

sites5m <- read.csv('Spatial machair.csv', header = TRUE)
sites5m <- sites5m[,c(1, 4, 6)]
names(sites5m) <- c('plot.label', 'Eastings', 'Northings')

commondiv4cm <- merge(sites4cm, commondiv4cm, by = 'plot.label', all.y = TRUE)
commondiv20cm <- merge(sites20cm, commondiv20cm, by = 'plot.label', all.y = TRUE)
commondiv1m <- merge(sites1m, commondiv1m, by = 'plot.label', all.y = TRUE)
commondiv5m <- merge(sites5m, commondiv5m, by = 'plot.label', all.y = TRUE)


### Models
## SR
## neighbourhood matrix
coords4cm <- commondiv4cm[,2:3]
dists4cm <- as.matrix(dist(coords4cm, method = 'maximum'))
w4cm <- array(0, dim = c(475, 475))
w4cm[dists4cm <= 16 & dists4cm > 0] <- 1 # gets 5 x 5 grid

coords20cm <- commondiv20cm[,2:3]
dists20cm <- as.matrix(dist(coords20cm, method = 'maximum'))
w20cm <- array(0, dim = c(475, 475))
w20cm[dists20cm <= 16 & dists20cm > 0] <- 1 # gets 5 x 5 grid

coords1m <- commondiv1m[,2:3]
dists1m <- as.matrix(dist(coords1m, method = 'maximum'))
w1m <- array(0, dim = c(475, 475))
w1m[dists1m <= 16 & dists1m > 0] <- 1 # gets 5 x 5 grid

## assume no autocorrelation for 5m
set.seed(32)
mod.sr4cm <- S.CARleroux(log(sr4cm) ~ common4cm, family = 'gaussian', data = commondiv4cm, W = w4cm, 
                         burnin = 10000, n.sample = 50000, prior.tau2 = c(0.5, 0.0005),
                         rho = 1) # rho set to 1 fits an intrinsic CAR prior
hist(residuals(mod.sr4cm))
print(mod.sr4cm) 
plot(mod.sr4cm$samples$beta)

set.seed(32)
mod.sr20cm <- S.CARleroux(log(sr20cm) ~ common20cm, family = 'gaussian', data = commondiv20cm, W = w20cm, 
                          burnin = 10000, n.sample = 50000, prior.tau2 = c(0.5, 0.0005),
                          rho = 1) # rho set to 1 fits an intrinsic CAR prior
hist(residuals(mod.sr20cm))
print(mod.sr20cm)
plot(mod.sr20cm$samples$beta)

set.seed(32)
mod.sr1m <- S.CARleroux(log(sr1m) ~ common1m, family = 'gaussian', data = commondiv1m, W = w1m, 
                        burnin = 10000, n.sample = 50000, prior.tau2 = c(0.5, 0.0005),
                        rho = 1) # rho set to 1 fits an intrinsic CAR prior
hist(residuals(mod.sr1m))
print(mod.sr1m)
plot(mod.sr1m$samples$beta)  

set.seed(32)
mod.sr5m <- MCMCglmm(log(sr5m) ~ common5m, family = 'gaussian', data = commondiv5m,
                     burnin = 10000, nitt = 50000, thin =1)
plot(mod.sr5m)
summary(mod.sr5m)
geweke.diag(mod.sr5m$Sol[,2])
median(mod.sr5m$Sol[,2])

### P + G
set.seed(32)
mod.pg4cm <- S.CARleroux(full.fd4cm ~ common4cm, family = 'gaussian', data = commondiv4cm, W = w4cm, 
                         burnin = 10000, n.sample = 50000, prior.tau2 = c(0.5, 0.0005),
                         rho = 1) # rho set to 1 fits an intrinsic CAR prior
hist(residuals(mod.pg4cm))
print(mod.pg4cm)
plot(mod.pg4cm$samples$beta)

set.seed(32)
mod.pg20cm <- S.CARleroux(full.fd20cm ~ common20cm, family = 'gaussian', data = commondiv20cm, W = w20cm, 
                          burnin = 10000, n.sample = 50000, prior.tau2 = c(0.5, 0.0005),
                          rho = 1) # rho set to 1 fits an intrinsic CAR prior
hist(residuals(mod.pg20cm))
print(mod.pg20cm)
plot(mod.pg20cm$samples$beta)

set.seed(32)
mod.pg1m <- S.CARleroux(full.fd1m ~ common1m, family = 'gaussian', data = commondiv1m, W = w1m, 
                        burnin = 10000, n.sample = 50000, prior.tau2 = c(0.5, 0.0005),
                        rho = 1) # rho set to 1 fits an intrinsic CAR prior
hist(residuals(mod.pg1m))
print(mod.pg1m)
plot(mod.pg1m$samples$beta)

set.seed(32)
mod.pg5m <- MCMCglmm(full.fd5m ~ common5m, family = 'gaussian', data = commondiv5m,
                     burnin = 10000, nitt = 50000, thin = 1)
plot(mod.pg5m)
summary(mod.pg5m)
geweke.diag(mod.pg5m$Sol[,2])
median(mod.pg5m$Sol[,2]) 

#### Functional Dispersion
set.seed(32)
mod.fdis4cm <- S.CARleroux(fdis4cm ~ common4cm, family = 'gaussian', data = commondiv4cm, W = w4cm, 
                           burnin = 10000, n.sample = 50000, prior.tau2 = c(0.5, 0.0005),
                           rho = 1) # rho set to 1 fits an intrinsic CAR prior
hist(residuals(mod.fdis4cm))
print(mod.fdis4cm)
plot(mod.fdis4cm$samples$beta)

set.seed(32)
mod.fdis20cm <- S.CARleroux(fdis20cm ~ common20cm, family = 'gaussian', data = commondiv20cm, W = w20cm, 
                            burnin = 10000, n.sample = 50000, prior.tau2 = c(0.5, 0.0005),
                            rho = 1) # rho set to 1 fits an intrinsic CAR prior
hist(residuals(mod.fdis20cm))
print(mod.fdis20cm)
plot(mod.fdis20cm$samples$beta)

set.seed(32)
mod.fdis1m <- S.CARleroux(fdis1m ~ common1m, family = 'gaussian', data = commondiv1m, W = w1m, 
                          burnin = 10000, n.sample = 50000, prior.tau2 = c(0.5, 0.0005),
                          rho = 1) # rho set to 1 fits an intrinsic CAR prior
hist(residuals(mod.fdis1m))
print(mod.fdis1m)
plot(mod.fdis1m$samples$beta)

set.seed(32)
mod.fdis5m <- MCMCglmm(fdis5m ~ common5m, family = 'gaussian', data = commondiv5m,
                       burnin = 10000, nitt = 50000, thin = 1)
plot(mod.fdis5m)
summary(mod.fdis5m)
geweke.diag(mod.fdis5m$Sol[,2])
median(mod.fdis5m$Sol[,2])


### Abundance weighted functional dispersion
set.seed(32)
mod.fdisabund4cm <- S.CARleroux(fdisabund4cm ~ common4cm, family = 'gaussian', data = commondiv4cm, W = w4cm, 
                                burnin = 10000, n.sample = 50000, prior.tau2 = c(0.5, 0.0005),
                                rho = 1) # rho set to 1 fits an intrinsic CAR prior
hist(residuals(mod.fdisabund4cm))
print(mod.fdisabund4cm)
plot(mod.fdisabund4cm$samples$beta)

set.seed(32)
mod.fdisabund20cm <- S.CARleroux(fdisabund20cm ~ common20cm, family = 'gaussian', data = commondiv20cm, W = w20cm, 
                                 burnin = 10000, n.sample = 50000, prior.tau2 = c(0.5, 0.0005),
                                 rho = 1) # rho set to 1 fits an intrinsic CAR prior
hist(residuals(mod.fdisabund20cm))
print(mod.fdisabund20cm)
plot(mod.fdisabund20cm$samples$beta)

set.seed(32)
mod.fdisabund1m <- S.CARleroux(fdisabund1m ~ common1m, family = 'gaussian', data = commondiv1m, W = w1m, 
                               burnin = 10000, n.sample = 50000, prior.tau2 = c(0.5, 0.0005),
                               rho = 1) # rho set to 1 fits an intrinsic CAR prior
hist(residuals(mod.fdisabund1m))
print(mod.fdisabund1m)
plot(mod.fdisabund1m$samples$beta)

set.seed(32)
mod.fdisabund5m <- MCMCglmm(fdisabund5m ~ common5m, family = 'gaussian', data = commondiv5m,
                            burnin = 10000, nitt = 50000, thin = 1)
plot(mod.fdisabund5m)
summary(mod.fdisabund5m)
geweke.diag(mod.fdisabund5m$Sol[,2])
median(mod.fdisabund5m$Sol[,2])
