############################################################################
### Combined plot of sequential correlations at different spatial scales ###
############################################################################

## This code can be adapted to create the plot for any functional diversity measure - 
## it is currently set up for functional dispersion

### Hannah White 22.04.2021
### Edited 10.06.2021 to include confidence intervals
### Edited 17.12.2021 to tidy up legend

library(ggplot2)
library(reshape2)
library(cowplot)

### 4cm
load('pearFD4cm.RData')

# read in correlations of randomisations - generated from randomisation R code files
load('randomcorrs.c2r4cm.RData')
load('randomcorrs.r2c4cm.RData')

## get 95% quantiles

upper.c2r <- apply(fdis.corr, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
lower.c2r <- apply(fdis.corr, 1, function(x) quantile(x, 0.025, na.rm = TRUE))

upper.r2c <- apply(r2c.fdis.corr, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
upper.r2c[2] <- NA
lower.r2c <- apply(r2c.fdis.corr, 1, function(x) quantile(x, 0.025, na.rm = TRUE))
lower.r2c[2] <- NA

subassem <- pearFD4cm$subassem
sum.pq.c2r <- pearFD4cm$sum.pq.c2r
sum.pq.r2c <- pearFD4cm$sum.pq.r2c

c2r.fdis.rand <- data.frame(subassem, fdis.corr)
c2r.fdis.rand <- melt(c2r.fdis.rand, id.vars="subassem")

r2c.fdis.rand <- data.frame(subassem, r2c.fdis.corr)
r2c.fdis.rand <- melt(r2c.fdis.rand, id.vars='subassem')

c2r.fdis.rand.pq <- data.frame(sum.pq.c2r, fdis.corr)
c2r.fdis.rand.pq <- melt(c2r.fdis.rand.pq, id.vars="sum.pq.c2r")

r2c.fdis.rand.pq <- data.frame(sum.pq.r2c, r2c.fdis.corr)
r2c.fdis.rand.pq <- melt(r2c.fdis.rand.pq, id.vars='sum.pq.r2c')

## create dataframe for confidence intervals
upper.c2r <- data.frame(sum.pq.c2r, upper.c2r)
lower.c2r <- data.frame(sum.pq.c2r, lower.c2r)

upper.r2c <- data.frame(sum.pq.r2c, upper.r2c)
lower.r2c <- data.frame(sum.pq.r2c, lower.r2c)

p.infofdis4cm <- ggplot() + geom_line(data = c2r.fdis.rand.pq, col = 'gold3', size = 0.2, alpha = 0.1,
                                      aes(x = sum.pq.c2r, y = value, group = variable))
p.infofdis4cm <- p.infofdis4cm + geom_line(data = r2c.fdis.rand.pq, col = 'royalblue2', size = 0.2, alpha = 0.03,
                                           aes(x = sum.pq.r2c, y =value, group = variable))

p.infofdis4cm <- p.infofdis4cm + geom_line(data = upper.c2r, size = 0.8, col = 'gold3', linetype = 'dashed',
                                           aes(x = sum.pq.c2r, y = upper.c2r)) # add cI
p.infofdis4cm <- p.infofdis4cm + geom_line(data = lower.c2r, size = 0.8, col = 'gold3', linetype = 'dashed',
                                           aes(x = sum.pq.c2r, y = lower.c2r)) # add CI
p.infofdis4cm <- p.infofdis4cm + geom_line(data = upper.r2c, size = 0.8, col = 'royalblue1', linetype = 'dashed',
                                           aes(x = sum.pq.r2c, y = upper.r2c))
p.infofdis4cm <- p.infofdis4cm + geom_line(data = lower.r2c, size = 0.8, linetype = 'dashed',
                                           aes(x = sum.pq.r2c, y = lower.r2c, col = 'royalblue1'))

p.infofdis4cm <- p.infofdis4cm + geom_line(data = pearFD4cm, size = 1.5,
                                           aes(x = sum.pq.c2r, y = c2r.fdis, col = 'gold1')) + 
  geom_rug(data = pearFD4cm, col = 'gold1', aes(x = sum.pq.c2r))
p.infofdis4cm <- p.infofdis4cm + geom_line(data = pearFD4cm, size = 1.5,
                                           aes(x = sum.pq.r2c, y = r2c.fdis, col = 'royalblue4')) + 
  geom_rug(data = pearFD4cm, col = 'royalblue4', alpha = 0.7, aes(x = sum.pq.r2c))
p.infofdis4cm <- p.infofdis4cm + scale_color_manual(name = '', 
                                                    values = c('gold1' = 'gold1', 'royalblue4' = 'royalblue4', 'royalblue1' = 'royalblue1'),
                                                    labels = c('95% CI', 'rare to\ncommon', 'common to\nrare'),
                                                    guide = guide_legend(reverse = TRUE,
                                                                         override.aes = list(
                                                                           linetype = c('solid', 'solid', 'dashed'),
                                                                           colour = c('gold1', 'royalblue4', 'black'),
                                                                           size = c(1.5, 1.5, 0.8))))

p.infofdis4cm <- p.infofdis4cm + ggtitle('0.04 m x 0.04 m')
p.infofdis4cm <- p.infofdis4cm + xlab('Cumulative information') + ylab('Correlation') + 
  theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
                          plot.title = element_text(hjust = 0.5, size = 21))

rm(list=setdiff(ls(), "p.infofdis4cm"))


### 20cm

load('pearFD20cm.RData')

# read in correlations of randomisations
load('randomcorrs.c2r20cm.RData')
load('randomcorrs.r2c20cm.RData')

# get 95% quantiles
upper.c2r <- apply(fdis.corr, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
lower.c2r <- apply(fdis.corr, 1, function(x) quantile(x, 0.025, na.rm = TRUE))

upper.r2c <- apply(r2c.fdis.corr, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
upper.r2c[2] <- NA
lower.r2c <- apply(r2c.fdis.corr, 1, function(x) quantile(x, 0.025, na.rm = TRUE))
lower.r2c[2] <- NA

subassem <- pearFD20cm$subassem
sum.pq.c2r <- pearFD20cm$sum.pq.c2r
sum.pq.r2c <- pearFD20cm$sum.pq.r2c

c2r.fdis.rand.pq <- data.frame(sum.pq.c2r, fdis.corr)
c2r.fdis.rand.pq <- melt(c2r.fdis.rand.pq, id.vars="sum.pq.c2r")

r2c.fdis.rand.pq <- data.frame(sum.pq.r2c, r2c.fdis.corr)
r2c.fdis.rand.pq <- melt(r2c.fdis.rand.pq, id.vars='sum.pq.r2c')

## create dataframe for confidence intervals
upper.c2r <- data.frame(sum.pq.c2r, upper.c2r)
lower.c2r <- data.frame(sum.pq.c2r, lower.c2r)

upper.r2c <- data.frame(sum.pq.r2c, upper.r2c)
lower.r2c <- data.frame(sum.pq.r2c, lower.r2c)


p.infofdis20cm <- ggplot() + geom_line(data = c2r.fdis.rand.pq, col = 'gold3', size = 0.2, alpha = 0.1,
                                       aes(x = sum.pq.c2r, y = value, group = variable))
p.infofdis20cm <- p.infofdis20cm + geom_line(data = r2c.fdis.rand.pq, col = 'royalblue2', size = 0.2, alpha = 0.03,
                                             aes(x = sum.pq.r2c, y =value, group = variable))

p.infofdis20cm <- p.infofdis20cm + geom_line(data = upper.c2r, size = 0.8,col = 'gold3', linetype = 'dashed',
                                             aes(x = sum.pq.c2r, y = upper.c2r)) # add cI
p.infofdis20cm <- p.infofdis20cm + geom_line(data = lower.c2r, size = 0.8, col = 'gold3', linetype = 'dashed',
                                             aes(x = sum.pq.c2r, y = lower.c2r)) # add CI
p.infofdis20cm <- p.infofdis20cm + geom_line(data = upper.r2c, size = 0.8, col = 'royalblue1', linetype = 'dashed',
                                             aes(x = sum.pq.r2c, y = upper.r2c))
p.infofdis20cm <- p.infofdis20cm + geom_line(data = lower.r2c, size = 0.8, linetype = 'dashed',
                                             aes(x = sum.pq.r2c, y = lower.r2c, col = 'royalblue1'))


p.infofdis20cm <- p.infofdis20cm + geom_line(data = pearFD20cm, size = 1.5,
                                             aes(x = sum.pq.c2r, y = c2r.fdis, col = 'gold1')) + 
  geom_rug(data = pearFD20cm, col = 'gold1', aes(x = sum.pq.c2r))
p.infofdis20cm <- p.infofdis20cm + geom_line(data = pearFD20cm, size = 1.5,
                                             aes(x = sum.pq.r2c, y = r2c.fdis, col = 'royalblue4')) + 
  geom_rug(data = pearFD20cm, col = 'royalblue4', alpha = 0.7, aes(x = sum.pq.r2c))
p.infofdis20cm <- p.infofdis20cm + scale_color_manual(name = '', 
                                                      values = c('gold1' = 'gold1', 'royalblue4' = 'royalblue4', 'royalblue1' = 'royalblue1'),
                                                      labels = c('95% CI', 'rare to\ncommon', 'common to\nrare'),
                                                      guide = guide_legend(reverse = TRUE,
                                                                           override.aes = list(
                                                                             linetype = c('solid', 'solid', 'dashed'),
                                                                             colour = c('gold1', 'royalblue4', 'black'),
                                                                             size = c(1.5, 1.5, 0.8))))

p.infofdis20cm <- p.infofdis20cm + ggtitle('0.2 m x 0.2 m')
p.infofdis20cm <- p.infofdis20cm + xlab('Cumulative information') + ylab(NULL) + 
  theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
                          plot.title = element_text(hjust = 0.5, size = 21))

rm(list=setdiff(ls(), c("p.infofdis4cm", "p.infofdis20cm")))


### 1m
load('pearFD1m.RData')

## load randomisations
load('randomcorrs.c2r.RData')
load('randomcorrs.r2c.RData')

# get 95% quantiles
upper.c2r <- apply(fdis.corr, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
lower.c2r <- apply(fdis.corr, 1, function(x) quantile(x, 0.025, na.rm = TRUE))

upper.r2c <- apply(r2c.fdis.corr, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
lower.r2c <- apply(r2c.fdis.corr, 1, function(x) quantile(x, 0.025, na.rm = TRUE))

subassem <- pearFD1m$subassem
sum.pq.c2r <- pearFD1m$sum.pq.c2r
sum.pq.r2c <- pearFD1m$sum.pq.r2c

c2r.fdis.rand.pq <- data.frame(sum.pq.c2r = pearFD1m$sum.pq.c2r, fdis.corr)
c2r.fdis.rand.pq <- melt(c2r.fdis.rand.pq, id.vars="sum.pq.c2r")

r2c.fdis.rand.pq <- data.frame(sum.pq.r2c = pearFD1m$sum.pq.r2c, r2c.fdis.corr)
r2c.fdis.rand.pq <- melt(r2c.fdis.rand.pq, id.vars='sum.pq.r2c')

## create dataframe for confidence intervals
upper.c2r <- data.frame(sum.pq.c2r, upper.c2r)
lower.c2r <- data.frame(sum.pq.c2r, lower.c2r)

upper.r2c <- data.frame(sum.pq.r2c, upper.r2c)
lower.r2c <- data.frame(sum.pq.r2c, lower.r2c)


p.infofdis1m <- ggplot() + geom_line(data = c2r.fdis.rand.pq, col = 'gold3', size = 0.2, alpha = 0.1,
                                     aes(x = sum.pq.c2r, y = value, group = variable))
p.infofdis1m <- p.infofdis1m + geom_line(data = r2c.fdis.rand.pq, col = 'royalblue2', size = 0.2, alpha = 0.03,
                                         aes(x = sum.pq.r2c, y =value, group = variable))

p.infofdis1m <- p.infofdis1m + geom_line(data = upper.c2r, size = 0.8,col = 'gold3', linetype = 'dashed',
                                         aes(x = sum.pq.c2r, y = upper.c2r)) # add cI
p.infofdis1m <- p.infofdis1m + geom_line(data = lower.c2r, size = 0.8, col = 'gold3', linetype = 'dashed',
                                         aes(x = sum.pq.c2r, y = lower.c2r)) # add CI
p.infofdis1m <- p.infofdis1m + geom_line(data = upper.r2c, size = 0.8, col = 'royalblue1', linetype = 'dashed',
                                         aes(x = sum.pq.r2c, y = upper.r2c))
p.infofdis1m <- p.infofdis1m + geom_line(data = lower.r2c, size = 0.8, linetype = 'dashed',
                                         aes(x = sum.pq.r2c, y = lower.r2c, col = 'royalblue1'))

p.infofdis1m <- p.infofdis1m + geom_line(data = pearFD1m, size = 1.5,
                                         aes(x = sum.pq.c2r, y = c2r.fdis, col = 'gold1')) + 
  geom_rug(data = pearFD1m, col = 'gold1', aes(x = sum.pq.c2r))
p.infofdis1m <- p.infofdis1m + geom_line(data = pearFD1m, size = 1.5,
                                         aes(x = sum.pq.r2c, y = r2c.fdis, col = 'royalblue4')) + 
  geom_rug(data = pearFD1m, col = 'royalblue4', alpha = 0.7, aes(x = sum.pq.r2c))
p.infofdis1m <- p.infofdis1m + scale_color_manual(name = '', 
                                                  values = c('gold1' = 'gold1', 'royalblue4' = 'royalblue4', 'royalblue1' = 'royalblue1'),
                                                  labels = c('95% CI', 'rare to\ncommon', 'common to\nrare'),
                                                  guide = guide_legend(reverse = TRUE,
                                                                       override.aes = list(
                                                                         linetype = c('solid', 'solid', 'dashed'),
                                                                         colour = c('gold1', 'royalblue4', 'black'),
                                                                         size = c(1.5, 1.5, 0.8))))

p.infofdis1m <- p.infofdis1m + ggtitle('1 m x 1 m')
p.infofdis1m <- p.infofdis1m + xlab('Cumulative information') + ylab(NULL) + 
  theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
                          plot.title = element_text(hjust = 0.5, size = 21))

rm(list=setdiff(ls(), c("p.infofdis4cm", "p.infofdis20cm", "p.infofdis1m")))

### 5m

load('pearFD5m.RData')

## load randomisations
load('randomcorrs.c2r5m.RData')
load('randomcorrs.r2c5m.RData')

# get 95% quantiles
upper.c2r <- apply(fdis.corr, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
lower.c2r <- apply(fdis.corr, 1, function(x) quantile(x, 0.025, na.rm = TRUE))

upper.r2c <- apply(r2c.fdis.corr, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
lower.r2c <- apply(r2c.fdis.corr, 1, function(x) quantile(x, 0.025, na.rm = TRUE))

subassem <- pearFD5m$subassem
sum.pq.c2r <- pearFD5m$sum.pq.c2r
sum.pq.r2c <- pearFD5m$sum.pq.r2c

c2r.fdis.rand.pq <- data.frame(sum.pq.c2r = pearFD5m$sum.pq.c2r, fdis.corr)
c2r.fdis.rand.pq <- melt(c2r.fdis.rand.pq, id.vars="sum.pq.c2r")

r2c.fdis.rand.pq <- data.frame(sum.pq.r2c = pearFD5m$sum.pq.r2c, r2c.fdis.corr)
r2c.fdis.rand.pq <- melt(r2c.fdis.rand.pq, id.vars='sum.pq.r2c')

## create dataframe for confidence intervals
upper.c2r <- data.frame(sum.pq.c2r, upper.c2r)
lower.c2r <- data.frame(sum.pq.c2r, lower.c2r)

upper.r2c <- data.frame(sum.pq.r2c, upper.r2c)
lower.r2c <- data.frame(sum.pq.r2c, lower.r2c)

p.infofdis5m <- ggplot() + geom_line(data = c2r.fdis.rand.pq, size = 0.2, alpha = 0.1, col = 'gold3',
                                     aes(x = sum.pq.c2r, y = value, group = variable))
p.infofdis5m <- p.infofdis5m + geom_line(data = r2c.fdis.rand.pq, size = 0.2, alpha = 0.03, col = 'royalblue2',
                                         aes(x = sum.pq.r2c, y =value, group = variable))

p.infofdis5m <- p.infofdis5m + geom_line(data = upper.c2r, size = 0.8,col = 'gold3', linetype = 'dashed',
                                         aes(x = sum.pq.c2r, y = upper.c2r)) # add cI
p.infofdis5m <- p.infofdis5m + geom_line(data = lower.c2r, size = 0.8, col = 'gold3', linetype = 'dashed',
                                         aes(x = sum.pq.c2r, y = lower.c2r)) # add CI
p.infofdis5m <- p.infofdis5m + geom_line(data = upper.r2c, size = 0.8, col = 'royalblue1', linetype = 'dashed',
                                         aes(x = sum.pq.r2c, y = upper.r2c))
p.infofdis5m <- p.infofdis5m + geom_line(data = lower.r2c, size = 0.8, linetype = 'dashed', col = 'royalblue1',
                                         aes(x = sum.pq.r2c, y = lower.r2c))

p.infofdis5m <- p.infofdis5m + geom_line(data = pearFD5m, size = 1.5,
                                         aes(x = sum.pq.c2r, y = c2r.pear.fdis, col = 'gold1')) + 
  geom_rug(data = pearFD5m, col = 'gold1', aes(x = sum.pq.c2r))
p.infofdis5m <- p.infofdis5m + geom_line(data = pearFD5m, size = 1.5,
                                         aes(x = sum.pq.r2c, y = r2c.pear.fdis, col = 'royalblue4')) + 
  geom_rug(data = pearFD5m, col = 'royalblue4', alpha = 0.7, aes(x = sum.pq.r2c))
p.infofdis5m <- p.infofdis5m + scale_color_manual(name = '', 
                                                  values = c('gold1' = 'gold1', 'royalblue4' = 'royalblue4'),
                                                  labels = c( 'common to\nrare', 'rare to\ncommon'),
                                                  guide = guide_legend(override.aes = list(
                                                    linetype = c('solid', 'solid'),
                                                    colour = c('gold1', 'royalblue4'),
                                                    size = c(1.5, 1.5),
                                                    alpha = c(1, 1))))

p.infofdis5m <- p.infofdis5m + ggtitle('5 m x 5 m')
p.infofdis5m <- p.infofdis5m + xlab('Cumulative information') + ylab(NULL) + 
  theme_classic() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
                          plot.title = element_text(hjust = 0.5, size = 21),
                          legend.position = 'bottom',
                          legend.key.width = unit(0.5, 'in'))

rm(list=setdiff(ls(), c("p.infofdis4cm", "p.infofdis20cm", "p.infofdis1m", "p.infofdis5m")))

#### Combine plots

legend <- get_legend(
  p.infofdis5m + theme(legend.text = element_text(size = 17), legend.box.margin = margin(0, 0, 0, 0))
)

p.row <- plot_grid(p.infofdis4cm + theme(legend.position = 'none') ,
                   p.infofdis20cm + theme(legend.position = 'none') ,
                   p.infofdis1m + theme(legend.position = 'none') ,
                   p.infofdis5m + theme(legend.position = 'none') , 
                   labels = c('a)', 'b)', 'c)', 'd)'), nrow = 1)

plot_grid(p.row, legend, nrow = 2, rel_heights = c(4, 0.5))


