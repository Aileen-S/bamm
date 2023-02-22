#!/usr/bin/env Rscript

# Load packages
library('BAMMtools')
library(coda)
library(getopt)

# Get command line arguments
# column 1 = long flag
# column 2 = short flag
# column 3 = argument mask (0=no argument, 1=required argument, 2=optional argument)
# column 4 = argument type (logical, integer, double, complex, character)
# column 5 = optional, description/help message
spec <- matrix(c(
  'tree', 't', 1, 'character', 'Path to tree file',
  'event', 'e', 1, 'character', 'Path to event_data csv file',
  'burnin', 'b', 1, 'integer', 'Burnin as decimal fraction, after assessing mcmc plot',
  'output', 'o', 1, 'character', 'Prefix for output png file'
  ), byrow = T, ncol = 5)

opt <- getopt(spec)

# Create bammdata object
tree <- read.tree(opt$tree)

# Get data, choosing burnin
edata <- getEventData(tree, eventdata = opt$event, burnin = opt$burnin)

# Phylorate Plot
# lwd = line weighting
# cex = tip label size

# Export quantile plot
png(filename = paste(opt$output, '_quantile.png', sep = ''), width = 3000, height = 15000)
plot.bammdata(edata, lwd=2, legend=T, labels=TRUE, breaksmethod = 'quantile')
dev.off()

# Export linear plot
png(filename = paste(opt$output, '_linear.png', sep = ''), width = 3000, height = 15000)
plot.bammdata(edata, lwd=2, legend=T, labels=TRUE, breaksmethod = 'linear')
dev.off()

