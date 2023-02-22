library('BAMMtools')
library(coda)

getwd()
setwd('~/OneDrive - Imperial College London/TreeBuilding/230203filter/bamm/longrun/')

# Create bammdata object
tree <- read.tree("ultrametric_codon.tre")

# Assess MCMC convergence
mcmcout <- read.csv("codon_genus_mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

# Discard % burnin
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
plot(postburn$logLik ~ postburn$generation)


# Check the effective sample sizes of the log-likelihood and the number of shift 
# events present in each sample. Want both at least 200, and more for smaller datasets
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

# compute the posterior probabilities of models sampled using BAMM (higher better)
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)
shift_probs <- summary(edata)

# Assess prior effect
bfmat <- computeBayesFactors(mcmcout, expectedNumberOfShifts=5, burnin=0.25)
bfmat

# Get data, choosing burnin
edata <- getEventData(tree, eventdata = "codon_overall_event_data.txt", burnin=0.1)

# Phylorate Plot
summary(edata)
plot.bammdata(edata, lwd=2)
# lwd = line weighting
# cex = tip label size
# colour breaks: default is linear. Can also choose breaksmethod = "quantile" or "jenks"
plot.bammdata(edata, lwd=2, legend=T, labels=TRUE, breaksmethod = 'quantile')
plot.bammdata(edata, lwd=2, legend=T, labels=TRUE, breaksmethod = 'linear')
plot.bammdata(edata, lwd=2, legend=T, labels=TRUE, breaksmethod = 'jenks')


# Credible sets of shift configurations
css <- credibleShiftSet(edata, expectedNumberOfShifts=10, threshold=5, set.limit = 0.95)
css$number.distinct
summary(css)
# Plots most credible shift configurations. f is posterior probability (higher better)
plot.credibleshiftset(css)

# Finding the single best shift configuration
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=10)
plot.bammdata(best, lwd=2, legend=T, labels=TRUE, breaksmethod = 'quantile')
plot.bammdata(best, lwd=2, legend=T, labels=TRUE, breaksmethod = 'linear')
plot.bammdata(best, lwd=2, legend=T, labels=TRUE, breaksmethod = 'jenks')addBAMMshifts(best, cex=2.5)
plot.bammdata(best, cex = 1, breaksmethod = 'jenks')
plot.bammshifts(best)

