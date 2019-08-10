# Load required libraries
require("vegan")
require("mvtnorm")

################################ H e l p e r   f u n c t i o n s  #################################

generate.population = function(nb.sites = 1e+6,
                                mu.species1 = c(0,0),
                                mu.species2 = c(1,1),
                                Sigma.species1 = matrix(c(0.5, 0, 0, 0.5), nrow=2),
                                Sigma.species2 = matrix(c(0.5, 0, 0, 0.5), nrow=2),
                                sigma.noise = 0.01) {
  # Generate population of nb.sites sites
  # mu.species1 and mu.species2 denote the maximum values of the bivariate response curves of the two species
  # Sigma.species1 and Sigma.species2 denote the covariance matrices of the bivariate response curves
  # sigma.noise is the SD of the white noise added to the relative abundance values
  
  sites = matrix(nrow=nb.sites, ncol=4)
  colnames(sites) = c("predictor1", "predictor2", "abundance1", "abundance2")
  
  # Draw value of environmental predictor variables
  sites[, "predictor1"] = runif(nb.sites, 0, 1)
  sites[, "predictor2"] = runif(nb.sites, 0, 1)
  
  # Simulate relative abundances
  sites[, "abundance1"] = dmvnorm(x = sites[, 1:2], mean = mu.species1, sigma = Sigma.species1) + rnorm(nb.sites, 0, sigma.noise)
  sites[, "abundance2"] = dmvnorm(x = sites[, 1:2], mean = mu.species2, sigma = Sigma.species2) + rnorm(nb.sites, 0, sigma.noise)
  
  # Abundances below zero will get replaced by zero
  sites[, "abundance1"][sites[, "abundance1"] < 0] = 0
  sites[, "abundance2"][sites[, "abundance2"] < 0] = 0
  
  # We need to make sure we only return samples, in which all community table row sums are greater than zero
  # We therefore add a pseudocount of 1e-5 to all abundance values
  sites[, 3:4] = sites[, 3:4] + 1e-5
  
  return (as.data.frame(sites))
}

sample.subset = function(population, nb.sampled.sites = 100, minimum.predictor2 = 0) {

  # Sample a subset of sites. Suitable for sampling are only those sites, where predictor #2
  # extends a certain minimum value, in oder to be able to analyse, how sampling range
  # may influence VP uncertainty.
  suitable.sites = population[population$predictor2 > minimum.predictor2, ]
  sampled.row.indices = sample((1:nrow(suitable.sites)), nb.sampled.sites, replace=FALSE)
  sampled.sites = suitable.sites[sampled.row.indices, ]
}

extract.total.variation = function(cca.obj) {
  return (summary(cca.obj)$tot.chi)
}

extract.residual.variation = function(cca.obj) {
  return (summary(cca.obj)$unconst.chi)
}

extract.constrained.variation = function(cca.obj) {
  return (summary(cca.obj)$constr.chi)
}

calculate.r2 = function(table) {
  # Fit CCA models
  cca.model.predictor1 = cca(table[, 3:4], table$predictor1)
  cca.model.predictor2 = cca(table[, 3:4], table$predictor2)
  cca.model.bothpredictors = cca(table[, 3:4], table$predictor1, table$predictor2)
  
  # Extract data from CCA models
  # http://pbil.univ-lyon1.fr/members/dray/files/articles/peres-neto2006a.pdf
  variation.total = extract.total.variation(cca.model.bothpredictors)
  variation.residual = extract.residual.variation(cca.model.bothpredictors)
  variation.predictor1 = extract.constrained.variation(cca.model.predictor1)
  variation.predictor2 = extract.constrained.variation(cca.model.predictor2)
  
  variation.union = variation.total - variation.residual
  variation.only.predictor1 = variation.union - variation.predictor2
  variation.only.predictor2 = variation.union - variation.predictor1
  variation.shared = variation.union - variation.only.predictor1 - variation.only.predictor2
  
  # Compute levels of relative explained variation
  r2.only.predictor1 = variation.only.predictor1 / variation.total
  r2.only.predictor2 = variation.only.predictor2 / variation.total
  r2.shared = variation.shared/variation.total
  r2.residual = variation.residual/variation.total
  
  return(c(r2.only.predictor1,
           r2.only.predictor2,
           r2.shared,
           r2.residual))
}

calculate.r2adj = function(table) {
  
  # Fit CCA models
  cca.model.predictor1 = cca(table[, 3:4], table$predictor1)
  cca.model.predictor2 = cca(table[, 3:4], table$predictor2)
  cca.model.bothpredictors = cca(table[, 3:4], table$predictor1, table$predictor2)
  
  # Extract data from CCA models
  # http://pbil.univ-lyon1.fr/members/dray/files/articles/peres-neto2006a.pdf
  variation.total = extract.total.variation(cca.model.bothpredictors)
  variation.residual = extract.residual.variation(cca.model.bothpredictors)
  variation.predictor1 = extract.constrained.variation(cca.model.predictor1)
  variation.predictor2 = extract.constrained.variation(cca.model.predictor2)
  variation.union = extract.constrained.variation(cca.model.bothpredictors)
  variation.shared = variation.union - variation.predictor1 - variation.predictor2
  
  # Compute levels of relative explained variation at population-level and store for later comparision
  r2.predictor1 = variation.predictor1 / variation.total
  r2.predictor2 = variation.predictor2 / variation.total
  r2.shared = variation.shared/variation.total
  r2.residual = variation.residual/variation.total
  
  # Calculate correction term for adjustment, following the paper cited above
  n = nrow(table)
  p = ncol(table) - 2
  correction = 1 - (n-1)/(n-p-1)
  
  return(c(r2.predictor1 + correction,
           r2.predictor2 + correction,
           r2.shared + correction,
           r2.residual + correction))
}

boot.cca = function(table, nb_bootstraps = 1000) {
  
  estimates = c()
  for (i in (1:nb_bootstraps)) {
    s = sample(c(1:nrow(table)), replace=T)
    r2 = calculate.r2(table[s, ])[2]
    estimates = c(estimates, r2)
  }
  
  ret = c(mean(estimates), sd(estimates), quantile(estimates, 0.025, na.rm=T), quantile(estimates, 0.975, na.rm=T))
  names(ret) = c("mu R^2", "s R^2", "CI.lower", "CI.upper")
  return(ret)
}

############################ S i m u l a t i n g   s c e n a r i o s  #############################

cat ("param, sd.noise, populationlevel.r2, mean.r2.resampling, sd.r2.resampling, bs.r2, bs.sd, bs.lower, bs.upper \n")

# 'param' will be changed throughout the simulated scenarios
#for (param in c(rep(0, 10), rep(0.1, 10), rep(0.2, 10), rep(0.3, 10), rep(0.4, 10), rep(0.5, 10), rep(0.6, 10), rep(0.7, 10), rep(0.8, 10), rep(0.9, 10), rep(1.0, 10))) {
for (param in c(rep(5,10), rep(10,10), rep(25,10), rep(50,10), rep(75,0), rep(100,10), rep(250, 10), rep(500, 10), rep(750, 10), rep(1000, 10))) {
  
  ########################## S i m u l a t i o n   p a r a m e t e r s  #############################
  
  # Number of the samples from the population of sites the analysis of a scenario is based on for assessing the variability of the estimates
  nb.iterations = 100
  
  # Size of every subsample drawn from the population of all sites
  nb.sampled.sites = param
    
  # Sample only sites where the second predictor exceeds this minimum
  minimum.predictor2 = 0
    
  # Difference in response curve optima of the two species for the secound environmental predictor
  responsecurve.difference = 1.0
      
  # SD of noise added to abundance data
  sigma.noise = 0.05
  

  ########################## S i m u l a t i o n   s t a r t s   h e r e  ###########################
  
  # Generate a population of sites
  population = generate.population(mu.species2 = c(1,0 + responsecurve.difference), sigma.noise = sigma.noise)
  
  # Calculate "true" VP estimates, based on the complete population of all sites
  populationlevel.vp = calculate.r2(population)
  populationlevel.r2.predictor1 = populationlevel.vp[1]
  populationlevel.r2.predictor2 = populationlevel.vp[2]
  populationlevel.r2.shared = populationlevel.vp[3]
  populationlevel.r2.residual = populationlevel.vp[4]
  
  # We store the estimates based on the sampled subsets here
  r2s.predictor1 = c()
  r2s.predictor2 = c()
  r2s.shared = c()
  r2s.residual = c()
  
  for (nb.iteration in (1:nb.iterations)) {
    
    # Sample nb.sampled.sites from it. Only sample sites, where predictor2 exceed given minimum value
    sample = sample.subset(population, nb.sampled.sites = nb.sampled.sites, minimum.predictor2 = minimum.predictor2)
    
    # Calculate VP estimates for the sample
    sample.vp = calculate.r2(sample)
    
    # And store them for later analysis
    r2s.predictor1 = c(r2s.predictor1, sample.vp[1])
    r2s.predictor2 = c(r2s.predictor2, sample.vp[2])
    r2s.shared = c(r2s.shared, sample.vp[3])
    r2s.residual = c(r2s.residual, sample.vp[4])
  }
  
    # Take another sample from the population
    sample = sample.subset(population, nb.sampled.sites = nb.sampled.sites, minimum.predictor2 = minimum.predictor2)
    
    # And bootstrap it
    bootstrap.results = boot.cca(sample)
    r2.bs = bootstrap.results[1]
    sd.bs = bootstrap.results[2]
    ci.lower.bs = bootstrap.results[3]
    ci.upper.bs = bootstrap.results[4]
  
  ######################## R e s u l t   o f   c u r r e n t   s c e n a r i o  ###################
  
  cat(param, ",",
      sigma.noise, ",",
      populationlevel.r2.predictor2, ",",
      
      mean(r2s.predictor2), ",", 
      sd(r2s.predictor2), ",",

      r2.bs,  ",", 
      sd.bs,  ",", 
      ci.lower.bs,  ",", 
      ci.upper.bs, "\n")
}
