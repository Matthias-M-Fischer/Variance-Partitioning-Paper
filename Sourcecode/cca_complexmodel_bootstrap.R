# Load required libraries
require("vegan")
require("mvtnorm")

################################ H e l p e r   f u n c t i o n s  #################################

generate.population = function(nb.sites = 1e+6,
                               mu.species1, mu.species2, mu.species3, mu.species4, mu.species5,
                               Sigma.species1, Sigma.species2, Sigma.species3, Sigma.species4, Sigma.species5,
                               sigma.noise = 0.01) {

  # Generate population of nb.sites sites
  sites = matrix(nrow=nb.sites, ncol=9)
  colnames(sites) = c("predictor1", "predictor2", "predictor3", "predictor4", "abundance1", "abundance2", "abundance3", "abundance4", "abundance5")
  
  # Draw value of environmental predictor variables
  sites[, "predictor1"] = runif(nb.sites, 0, 1)
  sites[, "predictor2"] = runif(nb.sites, 0, 1)
  sites[, "predictor3"] = runif(nb.sites, 0, 1)
  sites[, "predictor4"] = runif(nb.sites, 0, 1)
  
  # Simulate relative abundances
  sites[, "abundance1"] = dmvnorm(x = sites[, 1:4], mean = mu.species1, sigma = Sigma.species1) + rnorm(nb.sites, 0, sigma.noise)
  sites[, "abundance2"] = dmvnorm(x = sites[, 1:4], mean = mu.species2, sigma = Sigma.species2) + rnorm(nb.sites, 0, sigma.noise)
  sites[, "abundance3"] = dmvnorm(x = sites[, 1:4], mean = mu.species3, sigma = Sigma.species3) + rnorm(nb.sites, 0, sigma.noise)
  sites[, "abundance4"] = dmvnorm(x = sites[, 1:4], mean = mu.species4, sigma = Sigma.species4) + rnorm(nb.sites, 0, sigma.noise)
  sites[, "abundance5"] = dmvnorm(x = sites[, 1:4], mean = mu.species5, sigma = Sigma.species5) + rnorm(nb.sites, 0, sigma.noise)
  
  # Abundances below zero will get replaced by zero
  sites[, "abundance1"][sites[, "abundance1"] < 0] = 0
  sites[, "abundance2"][sites[, "abundance2"] < 0] = 0
  sites[, "abundance3"][sites[, "abundance3"] < 0] = 0
  sites[, "abundance4"][sites[, "abundance4"] < 0] = 0
  sites[, "abundance5"][sites[, "abundance5"] < 0] = 0
  
  # We need to make sure we only return samples, in which all community table row sums are greater than zero
  # We therefore add a pseudocount of 1e-5 to all abundance values
  sites[, 5:9] = sites[, 5:9] + 1e-5
  
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
    cca.model.predictor1 = cca(table[, 5:9], table[, 1:2])
    cca.model.predictor2 = cca(table[, 5:9], table[, 3:4])
    cca.model.bothpredictors = cca(table[, 5:9] , table[, 1:2], table[, 3:4])
    
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
########################## S i m u l a t i o n   p a r a m e t e r s  #############################

# Number of the samples from the population of sites the analysis of a scenario is based on for assessing the variability of the estimates
nb.iterations = 100

for (add_noise in c(0,0.01,0.05,0.1)) {

  # SD of noise added to abundance data
  sigma.noise = add_noise
  
  for (scenario in (1:100)) {
    
    ########################## S i m u l a t i o n   s t a r t s   h e r e  ###########################
    # Draw species parameters
    mu.species1 = runif(4, 0, 1)
    mu.species2 = runif(4, 0, 1)
    mu.species3 = runif(4, 0, 1)
    mu.species4 = runif(4, 0, 1)
    mu.species5 = runif(4, 0, 1)
    Sigma.species1 = matrix(c(runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1)), nrow=4)
    Sigma.species2 = matrix(c(runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1)), nrow=4)
    Sigma.species3 = matrix(c(runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1)), nrow=4)
    Sigma.species4 = matrix(c(runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1)), nrow=4)
    Sigma.species5 = matrix(c(runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1), 0, 0, 0, 0, runif(1, 0, 1)), nrow=4)
    
    # Generate a population of sites
    population = generate.population(nb.sites=10^6,
                                     mu.species1, mu.species2, mu.species3, mu.species4, mu.species5,
                                     Sigma.species1, Sigma.species2, Sigma.species3, Sigma.species4, Sigma.species5,
                                     sigma.noise = sigma.noise)
    
    # Calculate "true" VP estimates, based on the complete population of all sites
    populationlevel.vp = calculate.r2(population)
    populationlevel.r2.predictor12 = populationlevel.vp[1]
    populationlevel.r2.predictor34 = populationlevel.vp[2]
    populationlevel.r2.shared = populationlevel.vp[3]
    populationlevel.r2.residual = populationlevel.vp[4]
    
    # We store the estimates based on the sampled subsets here
    r2s.predictor12 = c()
    r2s.predictor34 = c()
    r2s.shared = c()
    r2s.residual = c()
    
    for (nb.iteration in (1:nb.iterations)) {
      
      # Sample nb.sampled.sites from it. Only sample sites, where predictor2 exceed given minimum value
      sample = sample.subset(population, nb.sampled.sites = 1000, minimum.predictor2 = 0)
      
      # Calculate VP estimates for the sample
      sample.vp = calculate.r2(sample)
      
      # And store them as relative values for later analysis
      r2s.predictor12 = c(r2s.predictor12, sample.vp[1])
      r2s.predictor34 = c(r2s.predictor34, sample.vp[2])
      r2s.shared = c(r2s.shared, sample.vp[3])
      r2s.residual = c(r2s.residual, sample.vp[4])
    }
    
    # Take another sample from the population
    sample = sample.subset(population, nb.sampled.sites = 1000, minimum.predictor2 = 0)
    
    # And bootstrap it
    bootstrap.results = boot.cca(sample)
    r2.bs = bootstrap.results[1]
    sd.bs = bootstrap.results[2]
    ci.lower.bs = bootstrap.results[3]
    ci.upper.bs = bootstrap.results[4]
    
    ######################## R e s u l t   o f   c u r r e n t   s c e n a r i o  ###################
    
    cat(scenario, ",",
        sigma.noise, ",",
        populationlevel.r2.predictor34, ",",
        
        mean(r2s.predictor34), ",", 
        sd(r2s.predictor34), ",",
        
        r2.bs,  ",", 
        sd.bs,  ",", 
        ci.lower.bs,  ",", 
        ci.upper.bs, "\n")
  }
}
