require("vegan")
data("BCI")
data("BCI.env")

# Helper functions
extract.total.variation = function(cca.obj) {
  return (summary(cca.obj)$tot.chi)
}

extract.residual.variation = function(cca.obj) {
  return (summary(cca.obj)$unconst.chi)
}

extract.constrained.variation = function(cca.obj) {
  return (summary(cca.obj)$constr.chi)
}

calculate.r2 = function(table, pred1, pred2) {
  # Fit CCA models
  cca.model.predictor1 = cca(table, pred1)
  cca.model.predictor2 = cca(table, pred2)
  cca.model.bothpredictors = cca(table, pred1, pred2)
  
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

# Import data
# We need to split the BCI.env data into predictor 1 (environmental data without coordinates) and predictor 2 (only the coordinates)
table = BCI
pred1 = BCI.env[, 3:ncol(BCI.env)]
pred2 = BCI.env[, 1:2]

# The first four of the remaining columns of the environmental predictors do not show any variation, so we remove them
pred1[, 1:4] = NULL

# Print column names
cat ("#nb.repetition, mean.R2.env, sd.R2.env, ci.lower.env, ci.upper.env, mean.R2.spat, sd.R2.spat, ci.lower.spat, ci.upper.spat, mean.R2.resid, sd.R2.resid, ci.lower.resid, ci.upper.resid \n")

# Repeat the bootstrap analysis 10 times
for (i in (1:10)) {
  
  # Store the M estimates of each analysis here
  est_env = c()
  est_spat = c()
  est_shared = c()
  est_resid = c()
  
  #M = 1000 draws per repetition
  for (M in (1:1000)) {
    # Which lines to sample to obtain the subsample
    s = sample(c(1:nrow(table)), replace=T)
    
    # Calculate the adjusted R2 values for the subsample
    estimates = calculate.r2(table[s, ], pred1[s, ], pred2[s, ])
    
    # Extract the adj. R2 values
    pred1_r2 = estimates[1]
    pred2_r2 = estimates[2]
    shared_r2 = estimates[3]
    resid_r2 = estimates[4]
    
    # And store them
    est_env = c(est_env, pred1_r2)
    est_spat = c(est_spat, pred2_r2)
    est_shared = c(est_shared, shared_r2)
    est_resid = c(est_resid, resid_r2)
  }
  
  # Print result for the i-th repetition
  cat(i, ",", mean(est_env, na.rm=T), ",", sd(est_env, na.rm=T) , ",", quantile(est_env, c(0.025), na.rm=T), ",", quantile(est_env, c(0.975), na.rm=T), ",",
      mean(est_spat, na.rm=T), ",", sd(est_spat, na.rm=T), ",", quantile(est_spat, c(0.025), na.rm=T), ",", quantile(est_spat, c(0.975), na.rm=T), ",",
      mean(est_shared, na.rm=T), ",", sd(est_shared, na.rm=T) , ",", quantile(est_shared, c(0.025), na.rm=T), ",", quantile(est_shared, c(0.975), na.rm=T), ",",
      mean(est_resid,na.rm=T), ",", sd(est_resid, na.rm=T) , ",", quantile(est_resid, c(0.025), na.rm=T), ",", quantile(est_resid, c(0.975), na.rm=T), "\n")
}