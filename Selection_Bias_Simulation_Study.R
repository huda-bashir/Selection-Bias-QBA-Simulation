###################################
## SELECTION BIAS SIMULAITION STUDY
## AUTHOR: ANONYMIZED FOR REVIEW
###################################


#### LOADING PACKAGES ######
pacman::p_load(dplyr, tidyverse, future,future.apply, triangulr, trapezoid, beta, ggpubr)

##### FIRST SIMULATION, UNIFORM DISTRIBUTION ############
## Writing a function (SIM_FUN_UNIF) for the simulation
SIM_FUN_UNIF <- function(niter, ssize, RR, prevE, baserisk, prD1E1, prD1E0, prD0E1, prD0E0, nsubiter) {
  
  #### COHORT SIMULATION DATA GENERATION
  E1 <- rbinom(niter, ssize, prevE) # total exposed
  E0 <- ssize - E1 # total unexposed
  D1E1 <- rbinom(niter, E1, RR*baserisk) # exposed cases
  D1E0 <- rbinom(niter, E0, baserisk) # unexposed cases
  
  #### CALCULATING TRUE RR ESTIMATE FOR COHORT
  RR1 <- (D1E1/E1)/(D1E0/E0)
  
  #### DATA GENERATION FOR CASE-COHORT ESTIMATE
  D0 <- (D1E1 + D1E0)*2 # 2 controls per case
  D0E1 <- rbinom(niter, D0, prevE) # exposed controls
  D0E0 <- D0 - D0E1 # unexposed controls
  
  #### CALCULATING OR ESTIMATE FOR CASE-COHORT
  OR <- (D1E1/D0E1)/(D1E0/D0E0)
  L_SERR <- sqrt(1/D1E1+1/D0E1+1/D1E0+1/D0E0)
  ME_OR <- qnorm(0.975)*L_SERR
  LCL95 <- OR*exp(-ME_OR)
  UCL95 <- OR*exp(ME_OR)
  MSE <- mean((log(OR) - log(RR))^2)
  coverage <- ifelse(LCL95 <= RR & UCL95 >= RR, 1, 0)
  
  #### EXPOSURE PREVALENCE CASES AND CONTROLS
  prevE_D1 <- D1E1/(D1E1+D1E0)
  prevE_D0 <- D0E1/(D0E1+D0E0)
  
  #### DATA GENERATION FOR BIASED ESTIMATE
  pD1E1 <- rbinom(niter, D1E1, prD1E1)
  pD1E0 <- rbinom(niter, D1E0, prD1E0)
  pD0E1 <- rbinom(niter, D0E1, prD0E1)
  pD0E0 <- rbinom(niter, D0E0, prD0E0)
  
  #### CALCULATING BIASED OR ESTIMATE
  OR_b <- (pD1E1/pD0E1)/(pD1E0/pD0E0) 
  median(OR_b)
  L_SERR_b <- sqrt(1/pD1E1+1/pD0E1+1/pD1E0+1/pD0E0)
  ME_OR_b <- qnorm(0.975)*L_SERR_b
  LCL95_b <- OR_b*exp(-ME_OR_b)
  UCL95_b <- OR_b*exp(ME_OR_b)
  MSE_b <- mean((log(OR_b) - log(RR))^2)
  coverage_b <- ifelse(LCL95_b <= RR & UCL95_b >= RR, 1, 0)
  
  #### DATA GENERATION FOR RANDOM UNIFORM ADJUSTMENT
  ## TOTAL NON-PARTICIPATING CASES & CONTROLS BY EXPOSURE STATUS
  npD1E1 <- D1E1-pD1E1
  npD1E0 <- D1E0-pD1E0
  npD0E1 <- D0E1-pD0E1
  npD0E0 <- D0E0-pD0E0
  
  ## TOTAL NON-PARTICIPATING CASES & CONTROLS
  npD1 <- npD1E1 + npD1E0
  npD0 <- npD0E1 + npD0E0
  
  ## EXPOSURE PREVALENCE NON-PARTICIPATING CASES & CONTROLS
  prevE_np_D1 <- npD1E1/npD1
  prevE_np_D0 <- npD0E1/npD0
  
  ## EXPOSURE PREVALENCE PARTICIPATING CASES & CONTROLS
  prevE_p_D1 <- pD1E1/(pD1E1 + pD1E0)
  prevE_p_D0 <- pD0E1/(pD0E1 + pD0E0)
  
  ## SET UP PARALLEL PROCESSING
  plan(multisession, workers = parallel::detectCores() - 1)
  
  #### BIAS-ADJUSTED ESTIMATE
  CL95_list_a <- future_lapply(1:niter, function(i) {
    
    #### EXPOSURE PREVALENCE - 0-1 RANDOM UNIFORM
    prevE_np_D1_r <- runif(nsubiter, min = 0, max = 1)
    prevE_np_D0_r <- runif(nsubiter, min = 0, max = 1)
    
    npD1E1_a <- rbinom(nsubiter, npD1[i], prevE_np_D1_r) # NONPARTICIPATING EXPOSED CASES
    npD1E0_a <- npD1[i] - npD1E1_a # NONPARTICIPATING UNEXPOSED CASES
    npD0E1_a <- rbinom(nsubiter, npD0[i], prevE_np_D0_r) # NONPARTICIPATING EXPOSED CONTROLS
    npD0E0_a <- npD0[i] - npD0E1_a # NONPARTICIPATING UNEXPOSED CONTROLS
    
    #### CALCULATING BIAS-ADJUSTED OR FOR RANDOM UNIFORM
    OR_a_vec <- ((pD1E1[i] + npD1E1_a) / (pD1E0[i] + npD1E0_a)) / 
      ((pD0E1[i] + npD0E1_a) / (pD0E0[i] + npD0E0_a))
    L_SERR_a <- sqrt((1/(pD1E1[i]))+(1/(pD1E0[i])) + 
                        (1/(pD0E1[i])) + (1/(pD0E0[i])))
    OR_a_vec <- exp(log(OR_a_vec) + rnorm(nsubiter) * L_SERR_a)
    
    #### CALCULATE 95% SI
    quantiles1 <- quantile(OR_a_vec, probs=c(0.025, 0.5, 0.975))
    c(quantiles1, ifelse(quantiles1[1] <= RR & quantiles1[3] >= RR, 1, 0))
    
  }, future.seed = TRUE)
  
  CL95_a <- do.call(rbind, CL95_list_a)
  
  ## SAVING RESULTS
  res <- data.frame(RR1, OR, LCL95, UCL95, MSE, coverage, OR_b, LCL95_b, UCL95_b, MSE_b, coverage_b, prevE_D1, prevE_D0, prevE_p_D1, prevE_p_D0, prevE_np_D1, prevE_np_D0)
  res_a <- as.data.frame(CL95_a)
  colnames(res_a) <- c("LCL95_a", "median_OR_a", "UCL95_a","coverage_a")
  
  return(list(res = res, res_a = res_a))
}

#### SIMULATION SCENARIO ONE:
# Cases: Marginal = 85%
# D1E1 = 75%
# D1E0 = 95%
# 
# Controls: Marginal = 85%
# D0E1 = 82%
# D0E0 = 88%

true_RR = 2
set.seed(seed = 23) 

sim1.1 <- SIM_FUN_UNIF(niter = 100000, ssize = 40000, RR = 2, prevE = 0.2, baserisk = 0.01, prD1E1 = 0.75, prD1E0 = 0.95, prD0E1 = 0.82, prD0E0 = 0.88, nsubiter = 10000)

## Extracting results
sim1.1_res <- sim1.1$res
sim1.1_res_a <- sim1.1$res_a

## true OR
quantile(sim1.1_res$OR, probs=c(0.025, 0.5, 0.975))
median(sim1.1_res$MSE)
sum(sim1.1_res$coverage)/length(sim1.1_res$coverage)

## biased OR
quantile(sim1.1_res$OR_b, probs=c(0.025, 0.5, 0.975))
median(sim1.1_res$MSE_b)
sum(sim1.1_res$coverage_b)/length(sim1.1_res$coverage_b)

## bias-adjusted OR
quantile(sim1.1_res_a$median_OR_a, probs=c(0.025, 0.5, 0.975))
sim1.1_MSE_a <- (mean((log(sim1.1_res_a$median_OR_a) - log(true_RR))^2))
sim1.1_MSE_a
sum(sim1.1_res_a$coverage_a)/length(sim1.1_res_a$coverage_a)

#### SIMULATION SCENARIO TWO:
# Cases: Marginal = 85%
# D1E1 = 82%
# D1E0 = 88%
# Controls: Marginal = 85%
# D0E1 = 75%
# D0E0 = 95%

set.seed(seed = 23) 
sim2.1 <- SIM_FUN_UNIF(niter = 100000, ssize = 40000, RR = 2, prevE = 0.2, baserisk = 0.01 , prD1E1 = 0.82, prD1E0 = 0.88, prD0E1 = 0.75, prD0E0 = 0.95, nsubiter = 10000)

sim2.1_res <- sim2.1$res
sim2.1_res_a <- sim2.1$res_a

## true OR
quantile(sim2.1_res$OR, probs=c(0.025, 0.5, 0.975))
median(sim2.1_res$MSE)
sum(sim2.1_res$coverage)/length(sim2.1_res$coverage)

## biased OR
quantile(sim2.1_res$OR_b, probs=c(0.025, 0.5, 0.975))
median(sim2.1_res$MSE_b)
sum(sim2.1_res$coverage_b)/length(sim2.1_res$coverage_b)

## bias-adjusted OR
quantile(sim2.1_res_a$median_OR_a, probs=c(0.025, 0.5, 0.975))
sim2.1_MSE_a <- mean((log(sim2.1_res_a$median_OR_a) - log(true_RR))^2)
sim2.1_MSE_a
sum(sim2.1_res_a$coverage_a)/length(sim2.1_res_a$coverage_a)

##### SECOND SIMULATION, TRIANGULAR DISTRIBUTION ############
## Writing a function (SIM_FUN_TRI) for the simulation
SIM_FUN_TRI <- function(niter, ssize, RR, prevE, baserisk, prD1E1, prD1E0, prD0E1, prD0E0, nsubiter) {
  
  #### COHORT SIMULATION DATA GENERATION
  E1 <- rbinom(niter, ssize, prevE) # total exposed
  E0 <- ssize - E1 # total unexposed
  D1E1 <- rbinom(niter, E1, RR*baserisk) # exposed cases
  D1E0 <- rbinom(niter, E0, baserisk) # unexposed cases
  
  #### CALCULATING TRUE RR ESTIMATE FOR COHORT
  RR1 <- (D1E1/E1)/(D1E0/E0)
  
  #### DATA GENERATION FOR CASE-COHORT ESTIMATE
  D0 <- (D1E1 + D1E0)*2 # 2 controls per case
  D0E1 <- rbinom(niter, D0, prevE) # exposed controls
  D0E0 <- D0 - D0E1 # unexposed controls
  
  #### CALCULATING OR ESTIMATE FOR CASE-COHORT
  OR <- (D1E1/D0E1)/(D1E0/D0E0)
  L_SERR <- sqrt(1/D1E1+1/D0E1+1/D1E0+1/D0E0)
  ME_OR <- qnorm(0.975)*L_SERR
  LCL95 <- OR*exp(-ME_OR)
  UCL95 <- OR*exp(ME_OR)
  MSE <- mean((log(OR) - log(RR))^2)
  coverage <- ifelse(LCL95 <= RR & UCL95 >= RR, 1, 0)
  
  #### EXPOSURE PREVALENCE CASES AND CONTROLS
  prevE_D1 <- D1E1/(D1E1+D1E0)
  prevE_D0 <- D0E1/(D0E1+D0E0)
  
  #### DATA GENERATION FOR BIASED ESTIMATE
  pD1E1 <- rbinom(niter, D1E1, prD1E1)
  pD1E0 <- rbinom(niter, D1E0, prD1E0)
  pD0E1 <- rbinom(niter, D0E1, prD0E1)
  pD0E0 <- rbinom(niter, D0E0, prD0E0)
  
  #### CALCULATING BIASED OR ESTIMATE
  OR_b <- (pD1E1/pD0E1)/(pD1E0/pD0E0) 
  median(OR_b)
  L_SERR_b <- sqrt(1/pD1E1+1/pD0E1+1/pD1E0+1/pD0E0)
  ME_OR_b <- qnorm(0.975)*L_SERR_b
  LCL95_b <- OR_b*exp(-ME_OR_b)
  UCL95_b <- OR_b*exp(ME_OR_b)
  MSE_b <- mean((log(OR_b) - log(RR))^2)
  coverage_b <- ifelse(LCL95_b <= RR & UCL95_b >= RR, 1, 0)
  
  #### DATA GENERATION FOR RANDOM TRIANGULAR ADJUSTMENT
  ## TOTAL NON-PARTICIPATING CASES & CONTROLS BY EXPOSURE STATUS
  npD1E1 <- D1E1-pD1E1
  npD1E0 <- D1E0-pD1E0
  npD0E1 <- D0E1-pD0E1
  npD0E0 <- D0E0-pD0E0
  
  ## TOTAL NON-PARTICIPATING CASES & CONTROLS
  npD1 <- npD1E1 + npD1E0
  npD0 <- npD0E1 + npD0E0
  
  ## EXPOSURE PREVALENCE NON-PARTICIPATING CASES & CONTROLS
  prevE_np_D1 <- npD1E1/npD1
  prevE_np_D0 <- npD0E1/npD0
  
  ## EXPOSURE PREVALENCE PARTICIPATING CASES & CONTROLS
  prevE_p_D1 <- pD1E1/(pD1E1 + pD1E0) # mode for rtri - cases
  prevE_p_D0 <- pD0E1/(pD0E1 + pD0E0)
  
  ## SET UP PARALLEL PROCESSING
  plan(multisession, workers = parallel::detectCores() - 1)
  
  #### BIAS-ADJUSTED ESTIMATE
  CL95_list_a <- future_lapply(1:niter, function(i) {
    
    #### EXPOSURE PREVALENCE - PARTICIPATING CASES - CENTERED ON PREVE_P_D1 & PREVE_P_D0
    #### EXPOSURE PREVALENCE - PARTICIPATING CASES - MIN = 0, MAX = 1
    prevE_np_D1_t <- rtri(nsubiter, min = 0, max = 1, mode = prevE_p_D1[i])
    prevE_np_D0_t <- rtri(nsubiter, min = 0, max = 1, mode = prevE_p_D0[i])
    
    npD1E1_a <- rbinom(nsubiter, npD1[i], prevE_np_D1_t) # NONPARTICIPATING EXPOSED CASES
    npD1E0_a <- npD1[i] - npD1E1_a # NONPARTICIPATING UNEXPOSED CASES
    npD0E1_a <- rbinom(nsubiter, npD0[i], prevE_np_D0_t) # NONPARTICIPATING EXPOSED CONTROLS
    npD0E0_a <- npD0[i] - npD0E1_a # NONPARTICIPATING UNEXPOSED CONTROLS
    
    #### CALCULATING BIAS-ADJUSTED OR FOR RANDOM TRAPEZOIDAL
    OR_a_vec <- ((pD1E1[i] + npD1E1_a) / (pD1E0[i] + npD1E0_a)) / 
      ((pD0E1[i] + npD0E1_a) / (pD0E0[i] + npD0E0_a))
    L_SERR_a <- sqrt((1/(pD1E1[i]))+(1/(pD1E0[i])) + 
                        (1/(pD0E1[i])) + (1/(pD0E0[i])))
    OR_a_vec <- exp(log(OR_a_vec) + rnorm(nsubiter) * L_SERR_a)
    
    ### CALCULATE 95% SI
    quantiles1 <- quantile(OR_a_vec, probs=c(0.025, 0.5, 0.975))
    c(quantiles1, ifelse(quantiles1[1] <= RR & quantiles1[3] >= RR, 1, 0))
    
  }, future.seed = TRUE)
  
  CL95_a <- do.call(rbind, CL95_list_a)
  
  ## SAVING RESULTS
  res <- data.frame(RR1, OR, LCL95, UCL95, MSE, coverage, OR_b, LCL95_b, UCL95_b, MSE_b, coverage_b, prevE_D1, prevE_D0, prevE_p_D1, prevE_p_D0, prevE_np_D1, prevE_np_D0)
  res_a <- as.data.frame(CL95_a)
  colnames(res_a) <- c("LCL95_a", "median_OR_a", "UCL95_a","coverage_a")

  return(list(res = res, res_a = res_a))
}

#### SIMULATION SCENARIO ONE:
set.seed(seed = 23) 
sim1.2 <- SIM_FUN_TRI(niter = 100000, ssize = 40000, RR = 2, prevE = 0.2, baserisk = 0.01, prD1E1 = 0.75, prD1E0 = 0.95, prD0E1 = 0.82, prD0E0 = 0.88, nsubiter = 10000)

sim1.2_res <- sim1.2$res
sim1.2_res_a <- sim1.2$res_a

## true OR
quantile(sim1.2_res$OR, probs=c(0.025, 0.5, 0.975))
median(sim1.2_res$MSE)
sum(sim1.2_res$coverage)/length(sim1.2_res$coverage)

## biased OR
quantile(sim1.2_res$OR_b, probs=c(0.025, 0.5, 0.975))
median(sim1.2_res$MSE_b)
sum(sim1.2_res$coverage_b)/length(sim1.2_res$coverage_b)

## bias-adjusted OR
quantile(sim1.2_res_a$median_OR_a, probs=c(0.025, 0.5, 0.975))
sim1.2_MSE_a <- (mean((log(sim1.2_res_a$median_OR_a) - log(true_RR))^2))
sim1.2_MSE_a
sum(sim1.2_res_a$coverage_a)/length(sim1.2_res_a$coverage_a)

#### SIMULATION SCENARIO TWO:
set.seed(seed = 23)
sim2.2 <- SIM_FUN_TRI(niter = 100000, ssize = 40000, RR = 2, prevE = 0.2, baserisk = 0.01 , prD1E1 = 0.82, prD1E0 = 0.88, prD0E1 = 0.75, prD0E0 = 0.95, nsubiter = 10000)

sim2.2_res <- sim2.2$res
sim2.2_res_a <- sim2.2$res_a

## true OR
quantile(sim2.2_res$OR, probs=c(0.025, 0.5, 0.975))
median(sim2.2_res$MSE)
sum(sim2.2_res$coverage)/length(sim2.2_res$coverage)

## biased OR
quantile(sim2.2_res$OR_b, probs=c(0.025, 0.5, 0.975))
median(sim2.2_res$MSE_b)
sum(sim2.2_res$coverage_b)/length(sim2.2_res$coverage_b)

## bias-adjusted OR
quantile(sim2.2_res_a$median_OR_a, probs=c(0.025, 0.5, 0.975))
sim2.2_MSE_a <- mean((log(sim2.2_res_a$median_OR_a) - log(true_RR))^2)
sim2.2_MSE_a
sum(sim2.2_res_a$coverage_a)/length(sim2.2_res_a$coverage_a)

##### THIRD SIMULATION, TRAPEZOIDAL DISTRIBUTION ############
## Writing a function (SIM_FUN_TRAP) for the simulation
SIM_FUN_TRAP <- function(niter, ssize, RR, prevE, baserisk, prD1E1, prD1E0, prD0E1, prD0E0, nsubiter) {
  
  #### COHORT SIMULATION DATA GENERATION
  E1 <- rbinom(niter, ssize, prevE) # total exposed
  E0 <- ssize - E1 # total unexposed
  D1E1 <- rbinom(niter, E1, RR*baserisk) # exposed cases
  D1E0 <- rbinom(niter, E0, baserisk) # unexposed cases
  
  #### CALCULATING TRUE RR ESTIMATE FOR COHORT
  RR1 <- (D1E1/E1)/(D1E0/E0)
  
  #### DATA GENERATION FOR CASE-COHORT ESTIMATE
  D0 <- (D1E1 + D1E0)*2 # 2 controls per case
  D0E1 <- rbinom(niter, D0, prevE) # exposed controls
  D0E0 <- D0 - D0E1 # unexposed controls
  
  #### CALCULATING OR ESTIMATE FOR CASE-COHORT
  OR <- (D1E1/D0E1)/(D1E0/D0E0)
  L_SERR <- sqrt(1/D1E1+1/D0E1+1/D1E0+1/D0E0)
  ME_OR <- qnorm(0.975)*L_SERR
  LCL95 <- OR*exp(-ME_OR)
  UCL95 <- OR*exp(ME_OR)
  MSE <- mean((log(OR) - log(RR))^2)
  coverage <- ifelse(LCL95 <= RR & UCL95 >= RR, 1, 0)
  
  #### EXPOSURE PREVALENCE CASES AND CONTROLS
  prevE_D1 <- D1E1/(D1E1+D1E0)
  prevE_D0 <- D0E1/(D0E1+D0E0)
  
  #### DATA GENERATION FOR BIASED ESTIMATE
  pD1E1 <- rbinom(niter, D1E1, prD1E1)
  pD1E0 <- rbinom(niter, D1E0, prD1E0)
  pD0E1 <- rbinom(niter, D0E1, prD0E1)
  pD0E0 <- rbinom(niter, D0E0, prD0E0)
  
  #### CALCULATING BIASED OR ESTIMATE
  OR_b <- (pD1E1/pD0E1)/(pD1E0/pD0E0) 
  median(OR_b)
  L_SERR_b <- sqrt(1/pD1E1+1/pD0E1+1/pD1E0+1/pD0E0)
  ME_OR_b <- qnorm(0.975)*L_SERR_b
  LCL95_b <- OR_b*exp(-ME_OR_b)
  UCL95_b <- OR_b*exp(ME_OR_b)
  MSE_b <- mean((log(OR_b) - log(RR))^2)
  coverage_b <- ifelse(LCL95_b <= RR & UCL95_b >= RR, 1, 0)
  
  #### DATA GENERATION FOR RANDOM TRAPEZOIDAL ADJUSTMENT
  ## TOTAL NON-PARTICIPATING CASES & CONTROLS BY EXPOSURE STATUS
  npD1E1 <- D1E1-pD1E1
  npD1E0 <- D1E0-pD1E0
  npD0E1 <- D0E1-pD0E1
  npD0E0 <- D0E0-pD0E0
  
  ## TOTAL NON-PARTICIPATING CASES & CONTROLS
  npD1 <- npD1E1 + npD1E0
  npD0 <- npD0E1 + npD0E0
  
  ## EXPOSURE PREVALENCE NON-PARTICIPATING CASES & CONTROLS
  prevE_np_D1 <- npD1E1/npD1
  prevE_np_D0 <- npD0E1/npD0
  
  ## EXPOSURE PREVALENCE PARTICIPATING CASES & CONTROLS
  prevE_p_D1 <- pD1E1/(pD1E1 + pD1E0)
  prevE_p_D0 <- pD0E1/(pD0E1 + pD0E0)
  
  ## SET UP PARALLEL PROCESSING
  plan(multisession, workers = parallel::detectCores() - 1)
  
  #### BIAS-ADJUSTED ESTIMATE
  CL95_list_a <- future_lapply(1:niter, function(i) {
    
    #### MODES FOR TRAPEZOIDAL
    mode1 <- ifelse(prevE_p_D1[i] < prevE_p_D0[i], prevE_p_D1[i], prevE_p_D0[i])
    mode2 <- ifelse(prevE_p_D0[i] > prevE_p_D1[i], prevE_p_D0[i], prevE_p_D1[i])
    
    #### EXPOSURE PREVALENCE - NON-PARTICIPATING CASES - CENTERED ON PREVE_P_D1 & PREVE_P_D0
    prevE_np_D1_tr <- rtrapezoid(nsubiter, min = 0, mode1 = mode1, mode2 = mode2, max = 1)
    prevE_np_D0_tr <- rtrapezoid(nsubiter, min = 0, mode1 = mode1, mode2 = mode2, max = 1)
    
    npD1E1_a <- rbinom(nsubiter, npD1[i], prevE_np_D1_tr) # NONPARTICIPATING EXPOSED CASES
    npD1E0_a <- npD1[i] - npD1E1_a # NONPARTICIPATING UNEXPOSED CASES
    npD0E1_a <- rbinom(nsubiter, npD0[i], prevE_np_D0_tr) # NONPARTICIPATING EXPOSED CONTROLS
    npD0E0_a <- npD0[i] - npD0E1_a # NONPARTICIPATING UNEXPOSED CONTROLS
    
    #### CALCULATING BIAS-ADJUSTED OR FOR RANDOM TRAPEZOIDAL
    OR_a_vec <- ((pD1E1[i] + npD1E1_a) / (pD1E0[i] + npD1E0_a)) / 
      ((pD0E1[i] + npD0E1_a) / (pD0E0[i] + npD0E0_a))
    L_SERR_a <- sqrt((1/(pD1E1[i])) + (1/(pD1E0[i])) + 
                        (1/(pD0E1[i])) + (1/(pD0E0[i])))
    OR_a_vec <- exp(log(OR_a_vec) + rnorm(nsubiter) * L_SERR_a)
    
    ### CALCULATE 95% SI
    quantiles1 <- quantile(OR_a_vec, probs = c(0.025, 0.5, 0.975))
    c(quantiles1, ifelse(quantiles1[1] <= RR & quantiles1[3] >= RR, 1, 0))
    
  }, future.seed = TRUE)
  
  CL95_a <- do.call(rbind, CL95_list_a)
  
  ## SAVING RESULTS
  res <- data.frame(RR1, OR, LCL95, UCL95, MSE, coverage, OR_b, LCL95_b, UCL95_b, MSE_b, coverage_b, prevE_D1, prevE_D0, prevE_p_D1, prevE_p_D0, prevE_np_D1, prevE_np_D0)
  res_a <- as.data.frame(CL95_a)
  colnames(res_a) <- c("LCL95_a", "median_OR_a", "UCL95_a","coverage_a")
  
  return(list(res = res, res_a = res_a))
}

#### SIMULATION SCENARIO ONE:
set.seed(seed = 23) 
sim1.3 <- SIM_FUN_TRAP(niter = 100000, ssize = 40000, RR = 2, prevE = 0.2, baserisk = 0.01, prD1E1 = 0.75, prD1E0 = 0.95, prD0E1 = 0.82, prD0E0 = 0.88, nsubiter = 10000)

sim1.3_res <- sim1.3$res
sim1.3_res_a <- sim1.3$res_a

## true OR
quantile(sim1.3_res$OR, probs=c(0.025, 0.5, 0.975))
median(sim1.3_res$MSE)
sum(sim1.3_res$coverage)/length(sim1.3_res$coverage)

## biased OR
quantile(sim1.3_res$OR_b, probs=c(0.025, 0.5, 0.975))
median(sim1.3_res$MSE_b)
sum(sim1.3_res$coverage_b)/length(sim1.3_res$coverage_b)

## bias-adjusted OR
quantile(sim1.3_res_a$median_OR_a, probs=c(0.025, 0.5, 0.975))
sim1.3_MSE_a <- (mean((log(sim1.3_res_a$median_OR_a) - log(true_RR))^2))
sim1.3_MSE_a
sum(sim1.3_res_a$coverage_a)/length(sim1.3_res_a$coverage_a)

#### SIMULATION SCENARIO ONE:
set.seed(seed = 23)
sim2.3 <- SIM_FUN_TRAP(niter = 100000, ssize = 40000, RR = 2, prevE = 0.2, baserisk = 0.01 , prD1E1 = 0.82, prD1E0 = 0.88, prD0E1 = 0.75, prD0E0 = 0.95, nsubiter = 10000)

sim2.3_res <- sim2.3$res
sim2.3_res_a <- sim2.3$res_a

## true OR
quantile(sim2.3_res$OR, probs=c(0.025, 0.5, 0.975))
median(sim2.3_res$MSE)
sum(sim2.3_res$coverage)/length(sim2.3_res$coverage)

## biased OR
quantile(sim2.3_res$OR_b, probs=c(0.025, 0.5, 0.975))
median(sim2.3_res$MSE_b)
sum(sim2.3_res$coverage_b)/length(sim2.3_res$coverage_b)

## bias-adjusted OR
quantile(sim2.3_res_a$median_OR_a, probs=c(0.025, 0.5, 0.975))
sim2.3_MSE_a <- mean((log(sim2.3_res_a$median_OR_a) - log(true_RR))^2)
sim2.3_MSE_a
sum(sim2.3_res_a$coverage_a)/length(sim2.3_res_a$coverage_a)

##### FOURTH SIMULATION, BETA DISTRIBUTION ############
## Writing a function (SIM_FUN_BETA) for the simulation
SIM_FUN_BETA <- function(niter, ssize, RR, prevE, baserisk, prD1E1, prD1E0, prD0E1, prD0E0, nsubiter) {
  
  #### COHORT SIMULATION DATA GENERATION
  E1 <- rbinom(niter, ssize, prevE) # total exposed
  E0 <- ssize - E1 # total unexposed
  D1E1 <- rbinom(niter, E1, RR*baserisk) # exposed cases
  D1E0 <- rbinom(niter, E0, baserisk) # unexposed cases
  
  #### CALCULATING TRUE RR ESTIMATE FOR COHORT
  RR1 <- (D1E1/E1)/(D1E0/E0)
  
  #### DATA GENERATION FOR CASE-COHORT ESTIMATE
  D0 <- (D1E1 + D1E0)*2 # 2 controls per case
  D0E1 <- rbinom(niter, D0, prevE) # exposed controls
  D0E0 <- D0 - D0E1 # unexposed controls
  
  #### CALCULATING OR ESTIMATE FOR CASE-COHORT
  OR <- (D1E1/D0E1)/(D1E0/D0E0)
  L_SERR <- sqrt(1/D1E1+1/D0E1+1/D1E0+1/D0E0)
  ME_OR <- qnorm(0.975)*L_SERR
  LCL95 <- OR*exp(-ME_OR)
  UCL95 <- OR*exp(ME_OR)
  MSE <- mean((log(OR) - log(RR))^2)
  coverage <- ifelse(LCL95 <= RR & UCL95 >= RR, 1, 0)
  
  #### EXPOSURE PREVALENCE CASES AND CONTROLS
  prevE_D1 <- D1E1/(D1E1+D1E0)
  prevE_D0 <- D0E1/(D0E1+D0E0)
  
  #### DATA GENERATION FOR BIASED ESTIMATE
  pD1E1 <- rbinom(niter, D1E1, prD1E1)
  pD1E0 <- rbinom(niter, D1E0, prD1E0)
  pD0E1 <- rbinom(niter, D0E1, prD0E1)
  pD0E0 <- rbinom(niter, D0E0, prD0E0)
  
  #### CALCULATING BIASED OR ESTIMATE
  OR_b <- (pD1E1/pD0E1)/(pD1E0/pD0E0) 
  median(OR_b)
  L_SERR_b <- sqrt(1/pD1E1+1/pD0E1+1/pD1E0+1/pD0E0)
  ME_OR_b <- qnorm(0.975)*L_SERR_b
  LCL95_b <- OR_b*exp(-ME_OR_b)
  UCL95_b <- OR_b*exp(ME_OR_b)
  MSE_b <- mean((log(OR_b) - log(RR))^2)
  coverage_b <- ifelse(LCL95_b <= RR & UCL95_b >= RR, 1, 0)
  
  #### DATA GENERATION FOR BETA ADJUSTMENT
  ## TOTAL NON-PARTICIPATING CASES & CONTROLS BY EXPOSURE STATUS
  npD1E1 <- D1E1-pD1E1
  npD1E0 <- D1E0-pD1E0
  npD0E1 <- D0E1-pD0E1
  npD0E0 <- D0E0-pD0E0
  
  ## TOTAL NON-PARTICIPATING CASES & CONTROLS
  npD1 <- npD1E1 + npD1E0
  npD0 <- npD0E1 + npD0E0
  
  ## EXPOSURE PREVALENCE NON-PARTICIPATING CASES & CONTROLS
  prevE_np_D1 <- npD1E1/npD1
  prevE_np_D0 <- npD0E1/npD0
  
  ## EXPOSURE PREVALENCE PARTICIPATING CASES & CONTROLS
  prevE_p_D1 <- pD1E1/(pD1E1 + pD1E0)
  prevE_p_D0 <- pD0E1/(pD0E1 + pD0E0)
  
  ## SET UP PARALLEL PROCESSING
  plan(multisession, workers = parallel::detectCores() - 1)
  
  #### BIAS-ADJUSTED ESTIMATE
  CL95_list_a <- future_lapply(1:niter, function(i) {
    
    #### EXPOSURE PREVALENCE - NON-PARTICIPATING CASES - CENTERED ON PD1E1/PD0E0 & PD0E1/PD0E0
    prevE_np_D1_b <- rbeta(nsubiter, pD1E1[i]+1, pD1E0[i]+1)
    prevE_np_D0_b <- rbeta(nsubiter, pD0E1[i]+1, pD0E0[i]+1)
    
    npD1E1_a <- rbinom(nsubiter, npD1[i], prevE_np_D1_b) # NONPARTICIPATING EXPOSED CASES
    npD1E0_a <- npD1[i] - npD1E1_a # NONPARTICIPATING UNEXPOSED CASES
    npD0E1_a <- rbinom(nsubiter, npD0[i], prevE_np_D0_b) # NONPARTICIPATING EXPOSED CONTROLS
    npD0E0_a <- npD0[i] - npD0E1_a # NONPARTICIPATING UNEXPOSED CONTROLS
    
    #### CALCULATING BIAS-ADJUSTED OR FOR RANDOM TRAPEZOIDAL
    OR_a_vec <- ((pD1E1[i] + npD1E1_a) / (pD1E0[i] + npD1E0_a)) / 
      ((pD0E1[i] + npD0E1_a) / (pD0E0[i] + npD0E0_a))
    L_SERR_a <- sqrt((1/(pD1E1[i]))+(1/(pD1E0[i])) + 
                        (1/(pD0E1[i])) + (1/(pD0E0[i])))
    OR_a_vec <- exp(log(OR_a_vec) + rnorm(nsubiter) * L_SERR_a)
    
    ### CALCULATE 95% SI
    quantiles1 <- quantile(OR_a_vec, probs=c(0.025, 0.5, 0.975))
    c(quantiles1, ifelse(quantiles1[1] <= RR & quantiles1[3] >= RR, 1, 0))
    
  }, future.seed = TRUE)
  
  CL95_a <- do.call(rbind, CL95_list_a)
  
  ## SAVING RESULTS
  res <- data.frame(RR1, OR, LCL95, UCL95, MSE, coverage, OR_b, LCL95_b, UCL95_b, MSE_b, coverage_b, prevE_D1, prevE_D0, prevE_p_D1, prevE_p_D0, prevE_np_D1, prevE_np_D0)
  res_a <- as.data.frame(CL95_a)
  colnames(res_a) <- c("LCL95_a", "median_OR_a", "UCL95_a","coverage_a")
  
  return(list(res = res, res_a = res_a))
}

#### SIMULATION SCENARIO ONE:
set.seed(seed = 23) 
sim1.4 <- SIM_FUN_BETA(niter = 100000, ssize = 40000, RR = 2, prevE = 0.2, baserisk = 0.01, prD1E1 = 0.75, prD1E0 = 0.95, prD0E1 = 0.82, prD0E0 = 0.88, nsubiter = 10000)

sim1.4_res <- sim1.4$res
sim1.4_res_a <- sim1.4$res_a

## true OR
quantile(sim1.4_res$OR, probs=c(0.025, 0.5, 0.975))
median(sim1.4_res$MSE)
sum(sim1.4_res$coverage)/length(sim1.4_res$coverage)

## biased OR
quantile(sim1.4_res$OR_b, probs=c(0.025, 0.5, 0.975))
median(sim1.4_res$MSE_b)
sum(sim1.4_res$coverage_b)/length(sim1.4_res$coverage_b)

## bias-adjusted OR
quantile(sim1.4_res_a$median_OR_a, probs=c(0.025, 0.5, 0.975))
sim1.4_MSE_a <- (mean((log(sim1.4_res_a$median_OR_a) - log(true_RR))^2))
sim1.4_MSE_a
sum(sim1.4_res_a$coverage_a)/length(sim1.4_res_a$coverage_a)

#### SIMULATION SCENARIO ONE:
set.seed(seed = 23)
sim2.4 <- SIM_FUN_BETA(niter = 100000, ssize = 40000, RR = 2, prevE = 0.2, baserisk = 0.01 , prD1E1 = 0.82, prD1E0 = 0.88, prD0E1 = 0.75, prD0E0 = 0.95, nsubiter = 10000)

sim2.4_res <- sim2.4$res
sim2.4_res_a <- sim2.4$res_a

## true OR
quantile(sim2.4_res$OR, probs=c(0.025, 0.5, 0.975))
median(sim2.4_res$MSE)
sum(sim2.4_res$coverage)/length(sim2.4_res$coverage)

## biased OR
quantile(sim2.4_res$OR_b, probs=c(0.025, 0.5, 0.975))
median(sim2.4_res$MSE_b)
sum(sim2.4_res$coverage_b)/length(sim2.4_res$coverage_b)

## beta (centered on npD1E1/npD1E0 & npD0E1/npD0E0)
quantile(sim2.4_res_a$median_OR_a, probs=c(0.025, 0.5, 0.975))
sim2.4_MSE_a <- mean((log(sim2.4_res_a$median_OR_a) - log(true_RR))^2)
sim2.4_MSE_a
sum(sim2.4_res_a$coverage_a)/length(sim2.4_res_a$coverage_a)