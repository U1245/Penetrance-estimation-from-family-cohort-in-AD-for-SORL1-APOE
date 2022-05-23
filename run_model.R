

run.SORL1_APOE_model <- function(dat, timeName, statusName, cuts, time_censoring = NA, maxIter = 30, cv.threshold = 1e-3, uAPOE) {
  
  ### Initialisations
  
  # number of individuals
  n <- dim(dat)[1]
  
  dat$status <- dat[, statusName]
  dat$time1 <- ifelse(is.na(as.numeric(dat[, timeName])) | as.numeric(dat[, timeName]) <= 0, 0.1, as.numeric(dat[, timeName])) # time should not equal 0
  dat$time <- dat$time1

  # Take into account ascertainment by removing index information
  dat$status <- ifelse(dat$index == 1, 0, dat$status)
  dat$time <- ifelse(dat$index == 1, 0.1, dat$time)

  # Censoring end time
  if (!is.na(time_censoring)) {
    dat$status <- ifelse(dat$time > time_censoring, 0, dat$status)
    dat$time <- ifelse(dat$time > time_censoring, time_censoring, dat$time)
    cat("Time censored at", time_censoring, "years old \n")
  } else{
    cat("Time not censored \n")
  }
  
  # beta names
  beta_names <- paste0("beta", 1:length(cuts))
  
  ### APOE
  
  source("Lambda_and_lambda_APOE.R")
  
  # Compute CumulL_nc (cumulative hazard for non carriers) at cut times (according to APOE genotype)
  CumulL_nc.i.cuts <- list()
  for (j in 1:length(cuts)) {
    CumulL_nc.i.0 <- mapply(CumulLambda.APOE, rep(0, n), cuts[j], u = uAPOE); CumulL_nc.i.1 <- mapply(CumulLambda.APOE, rep(1, n), cuts[j], u = uAPOE); CumulL_nc.i.2 <- mapply(CumulLambda.APOE, rep(2, n), cuts[j], u = uAPOE)
    CumulL_nc.i.cuts[[j]] <- matrix(rep(c(CumulL_nc.i.0, CumulL_nc.i.0, CumulL_nc.i.1, CumulL_nc.i.0, CumulL_nc.i.0, CumulL_nc.i.1, CumulL_nc.i.1, CumulL_nc.i.1, CumulL_nc.i.2), 4), ncol = 36)
  }
  # Compute l_nc (hazard for non carriers) at cut times (according to APOE genotype)
  l_nc.i.cuts <- list()
  for (j in 1:length(cuts)) {
    l_nc.i.0 <- mapply(lambda.APOE, rep(0, n), cuts[j], u = uAPOE); l_nc.i.1 <- mapply(lambda.APOE, rep(1, n), cuts[j], u = uAPOE); l_nc.i.2 <- mapply(lambda.APOE, rep(2, n), cuts[j], u = uAPOE)
    l_nc.i.cuts[[j]] <- matrix(rep(c(l_nc.i.0, l_nc.i.0, l_nc.i.1, l_nc.i.0, l_nc.i.0, l_nc.i.1, l_nc.i.1, l_nc.i.1, l_nc.i.2), 4), ncol = 36)
  }
  
  # Compute CumulL_nc (cumulative hazard for non carriers) for each APOE genotype (N x 36 matrix) at observed time
  CumulL_nc.i.0 <- mapply(CumulLambda.APOE, rep(0, n), dat$time, u = uAPOE); CumulL_nc.i.1 <- mapply(CumulLambda.APOE, rep(1, n), dat$time, u = uAPOE); CumulL_nc.i.2 <- mapply(CumulLambda.APOE, rep(2, n), dat$time, u = uAPOE)
  CumulL_nc.i <- matrix(rep(c(CumulL_nc.i.0, CumulL_nc.i.0, CumulL_nc.i.1, CumulL_nc.i.0, CumulL_nc.i.0, CumulL_nc.i.1, CumulL_nc.i.1, CumulL_nc.i.1, CumulL_nc.i.2), 4), ncol = 36)
  
  # Compute l_nc (hazard for non carriers) for each APOE genotype (N x 36 matrix) at observed time
  l_nc.i.0 <- mapply(lambda.APOE, rep(0, n), dat$time, u = uAPOE); l_nc.i.1 <- mapply(lambda.APOE, rep(1, n), dat$time, u = uAPOE); l_nc.i.2 <- mapply(lambda.APOE, rep(2, n), dat$time, u = uAPOE)
  l_nc.i <- matrix(rep(c(l_nc.i.0, l_nc.i.0, l_nc.i.1, l_nc.i.0, l_nc.i.0, l_nc.i.1, l_nc.i.1, l_nc.i.1, l_nc.i.2), 4), ncol = 36)
  
  # n x 36 matrix of status (delta_i)
  delta.i <- matrix(rep(dat$status, 36), ncol = 36)
  
  ### EM algorithm
  set.seed(123)
  # weights initialisation 
  w.ig <- matrix(runif(n*36, 0, 1), ncol = 36)
  x <- outer(c("w.22", "w.32", "w.42", "w.23", "w.33", "w.43", "w.24", "w.34", "w.44"), c(".00", ".01", ".10", ".11"), FUN = "paste0")
  dim(x) <- NULL; colnames(w.ig) <- x
  w.sum <- w.ig
  
  # beta initialisation
  beta_old <- rep(0, length(beta_names))
  
  # convergence initialisation
  convergence <- FALSE
  nbIter <- 0
  
  # Alternate M- and E-steps
  while (!convergence & nbIter < maxIter) {
    # M-step => compute beta
    beta <- rep(NA, length(beta_names))
    J <- length(beta)
    cuts.inf <- c(cuts[-1], Inf)
    cuts.0 <- cuts
    for (j in 1:J) {
      beta.num <- 0
      beta.den <- 0
      intervalle.sup <- cuts.inf[j]
      intervalle.min <- cuts.0[j]
      beta.num <- sum((w.sum*delta.i)[dat$time >= intervalle.min & dat$time < intervalle.sup, 10:36])
      beta.den <- sum((w.sum*(CumulL_nc.i - CumulL_nc.i.cuts[[which(cuts == intervalle.min)]]))[dat$time >= intervalle.min & dat$time < intervalle.sup, 10:36])
      if (intervalle.sup != Inf) {
        beta.den <- beta.den + sum((w.sum*(CumulL_nc.i.cuts[[which(cuts == intervalle.sup)]] - CumulL_nc.i.cuts[[which(cuts == intervalle.min)]]))[dat$time >= intervalle.sup, 10:36])
      }
      beta[j] <- log(beta.num/beta.den)
      
      if (beta[j] == -Inf) {
        beta[j] <- 0
      } 
    }
    # print current beta value
    cat(beta, "\n")
    
    if (sum(!(abs(beta_old - beta) < cv.threshold)) == 0) {
      convergence <- TRUE
    }
    beta_old <- beta
    
    # E-step => compute weights
    
    # compute evidence
    evidence <- matrix(NA, nrow = n, ncol = 36)
    
    # SORL1 WT
    evidence[dat$status == 1, 1:9] = (exp(-CumulL_nc.i) * l_nc.i)[dat$status == 1, 1:9]
    evidence[dat$status == 0, 1:9] = (exp(-CumulL_nc.i))[dat$status == 0, 1:9]
    
    # SORL1 +
    L_c.i.inf <- 0
    for (j in 1:length(cuts)) {
      cut.inf <- cuts[j]
      cut.sup <- cuts.inf[j]
      L_c.i <- L_c.i.inf + (CumulL_nc.i - CumulL_nc.i.cuts[[j]])*exp(beta[beta_names == beta_names[j]])
      evidence[dat$status == 0 & dat$time >= cut.inf & dat$time < cut.sup, 10:36] = (exp(-L_c.i))[dat$status == 0 & dat$time >= cut.inf & dat$time < cut.sup, 10:36]
      evidence[dat$status == 1 & dat$time >= cut.inf & dat$time < cut.sup, 10:36] = (exp(-L_c.i)*l_nc.i*exp(beta[beta_names == beta_names[j]]))[dat$status == 1 & dat$time >= cut.inf & dat$time < cut.sup, 10:36]
      if (cut.sup != Inf) {
        L_c.i.inf <- L_c.i.inf + (CumulL_nc.i.cuts[[j + 1]] - CumulL_nc.i.cuts[[j]])*exp(beta[beta_names == beta_names[j]])
      }
    }
    
    # Replace evidence by 0 when we know the truth
    # APOE
    evidence[dat$APOE == 22, c(2:9, 11:18, 20:27, 29:36)] <- 0
    evidence[dat$APOE %in% c(24, 42), c(1, 2, 4, 5, 6, 8:11, 13, 14, 15, 17:20, 22, 23, 24, 26:29, 31, 32, 33, 35, 36)] <- 0
    evidence[dat$APOE %in% c(34, 43), c(1:5, 7, 9:14, 16, 18:23, 25, 27:32, 34, 36)] <- 0
    evidence[dat$APOE %in% c(23, 32), c(1, 3, 5:10, 12, 14:19, 21, 23:28, 30, 32:36)] <- 0
    evidence[dat$APOE == 44, c(1:8, 10:17, 19:26, 28:35)] <- 0
    evidence[dat$APOE == 33, c(1:4, 6:13, 15:22, 24:31, 33:36)] <- 0
    # SORL1 known
    evidence[dat$SORL1 %in% c("01", "10"), c(1:9, 28:36)] <- 0
    evidence[dat$SORL1 %in% c("11"), c(1:27)] <- 0
    evidence[dat$SORL1 == "00", 10:36] <- 0
    
    
    ### Save data for bped algorithm
    write.table(file = "pedigree.ped", dat[,c("famid", "indid", "patid", "matid")], row.names = FALSE, col.names = FALSE, sep = "\t")
    write.table(file = "evidence.ev", cbind(dat[,c("famid", "indid")], evidence[,]), row.names = FALSE, col.names = FALSE, sep = "\t")
    
    # Run bped3alleles+2alleles
    system("/opt/Bped2/bped3alleles2alleles pedigree.ped evidence.ev 0.08452747 0.7641758 0.1512967 0.9999 0.0001 > margProb.out")
    
    marginalProbabilities <- read.delim("margProb.out", header = FALSE, stringsAsFactors = FALSE)
    mp <- marginalProbabilities[, 2:37]
    x <- outer(c("mp.22", "mp.32", "mp.42", "mp.23", "mp.33", "mp.43", "mp.24", "mp.34", "mp.44"), c(".00", ".01", ".10", ".11"), FUN = "paste0")
    dim(x) <- NULL
    colnames(mp) <- x
    
    # Update weights
    w.ig <- mp
    w.sum = w.ig
    
    # print(head(w.ig))
    nbIter <- nbIter + 1
  }
  
  names(beta) <- beta_names
  
  # Save weights
  dat.weights <- cbind(dat, mp)
  
  
  return(list(beta = beta, convergence = convergence, weights = dat.weights))
}
