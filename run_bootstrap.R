

run.SORL1_APOE_bootstrap <- function(dat, timeName, statusName, cuts, MAF, proband_info, time_censoring = NA, maxIter = 30, cv.threshold = 1e-3)
{

  library(foreach)

  datAll <- dat
  allFamilies <- unique(datAll$famID)
  
  resultsBootstrap <- foreach(iter = 1:B, .combine = "rbind", .multicombine=TRUE) %do% {
    
    set.seed(iter)
    # u.iter is a parameter to add the variability of the APOE estimation in the bootstrap analysis
    u.iter <- rnorm(1, 0, 1) 
    # select data for bootstrap iteration
    fam.boo <- sample(allFamilies, length(allFamilies), replace = TRUE)
    # attribute uniuqe number ID to each family
    datNew <- subset(datAll, famID == fam.boo[1])
    datNew$famID.num <- 1
    for (k in 2:length(fam.boo)) {
      datNew.k <- subset(datAll, famID == fam.boo[k])
      datNew.k$famID.num <- k
      datNew <- rbind(datNew, datNew.k)
    }
    dat <- datNew
 
    source("run_model.R")
    res.iter <- run.SORL1_APOE_model(dat = dat, timeName = timeName, statusName = statusName, cuts = cuts, time_censoring = time_censoring, 
                                     uAPOE = u.iter, maxIter = maxIter, cv.threshold = cv.threshold)

    
    Lambda <- res.iter[["Lambda"]]
    
    ### Saving the penetrance into table
    age <- 40:100
    
    penetrance0 <- 1 - exp(-mapply(Lambda, age, 0)); names(penetrance0) <- paste0("t", age)
    penetrance1 <- 1 - exp(-mapply(Lambda, age, 1)); names(penetrance1) <- paste0("t", age)
    penetrance2 <- 1 - exp(-mapply(Lambda, age, 2)); names(penetrance2) <- paste0("t", age)
    
    beta_cv <- c(res.iter[["beta"]], res.iter[["convergence"]])
    names(beta_cv) <- c(names(res.iter[["beta"]]), "convergence")
    return(list(beta_cv, penetrance0, penetrance1, penetrance2))
  }
  penetrance <- list(t(data.frame(resultsBootstrap[,2])), t(data.frame(resultsBootstrap[,3])), t(data.frame(resultsBootstrap[,4])))

  return(penetrance)
}