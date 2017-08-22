idingo <- function(dat,dat2=NULL,dat3=NULL,x,plats=NULL,rhoarray=NULL,diff.score=T,B=100,verbose=T,cores=1) {
#  Input:
####  - dat : n by p response data
####  - dat2 : second n by p response data
####  - dat3 : third n by p response data
####  - x : length n covariate vector
####  - plats : length r (number of response data sets) vector of platform names
####  - rhoarray : Candidate tuning parameters of glasso for fitting global network. If it is one value, then we use the value. If it is null, we set the candidates.
####  - diff.score : if TRUE, edge-wise differential scores are calculated from bootstrap standard error.
####  - B : the number of bootstrap samples. It is used if diff.score is TRUE
####  - verbose : if TRUE, the number of bootstrap samples are listed
####  - cores : the number of cores to use for parallel bootstrapping

  # Append platform names to gene symbols.
  if(is.null(plats)) plats <- c("plat1", "plat2", "plat3")
  colnames(dat) <- paste(plats[1], colnames(dat), sep = "||")
  
  if (is.null(dat2) & is.null(dat3)){

    # If only one data set is provided, run standard DINGO.
    return(dingo(dat,x,rhoarray,diff.score,B,verbose,cores))

  } else {

    # Append platform names to gene symbols.
    colnames(dat2) <- paste(plats[2], colnames(dat2), sep = "||")

    # Run DINGO for 1st platform.
    fit1 <- dingo(dat,x,rhoarray,diff.score,B,verbose,cores)
    if (verbose) cat("Model 1 complete","\n")

    # Run DINGO for 1st & 2nd platform.
    fit2 <- dingo(cbind(dat, dat2),x,rhoarray,diff.score,B,verbose,cores)
    if (verbose) cat("Model 2 complete","\n")

    # Remove platform1-platform1 edges from 2nd network.
    dups <- paste(fit2$genepair$gene1, fit2$genepair$gene2, sep = "-") %in%
      paste(fit1$genepair$gene1, fit1$genepair$gene2, sep = "-")
    fit2$genepair <- fit2$genepair[!dups,]
    fit2$R1 <- fit2$R1[!dups]
    fit2$R2 <- fit2$R2[!dups]
    fit2$diff.score <- fit2$diff.score[!dups]
    # Recalculate p-value based on updated diff.score distribution
    fit2$p.val <- getEfronp(fit2$diff.score, plotIt = FALSE)$correctp

    # Merge DINGO objects
    fit.merged <- list(genepair = rbind(fit1$genepair, fit2$genepair),
                       R1 = c(fit1$R1, fit2$R1),
                       R2 = c(fit1$R2, fit2$R2),
                       diff.score = c(fit1$diff.score, fit2$diff.score),
                       p.val = c(fit1$p.val, fit2$p.val))

    if (!is.null(dat3)) {

      # Append platform names to gene symbols for 3rd platform set.
      colnames(dat3) <- paste(plats[3], colnames(dat3), sep = "||")

      # Run DINGO for 1st, 2nd & 3rd platform.
      fit3 <- dingo(cbind(dat, dat2, dat3),x,rhoarray,diff.score,B,verbose,cores)
      if (verbose) cat("Model 3 complete","\n")

      # Remove plat1-plat1, plat1-plat2, and plat2-plat2 edges from 3rd network.
      dups <- paste(fit3$genepair$gene1, fit3$genepair$gene2, sep = "-") %in%
        paste(fit.merged$genepair$gene1, fit.merged$genepair$gene2, sep = "-")
      fit3$genepair <- fit3$genepair[!dups,]
      fit3$R1 <- fit3$R1[!dups]
      fit3$R2 <- fit3$R2[!dups]
      fit3$diff.score <- fit3$diff.score[!dups]
      # Recalculate p-value based on updated diff.score distribution
      fit3$p.val <- getEfronp(fit3$diff.score, plotIt = FALSE)$correctp

      # Merge the three networks.
      fit.merged <- list(genepair = rbind(fit.merged$genepair, fit3$genepair),
                         R1 = c(fit.merged$R1, fit3$R1),
                         R2 = c(fit.merged$R2, fit3$R2),
                         diff.score = c(fit.merged$diff.score, fit3$diff.score),
                         p.val = c(fit.merged$p.val, fit3$p.val))

    }

    return(list(genepair = fit.merged$genepair, levels.x = fit2$levels.x, R1 = fit.merged$R1, R2 = fit.merged$R2, diff.score = fit.merged$diff.score, p.val = fit.merged$p.val))

  }

}
