library(foreach)
library(doParallel)
library(parallel)
library(mixedCCA)
library(SPRING)
library(stabJGL)
stabENG = function(Y,var.thresh, subsample.ratio = NULL, labels = NULL,
                   rep.num,  nlambda1,lambda1.min,lambda1.max,nlambda2,lambda2.min,lambda2.max,lambda2.init,
                   ebic.gamma,verbose=F, retune.lambda1=F,parallelize=T,nCores,weights="equal"){

  if (length(unique(unlist(lapply(Y, ncol)))) !=1) {
    stop("dimension (no of variables) not equal for all data sets.")
  }
  if (nlambda1<2) {
    stop("more values of lambda1 should be considered. Try nlambda1=20. \n")
  }
  if (lambda1.min > lambda1.max) {
    stop("lambda1.min must be smaller than lambda1.max. \n")
  }
  if (nlambda2<2) {
    stop("more values of lambda2 should be considered. Try nlambda2=20. \n")
  }
  if (lambda2.min > lambda2.max) {
    stop("lambda2.min must be smaller than lambda2.max. \n")
  }
  if(!is.null(subsample.ratio)){
    if(subsample.ratio < 0 | subsample.ratio > 1) stop("subsample ratio must be between 0 and 1. \n")
  }
  if (var.thresh < 0 | var.thresh > 1) {
    stop("variability threshold should be between 0 and 1. \n")
  }
  if (ebic.gamma < 0) {
    stop("ebic.gamma cannot have a negative value. \n")
  }
  if (parallelize) {
    if (missing(nCores)) {
      nCores = parallel::detectCores() - 2
    } else if (nCores < 2) stop("if method is to be run in parallel, at least two threads should be initiated. Try nCores=2. \n")
    cat('Parallelizing with ', nCores, ' cores. \n')
  }
  if (rep.num < 1) {
    stop("Number of subsamplings must be positive. \n")
  }

  est = list()
  est.lambda1 = stabENG_select_lambda1(Y,labels=labels, weights=weights, stars.thresh=var.thresh,
                                   stars.subsample.ratio=subsample.ratio, rep.num=rep.num,nlambda1=nlambda1,lambda1.min=lambda1.min,
                                   lambda1.max=lambda1.max, lambda2=lambda2.init,verbose=verbose,parallelize=parallelize,nCores=nCores)
  est.lambda2 = stabENG_select_lambda2_eBIC(Y,labels=labels,weights=weights,
                                        nlambda2=nlambda2,lambda2.min=lambda2.min,lambda2.max=lambda2.max,
                                        lambda1=est.lambda1$opt.lambda1,gamma=ebic.gamma,verbose=verbose)
  est$opt.ebic = est.lambda2$opt.ebic
  est$ebic.vals = est.lambda2$ebic.vals
  est$lambda2s = est.lambda2$lambda2s
  est$opt.lambda2 = est.lambda2$opt.lambda2
  est$opt.fit = est.lambda2$opt.fit
  est$opt.sparsities = est.lambda2$opt.sparsities
  if(retune.lambda1){
    est.lambda1 = stabENG_select_lambda1(Y,labels=labels, weights=weights, stars.thresh=var.thresh,
                                     stars.subsample.ratio=subsample.ratio, rep.num=rep.num,nlambda1=nlambda1,lambda1.min=lambda1.min,
                                     lambda1.max=lambda1.max, lambda2=est.lambda2$opt.lambda2,verbose=verbose,parallelize=parallelize,nCores=nCores)
    est$opt.fit = est.lambda1$opt.fit
    est$opt.fit.lambda2 = est.lambda2$opt.fit
    est$opt.sparsities = est.lambda1$opt.sparsities
    est$opt.sparsities.lambda2 = est.lambda2$opt.sparsities
  }
  else {
    est$opt.fit.lambda1  = est.lambda1$opt.fit
    est$opt.sparsities.lambda1 = est.lambda1$opt.sparsities
  }
  est$lambda1s = est.lambda1$lambda1s
  est$total.variability = est.lambda1$total.variability
  est$variability = est.lambda1$variability
  est$opt.lambda1 = est.lambda1$opt.lambda1
  est$opt.fit.pcor = preprocess_and_estimate_network(Y,labels = labels, lambda1 = est.lambda1$opt.lambda1, lambda2 = est.lambda2$opt.lambda2)$pcor
  cat('Network estimation done. \n')
  return(est)
}

stabENG_select_lambda1 = function(Y, weights="equal",labels, stars.thresh, stars.subsample.ratio = NULL,rep.num,
                                  nlambda1,lambda1.min,lambda1.max, lambda2,verbose,parallelize,nCores){
  K = length(Y)
  n.vals = unlist(lapply(Y,nrow))
  p = ncol(Y[[1]])
  stars.subsample.ratios = rep(0,K)
  seeds=sample(1:1000,rep.num)
  est=list()
  lambda1s=seq(lambda1.max,lambda1.min,length.out=nlambda1)
  if(is.null(stars.subsample.ratio))
  {
    for(k in 1:K){
      if(n.vals[k]>144) stars.subsample.ratios[k] = 10*sqrt(n.vals[k])/n.vals[k]
      if(n.vals[k]<=144) stars.subsample.ratios[k] = 0.8
    }
  }
  else {
    if(length(stars.subsample.ratio)<=1){
      stars.subsample.ratios = rep(stars.subsample.ratio,K)
    }
  }
  est = list()
  est$merge=list()
  est$lambda1s=lambda1s
  for(i in 1:nlambda1)
  {
    est$merge[[i]] = array(0,dim=c(K,p,p))
  }
  if(!parallelize){
    if(verbose) cat('Tuning lambda1 unparallelized... \n ')
    for(i in 1:rep.num)
    {
      Y.sample = list()
      for(k in 1:K){
        ind.sample = sample(c(1:n.vals[k]), floor(n.vals[k]*stars.subsample.ratios[k]), replace=FALSE)
        Y.sample[[k]] = Y[[k]][ind.sample, ]
      }
      names(Y.sample) = names(Y)
      for(j in 1:nlambda1)
      {
        lambda = lambda1s[j]
        tmp = preprocess_and_estimate_network(Y.sample,labels = labels, lambda1 = lambda, lambda2 = lambda2)$prec
        for (k in 1:K){
          est$merge[[j]][k,,] = est$merge[[j]][k,,] + (tmp[[k]]!=0)
        }
      }
      rm(ind.sample,Y.sample,tmp)
      gc()

      if (verbose) {
        done <- round(100 * i / rep.num)
        done.next <- round(100 * (i + 1) / rep.num)
        if (i == rep.num | (done %% 5) == 0 & (done.next %% 5) != 0) cat('Tuning lambda1: ', done, ' % done \n')
      }
    }
  }
  else{
    if(verbose) cat('Tuning lambda1 parallelly... \n ')
    cl <- parallel::makeCluster(nCores)
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl,Sys.info())
    res.list = foreach::foreach(i=1:rep.num,
      .export = c('stabENG_select_lambda1_parallel', 'preprocess_and_estimate_network', 
      'myENG','myadmm.iters','myadmm.iters.unconnected','flsa2','flsa.general','soft',
    'penalty.as.matrix','dsgl'),
      .packages = c("stabJGL","igraph","qgraph")) %dopar% {
      stabENG_select_lambda1_parallel(Y,labels=labels, rep.num=rep.num,n.vals=n.vals,stars.subsample.ratios=stars.subsample.ratios,
                                      lambda1s=lambda1s,lambda2=lambda2,penalize.diagonal = penalize.diagonal,
                                      seed=seeds[i], array.list=est$merge, verbose=verbose)
    }
    parallel::stopCluster(cl)
    for(j in 1:length(lambda1s)){
      for(k in 1:K){
        est$merge[[j]][k,,] = Reduce('+',lapply(res.list,function(mat) mat[[j]][k,,]))
      }
    }
    rm(res.list)
  }
  est$variability = matrix(0,nlambda1,K)
  for(i in 1:nlambda1){
    for(k in 1:K){
      est$merge[[i]][k,,] = est$merge[[i]][k,,]/rep.num
      est$variability[i,k] = 4*sum(est$merge[[i]][k,,]*(1-est$merge[[i]][k,,]))/(p*(p-1))
    }
  }
  est$total.variability = rowMeans(est$variability)
  est$opt.index = max(which.max(est$total.variability >= stars.thresh)[1]-1,1)
  est$opt.lambda1 = est$lambda1s[est$opt.index]
  est$opt.fit = preprocess_and_estimate_network(Y,labels = labels, lambda1 = est$opt.lambda1, lambda2 = lambda2)$prec
  est$opt.sparsities = unlist(lapply(est$opt.fit,sparsity))
  return(est)
  if(verbose) cat('Tuning lambda1 done. \n')
}

stabENG_select_lambda1_parallel = function(Y,labels, rep.num,n.vals,stars.subsample.ratios,
                                           lambda1s,lambda2,penalize.diagonal,
                                           seed,array.list,verbose){
  set.seed(seed)
  Y.sample = list()
  K=length(Y)
  for(k in 1:K){
    ind.sample = sample(c(1:n.vals[k]), floor(n.vals[k]*stars.subsample.ratios[k]), replace=FALSE)
    Y.sample[[k]] = Y[[k]][ind.sample, ]
  }
  names(Y.sample) = names(Y)
  for(j in 1:length(lambda1s))
  {
    lambda = lambda1s[j]
    tmp = preprocess_and_estimate_network(Y.sample,labels = labels, lambda1 = lambda, lambda2 = lambda2)$prec
    for (k in 1:K){
      array.list[[j]][k,,] = (tmp[[k]]!=0)
    }
  }
  rm(ind.sample,Y.sample,tmp)
  gc()
  return(array.list)
}

stabENG_select_lambda2_eBIC = function(Y,labels, weights="equal",
                                       nlambda2,lambda2.min,lambda2.max, lambda1=NULL,gamma=NULL,verbose){
  ebic.vals = rep(0,nlambda2)
  lambda2.vals = seq(lambda2.min,lambda2.max,length.out= nlambda2)
  n.vals = unlist(lapply(Y,nrow))
  sample.cov = lapply(Y,FUN = function(s) stats::cov(s)) #? replace by Kcor
  mods.lam2=list()
  if(verbose) cat('Tuning lambda2...\n')
  for (i in 1:nlambda2){
    mods.lam2[[i]] = preprocess_and_estimate_network(Y,labels = labels, lambda1 = lambda1, lambda2 = lambda2.vals[i])$prec
    ebic.vals[i] = stabJGL::eBIC_adapted(mods.lam2[[i]],sample.cov=sample.cov,n.vals=n.vals,gamma=gamma)
    done <- round(100 * i / nlambda2)
    if(verbose) cat('Tuning lambda2: ', done, ' % done \n')
  }
  opt.ind = which.min(ebic.vals)
  res=list(opt.fit=mods.lam2[[opt.ind]], opt.lambda1 = lambda1, opt.lambda2 = lambda2.vals[opt.ind], opt.index = opt.ind,opt.ebic=ebic.vals[opt.ind],
           ebic.vals=ebic.vals, lambda2s=lambda2.vals, mods=mods.lam2, opt.sparsities = unlist(lapply(mods.lam2[[opt.ind]],sparsity)))
}

preprocess_and_estimate_network <- function(data_list, labels, lambda1, lambda2)
  {
    # Helper function for common scaling normalization
    norm_to_total <- function(x) x / sum(x)
    common_scaling <- function(data) {
      depths <- rowSums(data)
      data_normalized <- t(apply(data, 1, norm_to_total))
      common_depth <- min(depths)  # Calculate only once
      data_common_scaled <- round(data_normalized * common_depth) 
      return(data_common_scaled)
    }
    # Preprocessing steps for each data set in the list
    processed_data_list <- lapply(data_list, function(data) {
      scaled_data <- common_scaling(data)
      mclr_data <- SPRING::mclr(scaled_data) # Using SPRING's mclr
      Kcor <- mixedCCA::estimateR(mclr_data, type = "trunc", method = "approx", tol = 1e-6, verbose = FALSE)$R
      return(Kcor)
    })

    # Sample sizes for each dataset.
    n_samples <- sapply(data_list, nrow)

    # Estimate the network using the preprocessed data
    theta <- myENG(processed_data_list, n = n_samples, weights = "equal", penalty = "fused",
    lambda1 = lambda1, lambda2 = lambda2, maxiter=500, tol=1e-5, rho=1, truncate=1e-5)
    for(i in seq(length(theta$network)))
      {
        colnames(theta$network[[i]]) <- rownames(theta$network[[i]]) <-
        colnames(theta$concentrationMatrix[[i]]) <- rownames(theta$concentrationMatrix[[i]]) <-labels
      }
    out_prec <- lapply(theta$concentrationMatrix, as.matrix)
    names(out_prec) <- names(processed_data_list)
    out_pcor <- lapply(theta$network, as.matrix)
    names(out_pcor) <- names(processed_data_list)
    return(list(prec = out_prec, pcor = out_pcor))
  }
