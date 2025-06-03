#' Freedman-Lane Permutation Test for Partial Correlation Networks
#' @param data_list List of data matrices
#' @param pcor_matrix Initial partial correlation matrix 
#' @param n_perm Number of permutations
#' @param alpha Significance level
#' @return Permutation test threshold
source("Packages/stabENG.r")
source("Packages/MyENG.r")
permtest_freedman_lane <- function(data_list, 
                                 pcor_matrix,
                                 group_name,
                                 n_perm = 50,
                                 alpha = 0.1) {
  
  # 1. Fit full model using all covariates
  base_model <- preprocess_and_estimate_network(
    data_list = data_list,
    labels = colnames(pcor_matrix),
    lambda1 = stabENG_params$opt.lambda1,
    lambda2 = stabENG_params$opt.lambda2
  )
  
  # 2. Get residuals
  residuals <- lapply(data_list, function(d) {
    lm(d ~ 1)$residuals # 可根据需要加入其他协变量
  })
  
  # 3. Freedman-Lane 置换过程
  null_dist <- numeric(0)
  for(i in 1:n_perm) {
    # 3.1 置换残差
    perm_residuals <- lapply(residuals, function(r) {
      r[sample(nrow(r)),] 
    })
    
    # 3.2 重构数据  
    perm_data <- Map(function(res, orig) {
      fitted <- lm(orig ~ 1)$fitted
      perm <- fitted + res
      return(perm)
    }, perm_residuals, data_list)
    
    # 3.3 计算置换样本的偏相关
    perm_pcor <- preprocess_and_estimate_network(
      perm_data,
      labels = colnames(pcor_matrix),
      lambda1 = stabENG_params$opt.lambda1,
      lambda2 = stabENG_params$opt.lambda2
    )$pcor[[group_name]]
    
    null_dist <- c(null_dist, abs(perm_pcor[upper.tri(perm_pcor)])) 
  }
  
  # 4. 计算阈值
  threshold <- quantile(null_dist, 1-alpha)
  
  return(threshold)
}