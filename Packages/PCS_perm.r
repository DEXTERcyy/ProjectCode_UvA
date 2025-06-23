# 辅助函数：K折划分
cv.part_internal = function(n, k) {
  k_orig_for_warning <- k # 保存原始k用于警告信息

  if (n == 0) {
    return(list(trainIndices = list(), testIndices = list()))
  }

  # 调整k值以确保其有效性
  if (n == 1) {
    if (k_orig_for_warning != 1) {
      warning(paste0("For n=1, k must be 1. Adjusted k from ", k_orig_for_warning, " to 1."))
    }
    k <- 1 # 对于 n=1, k 必须是 1 (LOOCV)
  } else { # n > 1
    if (k_orig_for_warning <= 1) {
      k <- 2 # 如果k太小，默认为2折
      warning(paste0("Original k=", k_orig_for_warning, " was <=1 for n=", n, ". Adjusted k to ", k, "."))
    } else if (k_orig_for_warning > n) {
      k <- n # 如果k太大，默认为n折 (LOOCV)
      warning(paste0("Original k=", k_orig_for_warning, " > n=", n, ". Adjusted k to ", k, " (LOOCV)."))
    } else {
      k <- k_orig_for_warning # k值有效
    }
  }

  shuffled_indices = sample(n)
  fold_sizes = rep(floor(n/k), k) # k 是调整后的k
  remainder = n %% k
  if (remainder > 0) {
    fold_sizes[1:remainder] = fold_sizes[1:remainder] + 1
  }

  train_indices_list <- list()
  test_indices_list <- list()
  current_pos = 1
  for (j_fold in 1:k) { # k 是调整后的k
    fold_size_current = fold_sizes[j_fold]
    
    if(fold_size_current == 0) { 
        test_indices_list[[j_fold]] <- integer(0)
        train_indices_list[[j_fold]] <- shuffled_indices 
        next
    }
    
    end_pos = current_pos + fold_size_current - 1
    end_pos = min(end_pos, n) 

    test_sel_indices = shuffled_indices[current_pos:end_pos]
    train_sel_indices = setdiff(shuffled_indices, test_sel_indices)
    
    test_indices_list[[j_fold]] <- test_sel_indices
    train_indices_list[[j_fold]] <- train_sel_indices
    current_pos = current_pos + fold_size_current
  }
  return(list(trainIndices = train_indices_list, testIndices = test_indices_list))
}

library(glmnet)

Psi_Screen_Beta_Lasso = function(x, R_structure, tau) {
  # x 应该是已标准化的数值矩阵
  # R_structure 是偏相关矩阵，行列名与x的列名对应
  # tau 是筛选阈值

  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.matrix(R_structure)) R_structure <- as.matrix(R_structure)

  n_vars <- ncol(x)
  n_samples_current_x <- nrow(x)

  var_names <- colnames(x)
  if (is.null(var_names)) {
    var_names <- paste0("V", 1:n_vars)
    colnames(x) <- var_names
  }
  if (is.null(rownames(R_structure))) rownames(R_structure) <- var_names
  if (is.null(colnames(R_structure))) colnames(R_structure) <- var_names

  beta_matrix_final <- matrix(0, nrow = n_vars, ncol = n_vars,
                              dimnames = list(Predictors = var_names, Response = var_names))

  R_screened <- R_structure
  R_screened[abs(R_screened) <= tau] <- 0
  diag(R_screened) <- 0 # Ensure self-loops are zero for screening predictors

  for (i_response_idx in 1:n_vars) {
    response_var_name <- var_names[i_response_idx]
    # R_screened 列是响应变量，行是预测变量 (assuming R_structure[row_pred, col_resp])
    # To get predictors for a response, we look at the column corresponding to the response
    # and find non-zero rows (predictors).
    neighbor_names <- var_names[which(R_screened[, response_var_name] != 0)]

    if (length(neighbor_names) > 0 && n_samples_current_x >= 3) {
      y_response = x[, response_var_name, drop = FALSE] # Keep as a column matrix
      x_predictors = x[, neighbor_names, drop = FALSE]

      # Ensure y_response has variance
      if (stats::var(y_response, na.rm = TRUE) < .Machine$double.eps^0.5) {
          warning(paste("Response variable", response_var_name, "has zero or near-zero variance. Skipping Lasso."))
          next
      }
      # Ensure x_predictors have variance (glmnet handles this, but good to be aware)
      # col_vars_predictors <- apply(x_predictors, 2, stats::var, na.rm = TRUE)
      # if(any(col_vars_predictors < .Machine$double.eps^0.5)){
      #     warning(paste("Some predictors for response", response_var_name, "have zero variance. They will be handled by glmnet."))
      # }


      if(ncol(x_predictors) > 0){
        # nfolds must be >=3 and <= number of observations.
        current_nfolds <- if (n_samples_current_x < 3) { 0 } else { 
            if (n_samples_current_x <= length(neighbor_names)/2) { # For very small N in the fold, use LOOCV for glmnet
                n_samples_current_x
            } else { # For somewhat larger N, a fixed fold like 5 or 10 can be used
                min(n_samples_current_x, 5) # Or consider 10 if n_samples_current_x is moderately large
            }
        }
        # Ensure nfolds is at least 3 if it's > 0 (glmnet requirement)
        if (current_nfolds > 0 && current_nfolds < 3) { current_nfolds <- 3 }

        if (current_nfolds >=3 ){
            cv_fit_lasso <- tryCatch({
              # Using glmnet defaults: standardize=TRUE, intercept=TRUE
              # If x is already standardized, standardize=TRUE in glmnet is generally safe.
              glmnet::cv.glmnet(x_predictors, y_response[,1], alpha = 1, 
                                standardize = TRUE, intercept = TRUE,
                                thresh = 1e-7, # Default convergence threshold
                                nfolds = current_nfolds)
            }, error = function(e) {
                warning(paste("cv.glmnet error for response", response_var_name, ":", conditionMessage(e)))
                NULL
            })

            if(!is.null(cv_fit_lasso) && !is.null(cv_fit_lasso$lambda.min) && !is.null(cv_fit_lasso$glmnet.fit)){
              lambda_choice <- cv_fit_lasso$lambda.min
              lasso_coeffs_with_intercept <- stats::coef(cv_fit_lasso, s = lambda_choice)
              
              # coef returns a sparse matrix; convert to dense vector
              # First element is intercept, rest are coefficients for predictors
              lasso_slopes <- as.vector(lasso_coeffs_with_intercept)[-1]
              
              if (length(lasso_slopes) == length(neighbor_names)) {
                   beta_matrix_final[neighbor_names, response_var_name] <- lasso_slopes
              } else {
                  warning(paste("Mismatch in coefficient length for response", response_var_name,
                                ". Expected", length(neighbor_names), "got", length(lasso_slopes)))
              }
            }
        }
      }
    }
  }
  return(list(beta_matrix = beta_matrix_final)) # Ensure named list for consistency
}

# --- PCS 交叉验证阈值选择函数 (使用 Lasso 进行内部回归) ---
pcs_cv_threshold_stabENG_lasso_perm <- function(
  data_list_unscaled_full, 
  initial_pcor_matrix_target, 
  group_name_target, 
  other_group_name, 
  stabENG_params_list,
  otu_labels,
  fold = NULL, 
  plot_cv_curve = FALSE,
  plot_path_prefix = NULL,
  n_data_permutations = 8 
) {
  
  gc(); print(paste("Start of pcs_cv_threshold_stabENG_lasso_perm:", pryr::mem_used()))
  # --- 初始化和基本数据检查 (在重排循环外部执行一次) ---
  if (is.null(group_name_target) || !group_name_target %in% names(data_list_unscaled_full) || 
      is.null(data_list_unscaled_full[[group_name_target]])) {
    warning(paste0("Target group '", group_name_target, "' not found or is NULL. Returning default tau=0.05."))
    return(0.05)
  }
  original_data_matrix_target_unscaled <- data_list_unscaled_full[[group_name_target]]
  original_data_matrix_target_numeric <- apply(original_data_matrix_target_unscaled, 2, as.numeric)
  if(!is.matrix(original_data_matrix_target_numeric)) original_data_matrix_target_numeric <- as.matrix(original_data_matrix_target_numeric)
  
  p_vars <- ncol(original_data_matrix_target_numeric)
  n_samples_target_orig <- nrow(original_data_matrix_target_numeric) # Store original sample size
  
  # --- Determine the actual number of folds to use ---
  actual_fold_to_use <- 0 # Initialize
  if (!is.null(fold) && is.numeric(fold) && fold > 0 && is.finite(fold) && fold <= n_samples_target_orig) {
    actual_fold_to_use <- as.integer(fold)
    cat(paste0("Using user-specified fold: ", actual_fold_to_use, " for group '", group_name_target, "'.\n"))
  } else {
    if (!is.null(fold) && ( !is.numeric(fold) || fold <= 0 || !is.finite(fold) || fold > n_samples_target_orig)) {
        warning(paste0("Invalid 'fold' value (", fold, ") provided for group '", group_name_target, 
                       "'. Auto-determining fold. N=", n_samples_target_orig, ", P=", p_vars, "."))
    }
    # Auto-determine fold
    if (n_samples_target_orig == 0) {
      actual_fold_to_use <- 0 # No samples, no folds
      warning(paste0("No samples in target group '", group_name_target, "'. Using fold = 0 (no CV). Returning default tau=0.05.\n"))
    } else if (n_samples_target_orig < (p_vars / 2)) {
      actual_fold_to_use <- n_samples_target_orig # LOOCV
      cat(paste0("Fold not specified or invalid. n_samples (", n_samples_target_orig, ") < p_vars/2 (", round(p_vars/2,1), 
                 "). Using LOOCV (fold = ", actual_fold_to_use, ") for group '", group_name_target, "'.\n"))
    } else {
      actual_fold_to_use <- 5
      # Ensure 5-fold is not more than N
      if (actual_fold_to_use > n_samples_target_orig) {
        actual_fold_to_use <- n_samples_target_orig # Fallback to LOOCV if N < 5
      }
      cat(paste0("Fold not specified or invalid. n_samples (", n_samples_target_orig, ") >= p_vars/2 (", round(p_vars/2,1), 
                 "). Using ", actual_fold_to_use, "-fold CV for group '", group_name_target, "'.\n"))
    }
  }
  # Ensure fold is at least 1 if there are samples, to avoid issues with max(3,0)
  if (n_samples_target_orig > 0 && actual_fold_to_use < 1) {
      actual_fold_to_use <- 1 # Should be n_samples_target_orig if LOOCV was chosen for N=1
      warning(paste0("Adjusted actual_fold_to_use to 1 for group '", group_name_target, 
                     "' due to insufficient samples (N=", n_samples_target_orig, ").\n"))
  }


  min_samples_for_cv <- max(3, actual_fold_to_use) 
  if(nrow(original_data_matrix_target_numeric) < min_samples_for_cv || ncol(original_data_matrix_target_numeric) == 0) {
      warning(paste0("Insufficient data for target group '", group_name_target, 
                     "' (rows: ", nrow(original_data_matrix_target_numeric), ", cols: ", ncol(original_data_matrix_target_numeric), 
                     ", min_required_for_cv (max(3, actual_fold=", actual_fold_to_use,")): ", min_samples_for_cv,
                     "). Returning default tau=0.05.\n"))
      return(0.05)
  }
  colnames(original_data_matrix_target_numeric) <- otu_labels # 确保列名一致
  
  original_data_matrix_other_numeric <- NULL
  if (!is.null(other_group_name) && other_group_name %in% names(data_list_unscaled_full) && 
      !is.null(data_list_unscaled_full[[other_group_name]])) {
    original_data_matrix_other_unscaled <- data_list_unscaled_full[[other_group_name]]
    temp_other_numeric <- apply(original_data_matrix_other_unscaled, 2, as.numeric)
    # For other group, min_samples_for_cv also depends on its own actual_fold_to_use if we were to be fully symmetric,
    # but for simplicity, we use the target group's min_samples_for_cv or a fixed small number.
    # The current check is just nrow(temp_other_numeric) >= min_samples_for_cv (derived from target's fold)
    if(is.matrix(temp_other_numeric) && nrow(temp_other_numeric) >= min_samples_for_cv && ncol(temp_other_numeric) > 0) {
        original_data_matrix_other_numeric <- temp_other_numeric
        colnames(original_data_matrix_other_numeric) <- otu_labels
    } else {
        warning(paste0("Other group '", other_group_name, 
                         "' has insufficient data (rows: ", if(is.matrix(temp_other_numeric)) nrow(temp_other_numeric) else "N/A",
                         ", min_required_for_cv: ", min_samples_for_cv, 
                         "). It will not be included in internal stabENG calls.\n"))
    }
  }

  # p_vars and n_samples_target_orig are already defined above
  
  # 生成 taulist_cv (在重排循环外部执行一次)
  if (!is.null(initial_pcor_matrix_target) && sum(initial_pcor_matrix_target != 0, na.rm=TRUE) > 0) {
    non_zero_abs_pcors <- abs(initial_pcor_matrix_target[upper.tri(initial_pcor_matrix_target) & initial_pcor_matrix_target != 0])
    if (length(non_zero_abs_pcors) > 0) {
      # Example: Explore smaller tau values more densely
      custom_probs <- sort(unique(c(seq(0.001, 0.02, length.out = 5), seq(0.025, 0.1, length.out = 4), seq(0.15, 0.5, length.out = 3))))
      # Ensure probs are within [0,1] and there are some values
      custom_probs <- custom_probs[custom_probs >= 0 & custom_probs <= 1]
      if(length(custom_probs) == 0) custom_probs <- c(0.01, 0.05, 0.1, 0.2, 0.3) # Fallback

      taulist_cv <- stats::quantile(non_zero_abs_pcors, probs = custom_probs, na.rm = TRUE)
      taulist_cv <- sort(unique(round(taulist_cv, 4))) # Keep rounding to 4 decimal places
      # Ensure some minimal list if quantiles result in too few unique values or NAs
      if(length(taulist_cv) < 3 || all(is.na(taulist_cv))) {
         taulist_cv <- c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1) # More robust fallback
      }
    } else { # non_zero_abs_pcors is empty
      taulist_cv <- c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1) 
    }
  } else { # initial_pcor_matrix_target is NULL or all zeros
    warning(paste0("Initial pcor for group '", group_name_target, "' is NULL or all zeros. Using default taulist_cv."))
    taulist_cv <- c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1) 
  }
  # Final check to ensure taulist_cv is not empty and has reasonable values
  if(length(taulist_cv) == 0 || all(taulist_cv < 1e-5 & taulist_cv > -1e-5)) { # Check if all effectively zero
     taulist_cv <- c(0.001, 0.005, 0.01, 0.05)
  }
  taulist_cv <- taulist_cv[taulist_cv > 0] # Ensure positive tau values
  if(length(taulist_cv) == 0) taulist_cv <- c(0.001, 0.01, 0.05) # Ultimate fallback

  all_permutations_mean_cv_losses <- list()
  
  cat(paste0("Starting ", n_data_permutations, " data permutation(s) for PCS CV (Lasso) for group '", group_name_target,"' threshold tau...\n"))
  
  # --- 外层循环：数据重排 ---
  for (p_iter in 1:n_data_permutations) {
    gc(); print(paste("Start of permutation iter", p_iter, ":", pryr::mem_used()))
    
    # 对当前目标组数据进行行重排
    current_data_matrix_target_numeric <- original_data_matrix_target_numeric
    if (nrow(current_data_matrix_target_numeric) > 1) {
        perm_indices_target <- sample(nrow(current_data_matrix_target_numeric))
        current_data_matrix_target_numeric <- current_data_matrix_target_numeric[perm_indices_target, , drop = FALSE]
    }
    
    # 对“另一个组”的数据进行行重排（如果存在且使用）
    current_data_matrix_other_numeric <- original_data_matrix_other_numeric
    if (!is.null(current_data_matrix_other_numeric) && nrow(current_data_matrix_other_numeric) > 1) {
        perm_indices_other <- sample(nrow(current_data_matrix_other_numeric))
        current_data_matrix_other_numeric <- current_data_matrix_other_numeric[perm_indices_other, , drop = FALSE]
    }
    n_samples_target_current_perm <- nrow(current_data_matrix_target_numeric)


    # --- 内部CV逻辑 (与之前类似，但使用 current_data_matrix_target_numeric) ---
    # Ensure actual_fold_to_use is not zero if n_samples_target_current_perm > 0
    current_iter_fold <- actual_fold_to_use
    if (n_samples_target_current_perm > 0 && current_iter_fold == 0) {
        current_iter_fold <- 1 # Should not happen if logic above is correct
    }
    if (n_samples_target_current_perm == 0 && current_iter_fold > 0) {
        current_iter_fold <- 0 # No samples, no folds for cv.part_internal
    }


    cv_losses_all_folds_current_perm <- matrix(Inf, nrow = length(taulist_cv), ncol = if(current_iter_fold==0) 1 else current_iter_fold) # Handle ncol=0 case for matrix
    rownames(cv_losses_all_folds_current_perm) <- as.character(round(taulist_cv,4))
    cat(paste0("  Permutation ", p_iter, ": Candidate tau values: ", paste(round(taulist_cv,4), collapse=", "), "\n"))
    pb_cv_inner <- utils::txtProgressBar(min = 0, max = length(taulist_cv) * (if(current_iter_fold==0) 1 else current_iter_fold), style = 3) # Adjust max for progress bar
    progress_count_inner = 0

    part.list_cv_indices_target_current_perm <- cv.part_internal(n_samples_target_current_perm, current_iter_fold)
    
    part.list_cv_indices_other_current_perm <- NULL
    if (!is.null(current_data_matrix_other_numeric)) {
        n_samples_other_current_perm <- nrow(current_data_matrix_other_numeric)
        # min_samples_for_cv for other group should ideally be based on its own N and P,
        # but for simplicity, we use the target's min_samples_for_cv or a fixed threshold.
        # The current min_samples_for_cv is derived from target group's fold.
        if (n_samples_other_current_perm >= min_samples_for_cv) { # min_samples_for_cv is max(3, actual_fold_to_use for target)
            part.list_cv_indices_other_current_perm <- cv.part_internal(n_samples_other_current_perm, current_iter_fold) # Use same fold for simplicity
        } else { # not gonna happen
          warning(paste0("Other group '", other_group_name, 
                          "' has insufficient data (rows: ", n_samples_other_current_perm, 
                          ", min_required_for_cv: ", min_samples_for_cv, 
                          "). It will not be included in internal stabENG calls for this permutation.\n"))
            # warning for this specific permutation if other group becomes too small
            # data_matrix_other_numeric will be NULL for stabENG call below for this perm
        }
    }
    
    # --- k-折交叉验证循环 ---
    if (current_iter_fold > 0) {
      for (k_fold in 1:current_iter_fold) {
        cat(paste0("Permutation ", p_iter, "/", n_data_permutations, " for group '", group_name_target, "'...\n"))
        gc(); print(paste("Permutation", p_iter, "Fold", k_fold, "start:", pryr::mem_used()))

        train_indices_target <- part.list_cv_indices_target_current_perm$trainIndices[[k_fold]]
        test_indices_target <- part.list_cv_indices_target_current_perm$testIndices[[k_fold]]

        min_train_samples = max(3, if(!is.null(stabENG_params_list$min.obs)) stabENG_params_list$min.obs else 3) # min.obs for stabENG
        if(length(train_indices_target) < min_train_samples || length(test_indices_target) < 1){
            # warning(paste("Perm", p_iter, "Fold", k_fold, "target", group_name_target, "skip: insufficient train/test."))
            cv_losses_all_folds_current_perm[, k_fold] <- Inf 
            progress_count_inner <- progress_count_inner + length(taulist_cv) 
            utils::setTxtProgressBar(pb_cv_inner, progress_count_inner); next
        }

        x_train_unscaled_fold_target <- current_data_matrix_target_numeric[train_indices_target, , drop = FALSE]
      
        if (sum(apply(x_train_unscaled_fold_target, 2, stats::var, na.rm = TRUE) > 1e-6, na.rm = TRUE) < 2) {
          warning(paste("Perm", p_iter, "Fold", k_fold, "target", group_name_target, "skip: train data < 2 cols with variance."))
          cv_losses_all_folds_current_perm[, k_fold] <- Inf; progress_count_inner <- progress_count_inner + length(taulist_cv); utils::setTxtProgressBar(pb_cv_inner, progress_count_inner); next
        }
      
        train_means_target <- colMeans(x_train_unscaled_fold_target, na.rm = TRUE)
        train_sds_target <- apply(x_train_unscaled_fold_target, 2, stats::sd, na.rm = TRUE)
        train_sds_target[is.na(train_sds_target) | train_sds_target < 1e-6] <- 1
      
        x_train_scaled_fold_target <- sweep(x_train_unscaled_fold_target, 2, train_means_target, "-")
        x_train_scaled_fold_target <- sweep(x_train_scaled_fold_target, 2, train_sds_target, "/")
        x_train_scaled_fold_target[is.na(x_train_scaled_fold_target) | !is.finite(x_train_scaled_fold_target)] <- 0
      
        x_test_unscaled_fold_target <- current_data_matrix_target_numeric[test_indices_target, , drop = FALSE]
        x_test_scaled_fold_target <- sweep(x_test_unscaled_fold_target, 2, train_means_target, "-")
        x_test_scaled_fold_target <- sweep(x_test_scaled_fold_target, 2, train_sds_target, "/")
        x_test_scaled_fold_target[is.na(x_test_scaled_fold_target) | !is.finite(x_test_scaled_fold_target)] <- 0

        stabENG_data_list_train <- list()
        stabENG_data_list_train[[group_name_target]] <- as.data.frame(x_train_unscaled_fold_target)

        if (!is.null(current_data_matrix_other_numeric) && !is.null(part.list_cv_indices_other_current_perm)) {
          train_indices_other <- part.list_cv_indices_other_current_perm$trainIndices[[k_fold]]
          if (length(train_indices_other) >= min_train_samples) {
              x_train_unscaled_fold_other <- current_data_matrix_other_numeric[train_indices_other, , drop = FALSE]
              if (sum(apply(x_train_unscaled_fold_other, 2, stats::var, na.rm = TRUE) > 1e-6, na.rm = TRUE) >= 2) {
                   stabENG_data_list_train[[other_group_name]] <- as.data.frame(x_train_unscaled_fold_other)
              }
          }
        }
      
        current_stabENG_params <- stabENG_params_list
        current_stabENG_params$Y <- stabENG_data_list_train
        current_stabENG_params$labels <- otu_labels
      
        # 考虑不减少或较少地减少这些参数，以增强稳定性
        # 例如，使用原始值的75%或保持不变，而不是减半，但确保不低于一个合理的最小值（如10或15）
        # 示例：使用原始值的75%，但至少为15 (rep.num) 或 10 (nlambda)
        original_rep_num <- stabENG_params_list$rep.num
        original_nlambda1 <- stabENG_params_list$nlambda1
        original_nlambda2 <- stabENG_params_list$nlambda2

        # current_stabENG_params$rep.num <- max(15, round(original_rep_num * 0.75)) 
        # current_stabENG_params$nlambda1 <- max(10, round(original_nlambda1 * 0.75)) 
        # current_stabENG_params$nlambda2 <- max(10, round(original_nlambda2 * 0.75))
      
        # # 或者，如果计算资源允许，可以考虑完全不减少：
        current_stabENG_params$rep.num <- original_rep_num
        current_stabENG_params$nlambda1 <- original_nlambda1
        current_stabENG_params$nlambda2 <- original_nlambda2

        network_results_train <- tryCatch(do.call(stabENG, current_stabENG_params),
                                        error = function(e) { NULL })

        if(is.null(network_results_train) || is.null(network_results_train$opt.fit.pcor) || 
           is.null(network_results_train$opt.fit.pcor[[group_name_target]])){
          cv_losses_all_folds_current_perm[, k_fold] <- Inf; progress_count_inner <- progress_count_inner + length(taulist_cv); utils::setTxtProgressBar(pb_cv_inner, progress_count_inner); next
        }
        pcor_matrix_train_fold_target <- network_results_train$opt.fit.pcor[[group_name_target]]
        if(!is.matrix(pcor_matrix_train_fold_target) || !all(dim(pcor_matrix_train_fold_target) == c(p_vars, p_vars))){
          cv_losses_all_folds_current_perm[, k_fold] <- Inf; progress_count_inner <- progress_count_inner + length(taulist_cv); utils::setTxtProgressBar(pb_cv_inner, progress_count_inner); next
        }
        colnames(pcor_matrix_train_fold_target) <- rownames(pcor_matrix_train_fold_target) <- otu_labels

        # --- tau 循环 ---
        for (j_tau in 1:length(taulist_cv)) {
          tau_current <- taulist_cv[j_tau]
          
           beta_train_for_tau_list <- tryCatch(
             Psi_Screen_Beta_Lasso(x_train_scaled_fold_target, pcor_matrix_train_fold_target, tau_current),
             error = function(e) { NULL })

          if(!is.null(beta_train_for_tau_list) && !is.null(beta_train_for_tau_list$beta_matrix)){
              beta_train_for_tau = beta_train_for_tau_list$beta_matrix
              if(is.matrix(beta_train_for_tau) && all(dim(beta_train_for_tau) == c(p_vars, p_vars)) && 
                 !any(is.na(beta_train_for_tau)) && !any(is.infinite(beta_train_for_tau))){
                
                prediction_result <- tryCatch(
                    x_test_scaled_fold_target %*% beta_train_for_tau,
                    error = function(e) { matrix(Inf, nrow=nrow(x_test_scaled_fold_target), ncol=ncol(x_test_scaled_fold_target)) }
                )
                if (!is.null(prediction_result) && all(is.finite(prediction_result))) {
                    prediction_errors_sq <- (x_test_scaled_fold_target - prediction_result)^2
                    cv_losses_all_folds_current_perm[j_tau, k_fold] <- mean(prediction_errors_sq, na.rm = TRUE) 
                } else {
                    cv_losses_all_folds_current_perm[j_tau, k_fold] <- Inf 
                }
              } else { 
                cv_losses_all_folds_current_perm[j_tau, k_fold] <- Inf 
              }
          } else { 
            cv_losses_all_folds_current_perm[j_tau, k_fold] <- Inf 
          }
          progress_count_inner = progress_count_inner + 1
          utils::setTxtProgressBar(pb_cv_inner, progress_count_inner)
        } # End j_tau loop
        gc(); print(paste("Permutation", p_iter, "Fold", k_fold, "end:", pryr::mem_used()))
        progress_count_inner <- progress_count_inner + length(taulist_cv) 
        utils::setTxtProgressBar(pb_cv_inner, progress_count_inner)
      } # End k_fold loop
    } else { # current_iter_fold is 0 (e.g. n_samples_target_current_perm is 0)
        # All losses remain Inf for this permutation if no folds are run
        progress_count_inner = length(taulist_cv) * 1 # Assume 1 "fold" of doing nothing
        utils::setTxtProgressBar(pb_cv_inner, progress_count_inner)
    }
    close(pb_cv_inner)
    
    all_permutations_mean_cv_losses[[p_iter]] <- rowMeans(cv_losses_all_folds_current_perm, na.rm = TRUE)
    gc(); print(paste("End of permutation iter", p_iter, ":", pryr::mem_used()))
  } # End p_iter (data permutation) loop
  
  # --- 汇总结果并确定最终 optimal_tau ---
  final_mean_cv_loss_per_tau <- rep(Inf, length(taulist_cv))
  if (length(all_permutations_mean_cv_losses) == n_data_permutations && 
      all(sapply(all_permutations_mean_cv_losses, length) == length(taulist_cv))) {
      
      losses_matrix <- do.call(cbind, all_permutations_mean_cv_losses)
      final_mean_cv_loss_per_tau <- rowMeans(losses_matrix, na.rm = TRUE)
  } else {
      warning(paste0("Could not properly aggregate CV losses from all ",n_data_permutations," permutations for group '", group_name_target, "'. Fallback logic will be used."))
  }
  
  optimal_tau <- 0.05 # Initial default
  all_final_losses_invalid <- all(is.na(final_mean_cv_loss_per_tau) | is.infinite(final_mean_cv_loss_per_tau))
  
  if (all_final_losses_invalid) {
      warning(paste0("PCS CV (Lasso) for target group '", group_name_target, "' (across ",n_data_permutations," permutations): All aggregated mean CV losses are NA or Inf."))
      # Fallback logic using initial_pcor_matrix_target (original, unpermuted)
      if (!is.null(initial_pcor_matrix_target) && sum(initial_pcor_matrix_target != 0, na.rm=TRUE) > 0) {
          non_zero_pcors <- abs(initial_pcor_matrix_target[upper.tri(initial_pcor_matrix_target) & initial_pcor_matrix_target!=0])
          if(length(non_zero_pcors) > 0) {
            optimal_tau <- stats::quantile(non_zero_pcors, 0.05, na.rm = TRUE) 
            optimal_tau <- max(0.001, optimal_tau, na.rm=TRUE) 
            optimal_tau <- min(0.5, optimal_tau, na.rm=TRUE) 
            if(is.na(optimal_tau) || is.infinite(optimal_tau)) optimal_tau <- 0.05 
            warning(paste0("  Using quantile-based default tau = ", sprintf("%.4f", optimal_tau), " for group '", group_name_target, "' due to problematic aggregated CV losses.\n"))
          } else {
            warning(paste0("  Initial pcor for group '", group_name_target, "' has no non-zero off-diagonal. Using fixed default tau = 0.05 for aggregated CV.\n"))
            optimal_tau <- 0.05
          }
      } else {
        warning(paste0("  Initial pcor for group '", group_name_target, "' is NULL or all zeros. Using fixed default tau = 0.05 for aggregated CV.\n"))
        optimal_tau <- 0.05
      }
  } else {
      finite_losses_indices <- which(is.finite(final_mean_cv_loss_per_tau))
      if(length(finite_losses_indices) > 0){
          best_idx_in_finite <- which.min(final_mean_cv_loss_per_tau[finite_losses_indices])
          optimal_tau <- taulist_cv[finite_losses_indices[best_idx_in_finite]]
      } else {
          warning(paste0("PCS CV (Lasso) for target group '", group_name_target, "' (aggregated): No finite mean CV losses. Using fixed default tau = 0.05.\n"))
          optimal_tau <- 0.05
      }
  }

  # --- 绘图：包含所有重排的曲线 ---
  if (plot_cv_curve && !is.null(plot_path_prefix) && length(all_permutations_mean_cv_losses) > 0 && n_data_permutations > 0){
    png_filename <- paste0(plot_path_prefix, "_pcs_cv_curve_group_", group_name_target, "_with_permutations.png")
    tryCatch({
        grDevices::png(png_filename, width=1000, height=700) # Wider plot for legend
        
        # Determine Y-axis limits based on all finite loss values from all permutations
        all_finite_losses_values <- unlist(lapply(all_permutations_mean_cv_losses, function(perm_loss_vector) {
            if (is.numeric(perm_loss_vector)) perm_loss_vector[is.finite(perm_loss_vector)] else numeric(0)
        }))
        
        plot_ylim <- if(length(all_finite_losses_values) > 0) {
            range(all_finite_losses_values, na.rm=TRUE)
        } else {
            warning(paste0("No finite CV losses found for group '", group_name_target, "'. Using default Y-axis limits."))
            c(0,1) # Fallback if no finite values at all
        }
        if(any(!is.finite(plot_ylim))) plot_ylim <- c(0,1) # Further fallback if range results in Inf/NA

        # 使用预定义的颜色调色板
        # Ensure enough unique colors, fallback to rainbow if Okabe-Ito has fewer than needed
        perm_colors <- grDevices::palette.colors(n = max(1, n_data_permutations), palette = "Okabe-Ito")
        if(length(perm_colors) < n_data_permutations && n_data_permutations > 0) {
            perm_colors <- grDevices::rainbow(n_data_permutations)
        }
        if(n_data_permutations == 0) perm_colors <- "black" # Should not happen if length check passed

        # Plot the first permutation's curve
        if (length(all_permutations_mean_cv_losses) >= 1 && 
            !is.null(all_permutations_mean_cv_losses[[1]]) && 
            length(all_permutations_mean_cv_losses[[1]]) == length(taulist_cv)) {
            
            graphics::plot(taulist_cv, all_permutations_mean_cv_losses[[1]], type="l", 
                 xlab="Screening Threshold (tau)", ylab="Mean CV Prediction Error (MSE)", 
                 main=paste("PCS CV Curves (",n_data_permutations," Permutations) - Group:", group_name_target),
                 ylim=plot_ylim, col=perm_colors[1], lwd=1.5, xaxt = "n")
            graphics::axis(1, at = taulist_cv, labels = round(taulist_cv,3), las=2, cex.axis=0.8)
        } else {
             # Fallback plot if first permutation data is problematic
            graphics::plot(0, 0, type="n", 
                 xlab="Screening Threshold (tau)", ylab="Mean CV Prediction Error (MSE)", 
                 main=paste("PCS CV Curves (Error) - Group:", group_name_target),
                 xlim=range(taulist_cv, na.rm=TRUE), ylim=plot_ylim, xaxt = "n")
            graphics::axis(1, at = taulist_cv, labels = round(taulist_cv,3), las=2, cex.axis=0.8)
            graphics::text(mean(graphics::par("usr")[1:2]), mean(graphics::par("usr")[3:4]), "Error: Data for permutation 1 is invalid.", col="red")
        }

        # Plot curves for subsequent permutations
        if (n_data_permutations > 1) {
            for (p_idx in 2:n_data_permutations) {
                if(length(all_permutations_mean_cv_losses) >= p_idx && 
                   !is.null(all_permutations_mean_cv_losses[[p_idx]]) &&
                   length(all_permutations_mean_cv_losses[[p_idx]]) == length(taulist_cv)){
                     graphics::lines(taulist_cv, all_permutations_mean_cv_losses[[p_idx]], type="l", 
                                col=perm_colors[p_idx], lwd=1.5) # Direct indexing for color
                }
            }
        }
        
        # 绘制平均曲线
        if(any(is.finite(final_mean_cv_loss_per_tau))){
            graphics::lines(taulist_cv, final_mean_cv_loss_per_tau, type="l", col="black", lwd=2.5, lty=2)
        }

        # Prepare legend items
        legend_texts <- character(0)
        legend_cols <- character(0)
        legend_ltys <- numeric(0)
        legend_lwds <- numeric(0)

        if (n_data_permutations > 0) {
            legend_texts <- c(legend_texts, paste("Perm.", 1:n_data_permutations))
            legend_cols <- c(legend_cols, perm_colors[1:n_data_permutations])
            legend_ltys <- c(legend_ltys, rep(1, n_data_permutations))
            legend_lwds <- c(legend_lwds, rep(1.5, n_data_permutations))
        }
        
        if(any(is.finite(final_mean_cv_loss_per_tau))){
            legend_texts <- c(legend_texts, "Average Curve")
            legend_cols <- c(legend_cols, "black")
            legend_ltys <- c(legend_ltys, 2)
            legend_lwds <- c(legend_lwds, 2.5)
        }

        if(!all_final_losses_invalid && length(which(is.finite(final_mean_cv_loss_per_tau))) > 0){
            graphics::abline(v=optimal_tau, col="red", lty=2, lwd=2)
            legend_texts <- c(legend_texts, paste("Optimal tau =", round(optimal_tau,4)))
            legend_cols <- c(legend_cols, "red")
            legend_ltys <- c(legend_ltys, 2)
            legend_lwds <- c(legend_lwds, 2)
        } else if (n_data_permutations > 0) { # Only add this text if permutations were attempted
            graphics::text(mean(graphics::par("usr")[1:2]), mean(graphics::par("usr")[3:4]), "CV failed to find optimal tau.\nUsing default.", col="red", cex=1.2)
        }
        
        if(length(legend_texts) > 0){
            graphics::legend("topright", legend=legend_texts, col=legend_cols, lty=legend_ltys, lwd=legend_lwds, bg="white", cex=0.7)
        }
        
        grDevices::dev.off()
        cat(paste("\n  PCS CV curve plot (with permutations) saved to:", png_filename, "\n"))
    }, error = function(e){
        warning(paste0("Error plotting PCS CV curve (with permutations) for group '", group_name_target, "': ", conditionMessage(e)))
        if(grDevices::dev.cur() != 1) grDevices::dev.off()
    })
  }
  gc(); print(paste("End of pcs_cv_threshold_stabENG_lasso_perm:", pryr::mem_used()))
  cat(sprintf("\nFinal Optimal PCS (Lasso) screening threshold tau for group '%s' (from %d permutations) = %.4f\n", 
              group_name_target, n_data_permutations, optimal_tau))
  return(optimal_tau)
}

#—— 基于PCS论文的Permutation Screen——
# pcor_screen_pcs = function(data_list, labels,
#                            group_name,
#                            n_perm = 50,
#                            alpha = 0.1,
#                            lambda1, lambda2) {
#   # data_list: 原始 count 矩阵列表，data_list[[group_name]] 为目标组
#   # labels: 列名（OTU）
#   # lambda1/lambda2: stabENG 中选好的参数
#   null_pcors <- numeric(0)
#   for(p in 1:n_perm) {
#     # 1) 置换每列
#     perm_list <- lapply(data_list, function(df){
#       df2 <- as.matrix(df)
#       for(j in seq_len(ncol(df2))) df2[,j] <- sample(df2[,j])
#       as.data.frame(df2)
#     })
#     # 2) 用原始 stabENG+glasso 估计 pcor
#     res_p <- preprocess_and_estimate_network(perm_list,
#                  labels = labels,
#                  lambda1 = lambda1,
#                  lambda2 = lambda2)$pcor[[group_name]]
#     # 3) 收集 upper‐triangular 的绝对值
#     null_pcors <- c(null_pcors, abs(res_p[upper.tri(res_p)]))
#   }
#   # 4) 阈值为 1−α 分位数
#   quantile(null_pcors, probs = 1 - alpha, na.rm = TRUE)
# }

# # —— 在 pcs_cv_threshold_stabENG_lasso_perm 开头或最终 fallback 里调用 ——  
# # 假设你已经跑完 stabENG，拿到了 opt.lambda1 和 opt.lambda2
# tau_pcs_paper <- pcor_screen_pcs(
#     data_list         = data_list_unscaled_full,
#     labels            = otu_labels,
#     group_name        = group_name_target,
#     n_perm            = 50,       # 置换次数
#     alpha             = 0.1,      # false‐positive rate
#     lambda1           = stabENG_params_list$opt.lambda1,
#     lambda2           = stabENG_params_list$opt.lambda2
# )
# cat("PCS paper permutation get τ =", round(tau_pcs_paper,5), "\n")
# return(tau_pcs_paper)