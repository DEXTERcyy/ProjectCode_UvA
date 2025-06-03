# 辅助函数：K折划分
cv.part_internal = function(n, k) {
  k_orig_for_warning <- k # 保存原始k用于警告信息

  if (n == 0) {
    return(list(trainIndices = list(), testIndices = list()))
  }

  # 调整k值以确保其有效性
  if (n == 1) {
    if (k_orig_for_warning != 1) {
      # warning(paste0("For n=1, k must be 1. Adjusted k from ", k_orig_for_warning, " to 1."))
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
          # warning(paste("Response variable", response_var_name, "has zero or near-zero variance. Skipping Lasso."))
          next
      }
      # Ensure x_predictors have variance (glmnet handles this, but good to be aware)
      # col_vars_predictors <- apply(x_predictors, 2, stats::var, na.rm = TRUE)
      # if(any(col_vars_predictors < .Machine$double.eps^0.5)){
      #     warning(paste("Some predictors for response", response_var_name, "have zero variance. They will be handled by glmnet."))
      # }


      if(ncol(x_predictors) > 0){
        # nfolds must be >=3 and <= number of observations.
        current_nfolds <- if (n_samples_current_x < 3) { 0 } else { min(n_samples_current_x, 5) }
        if (current_nfolds < 3 && n_samples_current_x >=3) { current_nfolds <- 3 } # Ensure at least 3 if N allows

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
        } else {
          # warning(paste("Not enough samples (",n_samples_current_x,") for cv.glmnet (nfolds=",current_nfolds,") for response", response_var_name, ". Skipping Lasso."))
        }
      }
    }
  }
  return(list(beta_matrix = beta_matrix_final)) # Ensure named list for consistency
}

# --- PCS 交叉验证阈值选择函数 (使用 Lasso 进行内部回归) ---
pcs_cv_threshold_stabENG_lasso <- function(
  data_list_unscaled_full, 
  initial_pcor_matrix_target, 
  group_name_target, 
  other_group_name, 
  stabENG_params_list,
  otu_labels,
  fold = 5, 
  plot_cv_curve = FALSE,
  plot_path_prefix = NULL
) {
  
  if (is.null(group_name_target) || !group_name_target %in% names(data_list_unscaled_full) || 
      is.null(data_list_unscaled_full[[group_name_target]])) {
    warning(paste0("Target group '", group_name_target, "' not found or is NULL. Returning default tau=0.05."))
    return(0.05)
  }
  data_matrix_target_unscaled <- data_list_unscaled_full[[group_name_target]]
  data_matrix_target_numeric <- apply(data_matrix_target_unscaled, 2, as.numeric)
  if(!is.matrix(data_matrix_target_numeric)) data_matrix_target_numeric <- as.matrix(data_matrix_target_numeric)
  
  min_samples_for_cv <- max(3, fold) # Minimum samples needed for meaningful CV
  if(nrow(data_matrix_target_numeric) < min_samples_for_cv || ncol(data_matrix_target_numeric) == 0) {
      warning(paste0("Insufficient data for target group '", group_name_target, 
                     "' (rows: ", nrow(data_matrix_target_numeric), ", cols: ", ncol(data_matrix_target_numeric), 
                     ", min_required_for_cv: ", min_samples_for_cv,
                     "). Returning default tau=0.05.\n"))
      return(0.05)
  }
  colnames(data_matrix_target_numeric) <- otu_labels
  
  data_matrix_other_numeric <- NULL
  if (!is.null(other_group_name) && other_group_name %in% names(data_list_unscaled_full) && 
      !is.null(data_list_unscaled_full[[other_group_name]])) {
    data_matrix_other_unscaled <- data_list_unscaled_full[[other_group_name]]
    temp_other_numeric <- apply(data_matrix_other_unscaled, 2, as.numeric)
    if(is.matrix(temp_other_numeric) && nrow(temp_other_numeric) >= min_samples_for_cv && ncol(temp_other_numeric) > 0) {
        data_matrix_other_numeric <- temp_other_numeric
        colnames(data_matrix_other_numeric) <- otu_labels
    } else {
        warning(paste0("Other group '", other_group_name, 
                         "' has insufficient data for CV. It will not be included in internal stabENG calls.\n"))
    }
  }

  p_vars <- ncol(data_matrix_target_numeric)
  n_samples_target <- nrow(data_matrix_target_numeric)
  
  if (!is.null(initial_pcor_matrix_target) && sum(initial_pcor_matrix_target != 0, na.rm=TRUE) > 0) {
    non_zero_abs_pcors <- abs(initial_pcor_matrix_target[upper.tri(initial_pcor_matrix_target) & initial_pcor_matrix_target != 0])
    if (length(non_zero_abs_pcors) > 0) {
      taulist_cv <- stats::quantile(non_zero_abs_pcors, probs = seq(0.01, 0.5, length.out = 10), na.rm = TRUE)
      taulist_cv <- sort(unique(round(taulist_cv, 4)))
      if(length(taulist_cv) == 0) taulist_cv <- c(0.005, 0.01, 0.05, 0.1) 
    } else {
      taulist_cv <- c(0.005, 0.01, 0.05, 0.1) 
    }
  } else {
    warning(paste0("Initial pcor for group '", group_name_target, "' is NULL or all zeros. Using default taulist_cv."))
    taulist_cv <- c(0.005, 0.01, 0.05, 0.1) 
  }
  if(length(taulist_cv) == 0 || all(taulist_cv < 1e-5)) taulist_cv <- c(0.005, 0.01, 0.05) # Ensure some reasonable taus

  cv_losses_all_folds <- matrix(Inf, nrow = length(taulist_cv), ncol = fold)
  rownames(cv_losses_all_folds) <- as.character(round(taulist_cv,4))
  
  cat(paste0("  Starting PCS CV (Lasso) for group '", group_name_target,"' threshold tau...\n"))
  cat(paste0("  Candidate tau values: ", paste(round(taulist_cv,4), collapse=", "), "\n"))
  pb_cv <- utils::txtProgressBar(min = 0, max = length(taulist_cv) * fold, style = 3)
  progress_count = 0

  part.list_cv_indices_target <- cv.part_internal(n_samples_target, fold)
  
  part.list_cv_indices_other <- NULL
  if (!is.null(data_matrix_other_numeric)) {
      n_samples_other <- nrow(data_matrix_other_numeric)
      if (n_samples_other >= min_samples_for_cv) {
          part.list_cv_indices_other <- cv.part_internal(n_samples_other, fold)
      } else {
          warning(paste0("Other group '", other_group_name, 
                           "' has too few samples (", n_samples_other, ") for ", fold, 
                           "-fold CV. It will not be used in stabENG calls.\n"))
          data_matrix_other_numeric <- NULL 
      }
  }

  for (k_fold in 1:fold) {
    train_indices_target <- part.list_cv_indices_target$trainIndices[[k_fold]]
    test_indices_target <- part.list_cv_indices_target$testIndices[[k_fold]]

    min_train_samples = max(3, stabENG_params_list$min.obs %||% 3) # min.obs might be a stabENG param
    if(length(train_indices_target) < min_train_samples || length(test_indices_target) < 1){
        warning(paste("Fold", k_fold, "target group", group_name_target, "skip: insufficient train/test (train:", length(train_indices_target), ", test:", length(test_indices_target), ")."))
        cv_losses_all_folds[, k_fold] <- Inf # Already initialized to Inf
        progress_count = progress_count + length(taulist_cv); utils::setTxtProgressBar(pb_cv, progress_count); next
    }

    x_train_unscaled_fold_target <- data_matrix_target_numeric[train_indices_target, , drop = FALSE]
    
    if (sum(apply(x_train_unscaled_fold_target, 2, stats::var, na.rm = TRUE) > 1e-6, na.rm = TRUE) < 2) { # Check for at least 2 varying columns
        warning(paste("Fold", k_fold, "target group", group_name_target, "skip: training data has < 2 columns with variance."))
        cv_losses_all_folds[, k_fold] <- Inf; progress_count = progress_count + length(taulist_cv); utils::setTxtProgressBar(pb_cv, progress_count); next
    }
    
    # Standardize target group data
    # Training data scaling
    train_means_target <- colMeans(x_train_unscaled_fold_target, na.rm = TRUE)
    train_sds_target <- apply(x_train_unscaled_fold_target, 2, stats::sd, na.rm = TRUE)
    train_sds_target[is.na(train_sds_target) | train_sds_target < 1e-6] <- 1 # Avoid division by zero/small SD
    
    x_train_scaled_fold_target <- sweep(x_train_unscaled_fold_target, 2, train_means_target, "-")
    x_train_scaled_fold_target <- sweep(x_train_scaled_fold_target, 2, train_sds_target, "/")
    x_train_scaled_fold_target[is.na(x_train_scaled_fold_target) | !is.finite(x_train_scaled_fold_target)] <- 0 # Impute remaining NAs/Infs
    
    # Test data scaling using training parameters
    x_test_unscaled_fold_target <- data_matrix_target_numeric[test_indices_target, , drop = FALSE]
    x_test_scaled_fold_target <- sweep(x_test_unscaled_fold_target, 2, train_means_target, "-")
    x_test_scaled_fold_target <- sweep(x_test_scaled_fold_target, 2, train_sds_target, "/")
    x_test_scaled_fold_target[is.na(x_test_scaled_fold_target) | !is.finite(x_test_scaled_fold_target)] <- 0

    stabENG_data_list_train <- list()
    stabENG_data_list_train[[group_name_target]] <- as.data.frame(x_train_unscaled_fold_target) # stabENG might prefer unscaled

    if (!is.null(data_matrix_other_numeric) && !is.null(part.list_cv_indices_other)) {
        train_indices_other <- part.list_cv_indices_other$trainIndices[[k_fold]]
        if (length(train_indices_other) >= min_train_samples) {
            x_train_unscaled_fold_other <- data_matrix_other_numeric[train_indices_other, , drop = FALSE]
            if (sum(apply(x_train_unscaled_fold_other, 2, stats::var, na.rm = TRUE) > 1e-6, na.rm = TRUE) >= 2) {
                 stabENG_data_list_train[[other_group_name]] <- as.data.frame(x_train_unscaled_fold_other)
            } else {
                 warning(paste("Fold", k_fold, ": Other group '", other_group_name, "' training data has < 2 cols with variance. stabENG runs with K=1."))
            }
        } else {
            warning(paste("Fold", k_fold, ": Other group '", other_group_name, "' has <",min_train_samples,"samples in train fold. stabENG runs with K=1."))
        }
    }
    
    current_stabENG_params <- stabENG_params_list
    current_stabENG_params$Y <- stabENG_data_list_train
    current_stabENG_params$labels <- otu_labels
    # Potentially reduce complexity for inner CV loop stabENG calls
    current_stabENG_params$rep.num <- max(10, round(stabENG_params_list$rep.num / 2))
    current_stabENG_params$nlambda1 <- max(10, round(stabENG_params_list$nlambda1 / 2))
    current_stabENG_params$nlambda2 <- max(10, round(stabENG_params_list$nlambda2 / 2))
    
    network_results_train <- tryCatch(do.call(stabENG, current_stabENG_params),
                                      error = function(e) {
                                          warning(paste("Error stabENG target", group_name_target, "fold", k_fold, ":", conditionMessage(e)))
                                          NULL
                                      })

    if(is.null(network_results_train) || is.null(network_results_train$opt.fit.pcor) || 
       is.null(network_results_train$opt.fit.pcor[[group_name_target]])){
        warning(paste("stabENG NULL/invalid pcor target", group_name_target, "fold", k_fold))
        cv_losses_all_folds[, k_fold] <- Inf; progress_count = progress_count + length(taulist_cv); utils::setTxtProgressBar(pb_cv, progress_count); next
    }
    pcor_matrix_train_fold_target <- network_results_train$opt.fit.pcor[[group_name_target]]
    if(!is.matrix(pcor_matrix_train_fold_target) || !all(dim(pcor_matrix_train_fold_target) == c(p_vars, p_vars))){
        warning(paste("pcor_matrix from stabENG invalid dim for target", group_name_target, "fold", k_fold))
        cv_losses_all_folds[, k_fold] <- Inf; progress_count = progress_count + length(taulist_cv); utils::setTxtProgressBar(pb_cv, progress_count); next
    }
    colnames(pcor_matrix_train_fold_target) <- rownames(pcor_matrix_train_fold_target) <- otu_labels

    for (j_tau in 1:length(taulist_cv)) {
      tau_current <- taulist_cv[j_tau]
      
       beta_train_for_tau_list <- tryCatch(
         Psi_Screen_Beta_Lasso(x_train_scaled_fold_target, pcor_matrix_train_fold_target, tau_current),
         error = function(e) {
             warning(paste("Error Psi_Screen_Beta_Lasso target", group_name_target, "fold", k_fold, "tau", tau_current, ":", conditionMessage(e)))
             NULL
         })

      if(!is.null(beta_train_for_tau_list) && !is.null(beta_train_for_tau_list$beta_matrix)){
          beta_train_for_tau = beta_train_for_tau_list$beta_matrix # Access the named element
          if(is.matrix(beta_train_for_tau) && all(dim(beta_train_for_tau) == c(p_vars, p_vars)) && 
             !any(is.na(beta_train_for_tau)) && !any(is.infinite(beta_train_for_tau))){
              
              prediction_result <- tryCatch(
                  x_test_scaled_fold_target %*% beta_train_for_tau,
                  error = function(e) {
                      warning(paste("Error matrix mult CV loss (target", group_name_target, "fold", k_fold, "tau", tau_current, "):", conditionMessage(e)))
                      matrix(Inf, nrow=nrow(x_test_scaled_fold_target), ncol=ncol(x_test_scaled_fold_target)) # Return Inf matrix on error
                  }
              )
              if (!is.null(prediction_result) && all(is.finite(prediction_result))) {
                  prediction_errors_sq <- (x_test_scaled_fold_target - prediction_result)^2
                  # Use Mean Squared Error
                  cv_losses_all_folds[j_tau, k_fold] <- mean(prediction_errors_sq, na.rm = TRUE) 
              } else {
                  cv_losses_all_folds[j_tau, k_fold] <- Inf 
              }
          } else { 
              warning(paste("beta_train_for_tau invalid/NA/Inf for target", group_name_target, "fold", k_fold, "tau", tau_current))
              cv_losses_all_folds[j_tau, k_fold] <- Inf 
          }
      } else { 
          warning(paste("Psi_Screen_Beta_Lasso returned NULL/invalid for target", group_name_target, "fold", k_fold, "tau", tau_current))
          cv_losses_all_folds[j_tau, k_fold] <- Inf 
      }
      progress_count = progress_count + 1
      utils::setTxtProgressBar(pb_cv, progress_count)
    }
    gc() 
  } # End k_fold loop
  close(pb_cv)
  
  mean_cv_loss_per_tau <- rowMeans(cv_losses_all_folds, na.rm = TRUE) # na.rm=TRUE is important if some folds consistently fail for a tau
  
  optimal_tau <- 0.05 # Initial default
  all_losses_invalid <- all(is.na(mean_cv_loss_per_tau) | is.infinite(mean_cv_loss_per_tau))
  
  if (all_losses_invalid) {
      warning(paste0("PCS CV (Lasso) for target group '", group_name_target, "': All mean CV losses are NA or Inf."))
      if (!is.null(initial_pcor_matrix_target) && sum(initial_pcor_matrix_target != 0, na.rm=TRUE) > 0) {
          non_zero_pcors <- abs(initial_pcor_matrix_target[upper.tri(initial_pcor_matrix_target) & initial_pcor_matrix_target!=0])
          if(length(non_zero_pcors) > 0) {
            optimal_tau <- stats::quantile(non_zero_pcors, 0.05, na.rm = TRUE) 
            optimal_tau <- max(0.001, optimal_tau, na.rm=TRUE) # Ensure it's not too small/zero/NA
            optimal_tau <- min(0.5, optimal_tau, na.rm=TRUE) # Ensure it's not too large
            if(is.na(optimal_tau) || is.infinite(optimal_tau)) optimal_tau <- 0.05 # Ultimate fallback
            warning(paste0("  Using quantile-based default tau = ", sprintf("%.4f", optimal_tau), " for group '", group_name_target, "' due to problematic CV losses.\n"))
          } else {
            warning(paste0("  Initial pcor for group '", group_name_target, "' has no non-zero off-diagonal. Using fixed default tau = 0.05.\n"))
            optimal_tau <- 0.05
          }
      } else {
        warning(paste0("  Initial pcor for group '", group_name_target, "' is NULL or all zeros. Using fixed default tau = 0.05.\n"))
        optimal_tau <- 0.05
      }
  } else {
      # Find minimum among finite losses
      finite_losses_indices <- which(is.finite(mean_cv_loss_per_tau))
      if(length(finite_losses_indices) > 0){
          best_idx_in_finite <- which.min(mean_cv_loss_per_tau[finite_losses_indices])
          optimal_tau <- taulist_cv[finite_losses_indices[best_idx_in_finite]]
      } else {
          # This case should be caught by all_losses_invalid, but as a safeguard:
          warning(paste0("PCS CV (Lasso) for target group '", group_name_target, "': No finite mean CV losses found, though not all were NA/Inf. Using fixed default tau = 0.05.\n"))
          optimal_tau <- 0.05
      }
  }

  if (plot_cv_curve && !is.null(plot_path_prefix)){
    png_filename <- paste0(plot_path_prefix, "_pcs_cv_curve_group_", group_name_target, ".png")
    tryCatch({
        grDevices::png(png_filename, width=800, height=600)
        finite_mean_losses <- mean_cv_loss_per_tau[is.finite(mean_cv_loss_per_tau)]
        plot_ylim <- if(length(finite_mean_losses) > 0) range(finite_mean_losses, na.rm=TRUE) else c(0,1) # Fallback ylim

        graphics::plot(taulist_cv, mean_cv_loss_per_tau, type="b", 
             xlab="Screening Threshold (tau)", ylab="Mean CV Prediction Error (MSE)", 
             main=paste("PCS CV Curve for tau selection - Group:", group_name_target),
             ylim=plot_ylim, pch=19)
        if(!all_losses_invalid && length(which(is.finite(mean_cv_loss_per_tau))) > 0){ # Only add line if optimal_tau was found from CV
            graphics::abline(v=optimal_tau, col="red", lty=2)
            graphics::legend("topright", legend=paste("Optimal tau =", round(optimal_tau,4)), col="red", lty=2, bg="white")
        } else {
            graphics::text(mean(graphics::par("usr")[1:2]), mean(graphics::par("usr")[3:4]), "CV failed to find optimal tau.\nUsing default.", col="red", cex=1.2)
        }
        grDevices::dev.off()
        cat(paste("\n  PCS CV curve plot saved to:", png_filename, "\n"))
    }, error = function(e){
        warning(paste0("Error plotting PCS CV curve for group '", group_name_target, "': ", conditionMessage(e)))
        if(grDevices::dev.cur() != 1) grDevices::dev.off() # Close device if error occurred mid-plot
    })
  }
  
  cat(sprintf("\nOptimal PCS (Lasso) screening threshold tau for group '%s' = %.4f\n", group_name_target, optimal_tau))
  return(optimal_tau)
}