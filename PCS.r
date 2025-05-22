# [[在这里定义辅助函数]]

# 辅助函数：K折划分
cv.part_internal = function(n, k) {
  if (k <= 1 || k > n) { # 确保k是合理的
    warning("Number of folds k is not appropriate. Setting k to min(n, 10) if n > 1, else k=n.")
    k = if (n > 1) min(n, 10) else n
    if (k <= 1 && n > 1) k = 2 # 至少2折如果样本大于1
    else if (k < 1) k = n # 如果样本为1，则k为1
  }
  ntest = floor(n/k)
  if (ntest == 0 && k < n) { # 如果每折测试样本为0，但k<n，则调整
      ntest = 1 
      k = n # 实质上变成 LOOCV 或接近 LOOCV
      warning(paste0("Adjusted k to ", k, " (LOOCV or similar) as n/k was < 1."))
  } else if (ntest == 0 && k == n) { # LOOCV
      ntest = 1
  }


  ind = sample(n)
  train_indices_list <- list()
  test_indices_list <- list()

  for (j_fold in 1:k) {
    if (k == n) { # LOOCV case
        test_sel_indices = ind[j_fold]
    } else {
        start_idx = (j_fold - 1) * ntest + 1
        end_idx = j_fold * ntest
        if (j_fold == k) { # 最后一折可能包含余下的所有样本
             end_idx = n
        }
        test_sel_indices = ind[start_idx:end_idx]
    }
    
    train_sel_indices = setdiff(ind, test_sel_indices)

    train_indices_list[[j_fold]] <- train_sel_indices
    test_indices_list[[j_fold]] <- test_sel_indices
  }
  # 为了与原始Lafit脚本的返回结构兼容，我们仍然创建矩阵
  # 但要注意，如果每折的训练/测试样本数不同（例如最后一折），矩阵填充会有问题
  # 因此，返回列表通常更安全。这里我们先尝试保持矩阵结构，但给出警告。
  # 实际中，直接使用 train_indices_list 和 test_indices_list 会更好。
  
  #  Lafit 脚本的矩阵形式更适合固定大小的折
  #  对于可能不均等大小的折（特别是最后一折），直接用列表更好
  #  这里我们为了兼容性，还是尝试创建矩阵，但要知道这可能不完美
  #  pcs_cv_threshold 中我们直接用列表索引
  
  # 为了与Lafit的 cv_part 兼容，我们还是返回 trainMat 和 testMat
  # 但在 pcs_cv_threshold 中我们应该使用列表索引更安全
  # 这里我们简化，因为 pcs_cv_threshold 已经改为用列表索引了。
  # 所以cv.part_internal 可以直接返回列表。
  # 不过为了和你之前的代码兼容，我们还是返回矩阵形式（但要注意不均等折的问题）
  # 我们假设折大小基本一致，除了最后一折可能稍大
  
  # 修正：直接返回索引列表，因为 pcs_cv_threshold 已适配
  return(list(trainIndices = train_indices_list, testIndices = test_indices_list))
}


# 辅助函数：Psi_Screen_Beta (确保与预测误差计算兼容)
# 输入: x (已标准化的数值矩阵), R (偏相关矩阵，行列名与x的列名对应), tau (筛选阈值)
# 输出: list(beta_matrix)，其中 beta_matrix[j,i] 是 Xj 预测 Xi 时的系数
Psi_Screen_Beta = function(x, R, tau) {
  # 确保x是矩阵
  if (!is.matrix(x)) x <- as.matrix(x)
  
  n_vars <- ncol(x)
  var_names <- colnames(x)
  if (is.null(var_names)) {
    var_names <- paste0("V", 1:n_vars)
    colnames(x) <- var_names
  }
  if (is.null(rownames(R))) rownames(R) <- var_names
  if (is.null(colnames(R))) colnames(R) <- var_names

  beta_matrix_final <- matrix(0, nrow = n_vars, ncol = n_vars,
                              dimnames = list(Predictors = var_names, Response = var_names))

  R_screened <- R
  R_screened[abs(R_screened) <= tau] <- 0
  diag(R_screened) <- 0 # 节点不预测自身

  for (i_response_idx in 1:n_vars) {
    response_var_name <- var_names[i_response_idx]
    
    # 找到节点 i_response_idx 的邻居 (基于筛选后的R_screened)
    # R_screened[predictor_row, response_col] != 0
    neighbor_names <- var_names[which(R_screened[, response_var_name] != 0)] # 列是响应，行是潜在预测

    if (length(neighbor_names) > 0) {
      # 构建数据框和公式
      df_regression <- as.data.frame(x)
      # 确保变量名在公式中被正确引用 (特别是包含特殊字符时)
      response_var_formula <- paste0("`", response_var_name, "`")
      predictor_vars_formula <- paste0("`", neighbor_names, "`", collapse = " + ")
      
      formula_str <- paste(response_var_formula, "~ 0 +", predictor_vars_formula)
      
      model_fit <- tryCatch({
        lm(as.formula(formula_str), data = df_regression, na.action = na.omit)
      }, error = function(e) {
        # warning(paste("Error in lm for response", response_var_name, ":", e$message))
        NULL
      })

      if (!is.null(model_fit) && length(stats::coef(model_fit)) > 0) {
        estimated_coeffs <- stats::coef(model_fit)
        # coef的名称可能因为lm内部处理特殊字符而被修改，要匹配回去
        # estimated_coeffs 的名称是 R 修改后的 predictor_names
        for(pred_name in neighbor_names){
            # 尝试直接匹配，或匹配R的make.names版本
            original_name_in_beta_matrix = pred_name
            coeff_val = NA
            
            if (paste0("`", pred_name, "`") %in% names(estimated_coeffs)) { # R的lm公式会自动加反引号
                coeff_val = estimated_coeffs[paste0("`", pred_name, "`")]
            } else if (pred_name %in% names(estimated_coeffs)) {
                coeff_val = estimated_coeffs[pred_name]
            } else {
                 r_compatible_name = make.names(pred_name)
                 if (r_compatible_name %in% names(estimated_coeffs)){
                     coeff_val = estimated_coeffs[r_compatible_name]
                 } else {
                    # warning(paste("Coefficient for predictor", pred_name, "not found in model for", response_var_name))
                 }
            }
            if(!is.na(coeff_val)){
                 beta_matrix_final[original_name_in_beta_matrix, response_var_name] <- coeff_val
            }
        }
      }
    }
  }
  return(list(beta_matrix_final))
}


# --- 新的 PCS 交叉验证阈值选择函数 (与上一回答中的版本类似，但确保调用了这里定义的辅助函数) ---
pcs_cv_threshold <- function(
  data_matrix,             # 单个条件组的数据矩阵 (例如 data_list$Nplus), 应该是未标准化的原始计数或丰度
  initial_pcor_matrix,     # 在完整 data_matrix 上用 stabENG 得到的基础偏相关矩阵
  stabENG_params,          # 一个列表，包含传递给 stabENG 的所有参数 (除了data和labels)
  shared_otu_labels,       # OTU标签
  fold = 10,               # 交叉验证折数
  plot_cv_curve = FALSE,   # 是否绘制CV误差曲线
  plot_path_prefix = NULL  # 绘制曲线的路径前缀
) {
  # 确保data_matrix是数值型
  data_matrix_numeric <- apply(data_matrix, 2, as.numeric)
  if(!is.matrix(data_matrix_numeric)) data_matrix_numeric <- as.matrix(data_matrix_numeric)
  colnames(data_matrix_numeric) <- shared_otu_labels

  data_matrix_scaled <- scale(data_matrix_numeric) # 标准化数据
  # 处理全是0的列，scale后会变成NaN，将其设为0
  data_matrix_scaled[is.na(data_matrix_scaled)] <- 0


  n_samples <- nrow(data_matrix_scaled)
  p_vars <- ncol(data_matrix_scaled)

  if (n_samples < fold || n_samples < 2) {
      warning(paste0("Not enough samples (", n_samples, ") for ", fold, "-fold CV. Returning tau=0.05 or a quantile-based tau."))
      if (!is.null(initial_pcor_matrix) && sum(initial_pcor_matrix != 0) > 0) {
          non_zero_pcors <- abs(initial_pcor_matrix[upper.tri(initial_pcor_matrix) & initial_pcor_matrix!=0])
          if(length(non_zero_pcors) > 0) return(as.numeric(quantile(non_zero_pcors, 0.1, na.rm=TRUE)))
      }
      return(0.05)
  }

  part.list_cv_indices <- cv.part_internal(n_samples, fold) # 使用这里定义的cv.part_internal

  max_abs_pcor <- 0
  if (!is.null(initial_pcor_matrix) && length(initial_pcor_matrix[upper.tri(initial_pcor_matrix)]) > 0){
      max_abs_pcor <- max(abs(initial_pcor_matrix[upper.tri(initial_pcor_matrix)]), na.rm = TRUE)
  }

  if (max_abs_pcor == 0 || is.na(max_abs_pcor)) {
      warning("Initial partial correlation matrix has no edges or max_abs_pcor is NA. PCS screening might not be meaningful. Returning tau=0.")
      return(0)
  }
  taulist_cv <- seq(0.0001, min(1, max_abs_pcor + 0.05), length.out = 50)

  cv_losses_all_folds <- matrix(NA, nrow = length(taulist_cv), ncol = fold)

  cat("Starting PCS CV for threshold tau...\n")
  pb_cv <- txtProgressBar(min = 0, max = fold * length(taulist_cv), style = 3)
  progress_count = 0

  for (k_fold in 1:fold) {
    train_indices <- part.list_cv_indices$trainIndices[[k_fold]]
    test_indices <- part.list_cv_indices$testIndices[[k_fold]]

    # 确保训练集和测试集至少有2个样本，且训练集样本数大于变量数（对于OLS）或足够stabENG运行
    if(length(train_indices) < max(2, p_vars + 1) || length(test_indices) < 1){ # p_vars+1 for OLS with intercept, or just p_vars for no intercept
        warning(paste("Skipping fold", k_fold, "due to insufficient samples in train/test split for stabENG/lm."))
        cv_losses_all_folds[, k_fold] <- Inf # 标记此折无效
        progress_count = progress_count + length(taulist_cv)
        setTxtProgressBar(pb_cv, progress_count)
        next
    }

    x_train_unscaled <- data_matrix_numeric[train_indices, , drop = FALSE] # stabENG用原始数据
    x_train_scaled <- scale(x_train_unscaled) # Psi_Screen_Beta用标准化数据
    x_train_scaled[is.na(x_train_scaled)] <- 0


    x_test_scaled <- scale(data_matrix_numeric[test_indices, , drop = FALSE]) # 测试集也需要基于其自身进行标准化或基于训练集参数标准化
    # 更严谨：用训练集的中心和尺度参数标准化测试集
    train_means <- colMeans(x_train_unscaled, na.rm=TRUE)
    train_sds <- apply(x_train_unscaled, 2, sd, na.rm=TRUE)
    train_sds[train_sds == 0] <- 1 # 避免除以零
    x_test_properly_scaled <- sweep(data_matrix_numeric[test_indices, , drop = FALSE], 2, train_means, "-")
    x_test_properly_scaled <- sweep(x_test_properly_scaled, 2, train_sds, "/")
    x_test_properly_scaled[is.na(x_test_properly_scaled)] <- 0
    x_test_scaled <- x_test_properly_scaled


    stabENG_args_train <- c(
        list(data_list = list(current_group = as.data.frame(x_train_unscaled)), # stabENG用未标准化的
             labels = shared_otu_labels),
        stabENG_params
    )
    
    network_results_train <- tryCatch(do.call(stabENG::stabENG, stabENG_args_train), # 显式调用 stabENG::stabENG
                                      error = function(e) {
                                        # warning(paste("Error in stabENG for fold", k_fold, ":", e$message))
                                        NULL
                                      })

    if(is.null(network_results_train) || is.null(network_results_train$opt.fit.pcor$current_group)){
        # warning(paste("stabENG returned NULL pcor_matrix for fold", k_fold, "- setting losses to Inf for this fold."))
        cv_losses_all_folds[, k_fold] <- Inf
        progress_count = progress_count + length(taulist_cv)
        setTxtProgressBar(pb_cv, progress_count)
        next
    }
    pcor_matrix_train <- network_results_train$opt.fit.pcor$current_group
    colnames(pcor_matrix_train) <- rownames(pcor_matrix_train) <- shared_otu_labels


    for (j_tau in 1:length(taulist_cv)) {
      tau_current <- taulist_cv[j_tau]
      
       beta_train_for_tau_list <- tryCatch({
         Psi_Screen_Beta(x_train_scaled, pcor_matrix_train, tau_current) # Psi_Screen_Beta用标准化的x_train
       }, error = function(e) {
         # warning(paste("Error in Psi_Screen_Beta for fold", k_fold, "tau", tau_current, ":", e$message))
         NULL
       })

      if(!is.null(beta_train_for_tau_list)){
          beta_train_for_tau = beta_train_for_tau_list[[1]]
          prediction_errors_sq <- (x_test_scaled - x_test_scaled %*% beta_train_for_tau)^2
          cv_losses_all_folds[j_tau, k_fold] <- sum(colSums(prediction_errors_sq, na.rm = TRUE), na.rm = TRUE)
      } else {
          cv_losses_all_folds[j_tau, k_fold] <- Inf
      }
      progress_count = progress_count + 1
      setTxtProgressBar(pb_cv, progress_count)
    }
  }
  close(pb_cv)

  mean_cv_loss_per_tau <- rowMeans(cv_losses_all_folds, na.rm = TRUE)

  if(all(is.na(mean_cv_loss_per_tau)) || all(is.infinite(mean_cv_loss_per_tau)) || length(which.min(mean_cv_loss_per_tau))==0 ){
      warning("All mean CV losses are NA or Inf, or no minimum found. Returning a default tau.")
      if (!is.null(initial_pcor_matrix) && sum(initial_pcor_matrix != 0, na.rm=TRUE) > 0) {
          non_zero_pcors <- abs(initial_pcor_matrix[upper.tri(initial_pcor_matrix) & initial_pcor_matrix!=0])
          if(length(non_zero_pcors) > 0) {
              default_tau = as.numeric(quantile(non_zero_pcors, 0.1, na.rm=TRUE))
              cat(sprintf("\nUsing quantile-based default tau = %.4f\n", default_tau))
              return(default_tau)
          }
      }
      cat("\nUsing fixed default tau = 0.05\n")
      return(0.05)
  }

  optimal_tau_index <- which.min(mean_cv_loss_per_tau)
  optimal_tau <- taulist_cv[optimal_tau_index]

  if (plot_cv_curve && !is.null(plot_path_prefix)) {
    # ... (绘图代码不变) ...
  }
  cat(sprintf("\nOptimal PCS screening threshold tau = %.4f\n", optimal_tau))
  return(optimal_tau)
}


# --- 在你的主循环中使用 ---
# (你的脚本的其余部分，调用 pcs_cv_threshold)
# ... (主脚本中循环 for (i in timestamps) 的部分) ...

  if (perm_filter == TRUE) { # perm_filter 现在控制是否执行PCS筛选
    cat('\nPerforming PCS CV for threshold selection on Day',i,'\n')
    
    pcs_thresholds <- list()

    # --- PCS筛选 Nplus 组 ---
    if (!is.null(data_list_current_timestamp$Nplus) && 
        nrow(data_list_current_timestamp$Nplus) > 1 && # 需要至少2个样本
        ncol(data_list_current_timestamp$Nplus) > 1 && # 需要至少2个变量
        !is.null(network_pcor_raw[[i]]$Nplus)) {
      cat("  PCS for Nplus group:\n")
      pcs_thresholds$Nplus <- pcs_cv_threshold(
        data_matrix = data_list_current_timestamp$Nplus, # 未标准化的数据
        initial_pcor_matrix = network_pcor_raw[[i]]$Nplus,
        stabENG_params = stabENG_parameters, 
        shared_otu_labels = shared_otu,
        fold = 5, 
        plot_cv_curve = TRUE,
        plot_path_prefix = paste0(plot_path, "_Nplus")
      )
    } else {
      warning(paste("Skipping PCS for Nplus on Day", i, "due to insufficient data or missing initial pcor."))
      pcs_thresholds$Nplus <- 0.05 # 或其他合理的默认值，或者基于初始网络强度
    }

    # --- PCS筛选 Nminus 组 ---
    if (!is.null(data_list_current_timestamp$Nminus) && 
        nrow(data_list_current_timestamp$Nminus) > 1 &&
        ncol(data_list_current_timestamp$Nminus) > 1 &&
        !is.null(network_pcor_raw[[i]]$Nminus)) {
      cat("  PCS for Nminus group:\n")
      pcs_thresholds$Nminus <- pcs_cv_threshold(
        data_matrix = data_list_current_timestamp$Nminus, # 未标准化的数据
        initial_pcor_matrix = network_pcor_raw[[i]]$Nminus,
        stabENG_params = stabENG_parameters,
        shared_otu_labels = shared_otu,
        fold = 5,
        plot_cv_curve = TRUE,
        plot_path_prefix = paste0(plot_path, "_Nminus")
      )
    } else {
      warning(paste("Skipping PCS for Nminus on Day", i, "due to insufficient data or missing initial pcor."))
      pcs_thresholds$Nminus <- 0.05 # 或其他合理的默认值
    }

    # 应用PCS筛选得到的阈值 tau
    network_pcor[[i]]$Nplus <- network_pcor_raw[[i]]$Nplus
    network_pcor[[i]]$Nplus[abs(network_pcor[[i]]$Nplus) < pcs_thresholds$Nplus] <- 0
    diag(network_pcor[[i]]$Nplus) <- 0

    network_pcor[[i]]$Nminus <- network_pcor_raw[[i]]$Nminus
    network_pcor[[i]]$Nminus[abs(network_pcor[[i]]$Nminus) < pcs_thresholds$Nminus] <- 0
    diag(network_pcor[[i]]$Nminus) <- 0
    
    # (可选) 应用同样的逻辑过滤精度矩阵 network_list
    network_list[[i]]$Nplus <- network_list_raw[[i]]$Nplus
    if(!is.null(network_pcor_raw[[i]]$Nplus)) { # 确保 network_pcor_raw 存在
        indices_to_zero_Nplus <- abs(network_pcor_raw[[i]]$Nplus) < pcs_thresholds$Nplus
        network_list[[i]]$Nplus[indices_to_zero_Nplus] <- 0
    }
    diag(network_list[[i]]$Nplus) <- 0
    
    network_list[[i]]$Nminus <- network_list_raw[[i]]$Nminus
    if(!is.null(network_pcor_raw[[i]]$Nminus)) {
        indices_to_zero_Nminus <- abs(network_pcor_raw[[i]]$Nminus) < pcs_thresholds$Nminus
        network_list[[i]]$Nminus[indices_to_zero_Nminus] <- 0
    }
    diag(network_list[[i]]$Nminus) <- 0
    
    perm_threshold_Nplus <- pcs_thresholds$Nplus
    perm_threshold_Nminus <- pcs_thresholds$Nminus

  } else { 
    # ... ( perm_filter == FALSE 分支不变 ) ...
  }

# ... (你脚本的其余部分) ...
# cv.part_internal 的健壮性： 修改了 cv.part_internal 以处理当 k 大于样本数或 n/k < 1 的情况，并使其返回索引列表，这在处理不均等折时更安全。pcs_cv_threshold 内部也做了相应调整来使用这些列表。
# Psi_Screen_Beta 的健壮性：
# 确保输入 x 是矩阵。
# 确保 x 和 R 都有列名和行名，如果缺失则创建默认的。
# 在 lm 调用中，使用 paste0("", var_name, "", ...) 的方式构建公式，以更好地处理可能包含特殊字符的变量名。
# 增加了 tryCatch 来捕获 lm 可能发生的错误（例如，由于共线性或样本不足导致模型无法拟合）。
# 修正了将 coef(model_fit) 填回 beta_matrix_final 的逻辑，以正确匹配系数名称。
# pcs_cv_threshold 的健壮性：
# 在开始CV前检查样本量是否足够。
# 在每一折CV中，检查训练集和测试集样本是否足够运行 stabENG 和 lm。
# 标准化：stabENG 通常在内部处理数据的标准化或使用原始数据（取决于其设计，你的代码中似乎用的是原始数据列表）。Psi_Screen_Beta 则明确对输入 x 进行 scale。在 pcs_cv_threshold 中，x_train_unscaled 传递给 stabENG，而 x_train_scaled (对 x_train_unscaled 标准化得到) 传递给 Psi_Screen_Beta。测试集 x_test_scaled 也需要被正确标准化（理想情况是使用训练集的均值和标准差）。
# 增加了对 stabENG 和 Psi_Screen_Beta 返回 NULL 或计算失败的错误处理，将对应折的损失设为 Inf。
# 如果所有平均CV损失都是 NA 或 Inf，则返回一个默认的 tau 值（例如，0.05或基于初始网络稀疏性的分位数），并给出警告。
# 数据类型： 确保传递给 stabENG 的 data_list 中的数据是 data.frame，传递给 pcs_cv_threshold 的 data_matrix 以及内部传递给 Psi_Screen_Beta 的 x 是数值矩阵。
# 列名一致性： 整个流程中，数据矩阵、偏相关矩阵的列名和行名保持一致非常重要，尤其是在 Psi_Screen_Beta 中构建回归公式和匹配系数时。
# stabENG::stabENG： 在 do.call 中明确使用 stabENG::stabENG，以确保调用的是正确的包里的函数（如果你的环境中加载了多个同名函数）。
# 模拟数据部分的过滤： 在模拟循环中，如果 sim_perm_filter == TRUE，现在默认会为每个模拟数据集重新运行 pcs_cv_threshold（使用更少的 fold 以加速）。这是一个更严谨但耗时的做法。你可以根据需要切换回使用为真实数据计算的 tau。