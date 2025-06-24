# #%% requirements
# # !ebic.gamma: if too sparse change to 0.6
cat(R.version.string, "\nR Home:", R.home(), "\nR lib:", .libPaths(), "\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
library(devtools, quietly = TRUE)
library(BiocManager, quietly = TRUE)
if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
if (!requireNamespace("SPRING", quietly = TRUE)) devtools::install_github("GraceYoon/SPRING")
if (!requireNamespace("SpiecEasi", quietly = TRUE)) devtools::install_github("zdk123/SpiecEasi")
if (!requireNamespace("stabJGL", quietly = TRUE)) devtools::install_github("camiling/stabJGL")
if (!requireNamespace("NetCoMi", quietly = TRUE)) devtools::install_github("stefpeschel/NetCoMi", ref = "develop",repos = c("https://cloud.r-project.org/",BiocManager::repositories()))
if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
if (!requireNamespace("MLmetrics", quietly = TRUE)) install.packages("MLmetrics")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr") # Likely installed with tidyverse
if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse") # Installs a suite of packages
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if (!requireNamespace("foreach", quietly = TRUE)) install.packages("foreach")
if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel") # Usually base R, no need to install
if (!requireNamespace("mixedCCA", quietly = TRUE)) install.packages("mixedCCA")
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("qgraph", quietly = TRUE)) install.packages("qgraph")
if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
start_time <- Sys.time()
library(caret)
library(phyloseq)
library(SPRING)
library(SpiecEasi)
library(pROC)
library(MLmetrics)
library(ggplot2)
library(tidyr)
library(ggraph)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
source("//home//14720078//ProjectCode//Packages//stabENG.r")
source("//home//14720078//ProjectCode//Packages//MyENG.r")
source("//home//14720078//ProjectCode//Packages//MyPlot.r")
source("//home//14720078//ProjectCode//Packages//PCS_perm.r") # Source the PCS functions
# source("//home//14720078//ProjectCode//Packages//MyPerm.r") 
cat("Packages loaded successfully\n")

# Helper function for Matthews Correlation Coefficient (MCC)
# Define this if it's not available from sourced files (e.g., PCS.r)
calculate_mcc_for_metrics <- function(tp, tn, fp, fn) {
  tp <- as.numeric(tp)
  tn <- as.numeric(tn)
  fp <- as.numeric(fp)
  fn <- as.numeric(fn)
  
  numerator <- (tp * tn) - (fp * fn)
  
  d1 <- tp + fp
  d2 <- tp + fn
  d3 <- tn + fp
  d4 <- tn + fn
  
  if (d1 == 0 || d2 == 0 || d3 == 0 || d4 == 0) {
    return(0) # Conventionally, MCC is 0 if any part of the denominator product is 0
  }
  
  denominator <- sqrt(as.numeric(d1) * as.numeric(d2) * as.numeric(d3) * as.numeric(d4))
  
  if (denominator == 0) {
    return(0)
  }
  
  return(numerator / denominator)
}

# Metrics calculation function based on your provided structure, with corrections and robustness
calculate_metrics <- function(true_adj, sim_adj) {
  # Ensure matrices are binary (0 or 1)
  true_adj_bin <- (true_adj != 0) * 1
  sim_adj_bin  <- (sim_adj != 0) * 1
  
  # Flatten the matrices into vectors (upper triangular part, excluding diagonal)
  true_edges <- as.vector(true_adj_bin[upper.tri(true_adj_bin)])
  sim_edges  <- as.vector(sim_adj_bin[upper.tri(sim_adj_bin)])

  if (length(true_edges) == 0) { # No off-diagonal edges to compare
    return(list(TPR=NA, FPR=NA, Precision=NA, Recall=NA, F1=NA, AUC=NA, MCC=NA, TP=0,FP=0,TN=0,FN=0))
  }

  expected_levels <- c("0", "1")
  true_edges_factor <- factor(true_edges, levels = expected_levels)
  sim_edges_factor  <- factor(sim_edges,  levels = expected_levels)

  tp <- NA; tn <- NA; fp <- NA; fn <- NA
  tpr <- NA; fpr <- NA; precision <- NA; recall <- NA; f1_val <- NA; auc_val <- NA; mcc_val <- NA

  cm_obj <- tryCatch({
    caret::confusionMatrix(data = sim_edges_factor, reference = true_edges_factor, positive = "1")
  }, error = function(e) {
    warning(paste("caret::confusionMatrix failed:", conditionMessage(e), ". Manually calculating TP/FP/TN/FN."))
    NULL
  })

  if (!is.null(cm_obj) && !is.null(cm_obj$table) && all(dim(cm_obj$table) == c(2,2))) {
    tn <- as.numeric(cm_obj$table[1,1]) # Predicted 0, True 0
    fn <- as.numeric(cm_obj$table[1,2]) # Predicted 0, True 1 (Corrected from your original fp)
    fp <- as.numeric(cm_obj$table[2,1]) # Predicted 1, True 0 (Corrected from your original fn)
    tp <- as.numeric(cm_obj$table[2,2]) # Predicted 1, True 1
  } else {
    if(is.null(cm_obj) || is.null(cm_obj$table) || !all(dim(cm_obj$table) == c(2,2))) {
        warning("Confusion matrix from caret was not in the expected 2x2 format or failed. Manually calculating TP/FP/TN/FN.")
    }
    tp <- sum(sim_edges == 1 & true_edges == 1)
    fp <- sum(sim_edges == 1 & true_edges == 0)
    tn <- sum(sim_edges == 0 & true_edges == 0)
    fn <- sum(sim_edges == 0 & true_edges == 1)
  }
  
  # Calculate TPR, FPR, Precision, Recall
  # Handle division by zero: if denominator is 0, result is NA unless numerator is also 0 (then 0 or NA) or non-zero (then Inf or NA).
  # For ratios like TP/(TP+FN), if TP+FN=0, it implies TP=0 and FN=0.
  # If TP=0 and TP+FN=0, then TPR is undefined (0/0), often NA. Some might define as 0 or 1.
  if (!is.na(tp) && !is.na(fn)) {
    tpr <- ifelse((tp + fn) == 0, NA, tp / (tp + fn))
  }
  if (!is.na(fp) && !is.na(tn)) {
    fpr <- ifelse((fp + tn) == 0, NA, fp / (fp + tn))
  }
  if (!is.na(tp) && !is.na(fp)) {
    precision <- ifelse((tp + fp) == 0, NA, tp / (tp + fp))
  }
  recall <- tpr

  f1_val <- tryCatch({
    MLmetrics::F1_Score(y_true = true_edges_factor, y_pred = sim_edges_factor, positive = "1")
  }, error = function(e) {NA})
  
  if (is.na(f1_val) && !is.na(precision) && !is.na(recall)) { # Manual fallback for F1
      if (precision == 0 && recall == 0) f1_val <- 0
      else if ((precision + recall) > 0) f1_val <- 2 * (precision * recall) / (precision + recall)
      # else f1_val remains NA
  }


  if (length(unique(true_edges)) > 1) {
    roc_obj <- tryCatch({
      pROC::roc(response = true_edges, predictor = sim_edges, quiet = TRUE, levels=c(0,1), direction="<")
    }, error = function(e) {NULL})
    if(!is.null(roc_obj)) auc_val <- as.numeric(pROC::auc(roc_obj))
  } else { # True edges have only one class
    if (length(unique(sim_edges)) == 1 && unique(true_edges)[1] == unique(sim_edges)[1]) auc_val <- 1.0 # Perfect prediction of constant
    else if (length(unique(sim_edges)) == 1 && unique(true_edges)[1] != unique(sim_edges)[1]) auc_val <- 0.0 # Perfectly wrong
    else auc_val <- 0.5 # Mixed predictions for a constant true value
  }

  if (!any(is.na(c(tp, tn, fp, fn)))) {
    mcc_val <- calculate_mcc_for_metrics(tp, tn, fp, fn) # Use the helper defined above
  }

  return(list(TPR = tpr, FPR = fpr, Precision = precision, Recall = recall,
              F1 = f1_val, AUC = auc_val, MCC = mcc_val, TP=tp, FP=fp, TN=tn, FN=fn))
}


#load("//home//14720078//ProjectCode//DataImage//big1226_Days_network_results_Big_Days_Filtered.RData")
# start_time <- Sys.time() # This was the second definition of start_time, commented out.
rawdata <- readRDS("//home//14720078//ProjectCode//data//DavarData1_substrate_phyloseq_1226_final_filtered.rds")
otu_raw <- otu_table(rawdata)
sam_info <- as.data.frame(sample_data(rawdata))
otu_tax_full <- as.data.frame(tax_table(rawdata)) # Keep full tax table initially, rename to avoid conflict
rm(rawdata); gc() # Remove rawdata after extracting necessary components

otu_RA <- transform_sample_counts(otu_raw, function(x) x / sum(x) )
otu_RA <- filter_taxa(otu_RA, function(x) mean(x) > 1e-3, TRUE)
shared_otu_initial <- rownames(otu_RA) # Rename to avoid confusion before taxa_filtering
rm(otu_RA); gc() # Remove otu_RA

# Ensure shared_otu_initial has valid OTUs before subsetting otu_raw
if(length(shared_otu_initial) == 0) stop("No OTUs remaining after initial abundance filtering (1e-3).")
otu_Ab <- as.data.frame(t(otu_raw[shared_otu_initial,]))
rm(otu_raw); gc() # Remove otu_raw

taxa_filtering <- TRUE
n_sim <- 15
PCS_filter <- TRUE
sim_PCS_filter <- FALSE
# %%split otu_Ab by condition group
# otu_Ab_Nplus <- otu_Ab[rownames(otu_Ab) %in% rownames(sam_info[sam_info$growthCondition=="plusN",]),]
# otu_Ab_Nminus <- otu_Ab[rownames(otu_Ab) %in% rownames(sam_info[sam_info$growthCondition=="minusN",]),]
# data_list <- list(Nplus = otu_Ab_Nplus, Nminus = otu_Ab_Nminus)
# by Group and Days
timestamps <- as.character(sort(as.integer(levels(sam_info$Days))))
otu_Ab_Nplus_times <- list()
otu_Ab_Nminus_times <- list()
data_list_times <- list()
common_taxa <- list()
for (i in timestamps) {
  otu_Ab_Nplus_times[[i]] <- otu_Ab[rownames(otu_Ab) %in% 
    rownames(sam_info[sam_info$growthCondition == "plusN" & 
      sam_info$Days == i, ]), ]
  otu_Ab_Nminus_times[[i]] <- otu_Ab[rownames(otu_Ab) %in% 
    rownames(sam_info[sam_info$growthCondition == "minusN" & 
      sam_info$Days == i, ]), ]
  if (taxa_filtering == TRUE)
  {
    # Sum the number of zeros in each column and sort by the number of zeros
    zero_counts_Nplus <- colSums(otu_Ab_Nplus_times[[i]] == 0)
    zero_counts_Nminus <- colSums(otu_Ab_Nminus_times[[i]] == 0)
    
    # Determine the number of columns to select (top 75%)
    num_cols_to_select_Nplus <- floor(ncol(otu_Ab_Nplus_times[[i]]) * 0.75)
    num_cols_to_select_Nminus <- floor(ncol(otu_Ab_Nminus_times[[i]]) * 0.75)

    # Select the top 75% columns with the least zeros for each group
    sorted_Nplus <- names(sort(zero_counts_Nplus, decreasing = FALSE))[1:num_cols_to_select_Nplus]
    sorted_Nminus <- names(sort(zero_counts_Nminus, decreasing = FALSE))[1:num_cols_to_select_Nminus]
    
    # Use the intersecting column names for Nplus & Nminus
    common_taxa[[i]] <- intersect(sorted_Nplus, sorted_Nminus)
  }
  data_list_times[[i]] <- list(Nplus = otu_Ab_Nplus_times[[i]], 
                               Nminus = otu_Ab_Nminus_times[[i]])
}
# common_taxa <- Filter(Negate(is.null), common_taxa)
# common_taxa <- Filter(function(x) length(x) > 0, common_taxa)
if (taxa_filtering == TRUE)
{
  # After the first loop that populates common_taxa and before common_taxa_intersection is calculated
  # (Assuming this loop is before the main "for (i in timestamps)" for network calculation)
  # If common_taxa is large and only common_taxa_intersection is needed later:
  common_taxa_intersection <- Reduce(intersect, common_taxa)
  rm(common_taxa); gc() # Remove the full list common_taxa

  # Filter the data_list_times to only include the common taxa
  for (i in timestamps) {
    if (length(common_taxa_intersection) > 0) {
      if(!is.null(data_list_times[[i]]$Nplus)) data_list_times[[i]]$Nplus <- data_list_times[[i]]$Nplus[, colnames(data_list_times[[i]]$Nplus) %in% common_taxa_intersection, drop = FALSE]
      if(!is.null(data_list_times[[i]]$Nminus)) data_list_times[[i]]$Nminus <- data_list_times[[i]]$Nminus[, colnames(data_list_times[[i]]$Nminus) %in% common_taxa_intersection, drop = FALSE]
    } else {
      stop("No intersection found in common_taxa on Day ", i, ". Check the initial filtering criteria.")
    }
  }
  # Filter the shared_otu to only include the common taxa
  shared_otu <- shared_otu_initial[shared_otu_initial %in% common_taxa_intersection]
  # Filter the otu_tax to only include the common taxa
  otu_tax <- otu_tax_full[rownames(otu_tax_full) %in% common_taxa_intersection, ]
  rm(shared_otu_initial, otu_tax_full, common_taxa_intersection); gc() # Remove intermediate objects
} else {
  shared_otu <- shared_otu_initial
  otu_tax <- otu_tax_full
  rm(shared_otu_initial, otu_tax_full); gc()
}
# 
network_list_raw <- list()
network_list <- list()
network_pcor_raw <- list()
network_pcor <- list()
network_pcor_pcs_screened <- list() # 初始化用于存储PCS筛选结果的列表
confusion_results_df <- list()
Sim_list <- list()
Res_sim <- list()
lambda1 <- list()
lambda2 <- list()
# Function to synthesize scaled data based on the provided network
synthesize_scaled_data <- function(dat, net)
  {
    graph <- (net != 0)*1     # SpiecEasi::make_graph('cluster', dat$n_OTUs, dat$n_edge)
    attr(graph, "class") <- "graph"
    #cat("    Calling SpiecEasi::graph2prec...\n") # Add print statement
    Prec <- SpiecEasi::graph2prec(graph)
    Cor <- cov2cor(SpiecEasi::prec2cov(Prec))
    X <- SpiecEasi::synth_comm_from_counts(dat, mar = 2, distr = 'zinegbin', Sigma = Cor, n = nrow(dat))
    return(X)
  }

# Define stabENG parameters
stabENG_parameters <- list(
  var.thresh = 0.1,
  rep.num = 25,
  nlambda1 = 20,
  lambda1.min = 0.01,
  lambda1.max = 1,
  nlambda2 = 20,
  lambda2.min = 0,
  lambda2.max = 0.1,
  lambda2.init = 0.01,
  ebic.gamma = 0.5,
  nCores = if (Sys.getenv("SLURM_CPUS_PER_TASK") != "") as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")) else 2 # Added a default for nCores
)

# Modified permtest_freedman_lane function
permtest_freedman_lane <- function(data_list_current_group, # Data for the specific group to be permuted (e.g., data_list_current_timestamp$Nplus)
                                 original_other_group_data, # Original data for the *other* group (e.g., data_list_current_timestamp$Nminus)
                                 pcor_matrix_target_group,  # Raw pcor matrix for the target group (used for colnames)
                                 group_name_to_permute,     # Name of the group being permuted ("Nplus" or "Nminus")
                                 n_perm = 50,
                                 alpha = 0.1,
                                 stabENG_general_params,
                                 opt_lambda1_current,
                                 opt_lambda2_current,
                                 shared_labels_for_stabENG) { # OTU labels for stabENG call

  if (is.null(data_list_current_group) || nrow(data_list_current_group) == 0 || ncol(data_list_current_group) == 0) {
    warning(paste("No data provided for group to permute:", group_name_to_permute, ". Returning threshold 0."))
    return(0)
  }
  # original_other_group_data can be NULL if that group has no data, stabENG should handle this.

  residuals_current_group <- lm(as.matrix(data_list_current_group) ~ 1)$residuals
  col_means_current_group <- colMeans(as.matrix(data_list_current_group))
  
  null_dist <- numeric(0)
  
  other_group_name <- if(group_name_to_permute == "Nplus") "Nminus" else "Nplus"

  for (p_iter in 1:n_perm) {
    perm_residuals <- residuals_current_group[sample(nrow(residuals_current_group)), , drop = FALSE]
    perm_data_matrix_current_group <- sweep(perm_residuals, 2, col_means_current_group, "+")
    perm_data_df_current_group <- as.data.frame(perm_data_matrix_current_group)

    # Construct the Y list for stabENG: one group permuted, the other original
    Y_for_stabENG_perm <- list()
    if (group_name_to_permute == "Nplus") {
      Y_for_stabENG_perm$Nplus <- perm_data_df_current_group
      Y_for_stabENG_perm$Nminus <- original_other_group_data # This could be NULL if Nminus has no data
    } else { # group_name_to_permute == "Nminus"
      Y_for_stabENG_perm$Nplus <- original_other_group_data # This could be NULL if Nplus has no data
      Y_for_stabENG_perm$Nminus <- perm_data_df_current_group
    }
    
    # Ensure labels match the data being passed (shared_labels_for_stabENG should be appropriate for both groups)
    # stabENG will internally use labels for the matrices it creates.

    perm_network_pcor_matrix <- tryCatch({
        # Calling stabENG with fixed lambdas (by setting min=max=opt_lambda)
        # and reduced rep.num/nlambda for speed in permutation.
        # The critical part is that stabENG should use these fixed lambdas.
        stabENG(Y = Y_for_stabENG_perm,
                labels = shared_labels_for_stabENG, # Use the broader set of labels
                var.thresh = stabENG_general_params$var.thresh,
                rep.num = 1, # Minimal repetitions for permutation
                nlambda1 = 5, # Minimal search for fixed lambda
                lambda1.min = opt_lambda1_current, lambda1.max = opt_lambda1_current,
                nlambda2 = 5, # Minimal search for fixed lambda
                lambda2.min = opt_lambda2_current, lambda2.max = opt_lambda2_current,
                ebic.gamma = stabENG_general_params$ebic.gamma,
                nCores = stabENG_general_params$nCores,
                lambda1.sel = opt_lambda1_current, # Explicitly pass selected lambdas
                lambda2.sel = opt_lambda2_current  # Explicitly pass selected lambdas
        )$opt.fit.pcor[[group_name_to_permute]] # Extract pcor for the permuted group
    }, error = function(e){
        warning(paste("Error in stabENG during permutation for group", group_name_to_permute, "iter", p_iter, ":", conditionMessage(e)))
        NULL
    })

    if (!is.null(perm_network_pcor_matrix) && nrow(perm_network_pcor_matrix) > 0 && ncol(perm_network_pcor_matrix) > 0) {
      # Ensure the matrix is not empty and has the correct labels before taking upper.tri
      if(all(shared_labels_for_stabENG %in% colnames(perm_network_pcor_matrix))) {
         # Reorder if necessary, though stabENG should return with consistent labels
         # perm_network_pcor_matrix_ordered <- perm_network_pcor_matrix[shared_labels_for_stabENG, shared_labels_for_stabENG]
         null_dist <- c(null_dist, abs(perm_network_pcor_matrix[upper.tri(perm_network_pcor_matrix)]))
      } else {
         warning(paste("Label mismatch in permuted pcor matrix for group", group_name_to_permute, "iter", p_iter))
      }
    }
  }
  
  if (length(null_dist) == 0) {
    warning(paste("No null pcor values generated for group", group_name_to_permute, ". Returning threshold 0."))
    return(0)
  }
  
  threshold <- quantile(null_dist, 1 - alpha, na.rm = TRUE)
  return(threshold)
}

# %% calculate network
# timestamps <- timestamps[1] # For testing a single timestamp
for (i in timestamps){
  plot_dir <- file.path("Plots", "BigDataDaysFreedmanLane") # Changed plot dir
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  plot_path <- file.path(plot_dir, paste0("Day_", i))
  data_list_current_timestamp <- data_list_times[[i]]
  cat('Calculating network on day ',i,'\n')
  
  stabENG_args_real_data <- c(
    list(Y = data_list_current_timestamp, labels = shared_otu),
    stabENG_parameters # General stabENG parameters
  )
  network_results <- do.call(stabENG, stabENG_args_real_data)
  cat('Network calculation on Day', i, 'completed.\n')
  
  network_list_raw[[i]]$Nplus <- network_results$opt.fit$Nplus
  network_list_raw[[i]]$Nminus <- network_results$opt.fit$Nminus
  if(!is.null(network_list_raw[[i]]$Nplus)) diag(network_list_raw[[i]]$Nplus) <- 0
  if(!is.null(network_list_raw[[i]]$Nminus)) diag(network_list_raw[[i]]$Nminus) <- 0
  
  network_pcor_raw[[i]]$Nplus <- network_results$opt.fit.pcor$Nplus
  network_pcor_raw[[i]]$Nminus <- network_results$opt.fit.pcor$Nminus
  if(!is.null(network_pcor_raw[[i]]$Nplus)) diag(network_pcor_raw[[i]]$Nplus) <- 0
  if(!is.null(network_pcor_raw[[i]]$Nminus)) diag(network_pcor_raw[[i]]$Nminus) <- 0
  
  # Store optimized lambdas for the current day i
  current_opt_lambda1 <- network_results$opt.lambda1
  current_opt_lambda2 <- network_results$opt.lambda2
  lambda1[[i]] <- current_opt_lambda1 # Storing for record
  lambda2[[i]] <- current_opt_lambda2 # Storing for record
  
  rm(network_results); gc() # Removed stabENG_args_real_data from rm as it's small

  # --- BEGIN PCS Screening Section ---
  # cat('\nPerforming PCS CV (Lasso based method, returning PCS paper tau) screening for network on Day',i,'\n')
  # # Define PCS parameters
  # n_data_permutations_pcs <- 3 # PCS CV 中的数据重排次数
  # tau_pcs_nplus <- NA # 使用 NA 进行初始化
  # tau_pcs_nminus <- NA # 使用 NA 进行初始化

  # # 为 PCS 调用准备 stabENG_params_list
  # # pcs_cv_threshold_stabENG_lasso 内部的 CV 部分会使用 stabENG_parameters 中的 nlambda, lambda.min/max 等
  # # 而其末尾调用的 pcor_screen_pcs (其返回值是 pcs_cv_threshold_stabENG_lasso 的最终返回值)
  # # 需要 opt.lambda1 和 opt.lambda2
  # stabENG_params_for_pcs <- stabENG_parameters
  # stabENG_params_for_pcs$opt.lambda1 <- current_opt_lambda1
  # stabENG_params_for_pcs$opt.lambda2 <- current_opt_lambda2

  # # PCS screening for Nplus
  # if (!is.null(data_list_current_timestamp$Nplus) && nrow(data_list_current_timestamp$Nplus) > 0 && 
  #     !is.null(network_pcor_raw[[i]]$Nplus) && ncol(network_pcor_raw[[i]]$Nplus) > 0 && nrow(network_pcor_raw[[i]]$Nplus) > 0) {
      
  #     tau_pcs_nplus <- pcs_cv_threshold_stabENG_lasso_perm(
  #         data_list_unscaled_full = data_list_current_timestamp,
  #         initial_pcor_matrix_target = network_pcor_raw[[i]]$Nplus,
  #         group_name_target = "Nplus",
  #         other_group_name = "Nminus",
  #         stabENG_params_list = stabENG_params_for_pcs,
  #         otu_labels = shared_otu, # 确保 shared_otu 与 network_pcor_raw[[i]]$Nplus 的维度匹配
  #         fold = nrow(data_list_current_timestamp$Nplus), # LOOCV
  #         plot_cv_curve = TRUE, 
  #         plot_path_prefix = plot_path, # 使用 plot_path 作为前缀
  #         n_data_permutations = n_data_permutations_pcs
  #     )
  # } else {
  #     warning(paste0("Skipping PCS CV for Nplus on Day ", i, " due to no data or no raw pcor matrix."))
  #     # tau_pcs_nplus 保持 NA 或可设为默认值如 0.05，但 NA 更能反映未执行
  # }

  # # PCS screening for Nminus
  # if (!is.null(data_list_current_timestamp$Nminus) && nrow(data_list_current_timestamp$Nminus) > 0 && 
  #     !is.null(network_pcor_raw[[i]]$Nminus) && ncol(network_pcor_raw[[i]]$Nminus) > 0 && nrow(network_pcor_raw[[i]]$Nminus) > 0) {
      
  #     tau_pcs_nminus <- pcs_cv_threshold_stabENG_lasso_perm(
  #         data_list_unscaled_full = data_list_current_timestamp,
  #         initial_pcor_matrix_target = network_pcor_raw[[i]]$Nminus,
  #         group_name_target = "Nminus",
  #         other_group_name = "Nplus",
  #         stabENG_params_list = stabENG_params_for_pcs,
  #         otu_labels = shared_otu, # 确保 shared_otu 与 network_pcor_raw[[i]]$Nminus 的维度匹配
  #         fold = nrow(data_list_current_timestamp$Nminus), # LOOCV
  #         plot_cv_curve = TRUE,
  #         plot_path_prefix = plot_path, # 使用 plot_path 作为前缀
  #         n_data_permutations = n_data_permutations_pcs
  #     )
  # } else {
  #     warning(paste0("Skipping PCS CV for Nminus on Day ", i, " due to no data or no raw pcor matrix."))
  #     # tau_pcs_nminus 保持 NA
  # }
  
  # cat(sprintf("\nPCS method optimal tau for Day %s:\n  Nplus: %.4f\n  Nminus: %.4f\n",
  #             i, ifelse(is.na(tau_pcs_nplus), NA, tau_pcs_nplus), ifelse(is.na(tau_pcs_nminus), NA, tau_pcs_nminus)))

  # # Apply PCS screening to network_pcor_raw
  # network_pcor_pcs_screened[[i]] <- list()
  # network_pcor_pcs_screened[[i]]$Nplus <- if(!is.null(network_pcor_raw[[i]]$Nplus)) network_pcor_raw[[i]]$Nplus else NULL
  # network_pcor_pcs_screened[[i]]$Nminus <- if(!is.null(network_pcor_raw[[i]]$Nminus)) network_pcor_raw[[i]]$Nminus else NULL

  # if(!is.null(network_pcor_pcs_screened[[i]]$Nplus) && !is.na(tau_pcs_nplus)) {
  #     network_pcor_pcs_screened[[i]]$Nplus[abs(network_pcor_pcs_screened[[i]]$Nplus) < tau_pcs_nplus] <- 0
  #     diag(network_pcor_pcs_screened[[i]]$Nplus) <- 0 
  # }
  # if(!is.null(network_pcor_pcs_screened[[i]]$Nminus) && !is.na(tau_pcs_nminus)) {
  #     network_pcor_pcs_screened[[i]]$Nminus[abs(network_pcor_pcs_screened[[i]]$Nminus) < tau_pcs_nminus] <- 0
  #     diag(network_pcor_pcs_screened[[i]]$Nminus) <- 0
  # }
  
  # cat('PCS Screened edge number on Day', i, 
  #     'Nplus (pcor):', if(!is.null(network_pcor_pcs_screened[[i]]$Nplus)) sum(network_pcor_pcs_screened[[i]]$Nplus[upper.tri(network_pcor_pcs_screened[[i]]$Nplus)] != 0) else 0,
  #     'Nminus (pcor):', if(!is.null(network_pcor_pcs_screened[[i]]$Nminus)) sum(network_pcor_pcs_screened[[i]]$Nminus[upper.tri(network_pcor_pcs_screened[[i]]$Nminus)] != 0) else 0, '\n')
  # # --- END PCS Screening Section ---

  cat('\nPerforming Freedman-Lane permutation test for network on Day',i,'\n')
  perm_thresholds_values <- list()

  # shared_otu is the set of labels used for the initial stabENG call for this day
  # It should be appropriate for permtest_freedman_lane's stabENG call too.

  # Permutation Test for Nplus group
  if (!is.null(data_list_current_timestamp$Nplus) && nrow(data_list_current_timestamp$Nplus) > 0 && !is.null(network_pcor_raw[[i]]$Nplus)) {
      perm_thresholds_values$Nplus <- permtest_freedman_lane(
        data_list_current_group = data_list_current_timestamp$Nplus,
        original_other_group_data = data_list_current_timestamp$Nminus, # Pass original Nminus data
        pcor_matrix_target_group = network_pcor_raw[[i]]$Nplus, # For colnames reference if needed, not directly used now
        group_name_to_permute = "Nplus",
        n_perm = 50, 
        alpha = 0.1,
        stabENG_general_params = stabENG_parameters,
        opt_lambda1_current = current_opt_lambda1,
        opt_lambda2_current = current_opt_lambda2,
        shared_labels_for_stabENG = shared_otu # Labels used in the initial stabENG call
      )
  } else {
      warning(paste0("Skipping Permutation Test for Nplus on Day ", i, " due to no data or no raw pcor matrix."))
      perm_thresholds_values$Nplus <- 0.05 # Default or handle as error
  }

  # Permutation Test for Nminus group
  if (!is.null(data_list_current_timestamp$Nminus) && nrow(data_list_current_timestamp$Nminus) > 0 && !is.null(network_pcor_raw[[i]]$Nminus)) {
      perm_thresholds_values$Nminus <- permtest_freedman_lane(
        data_list_current_group = data_list_current_timestamp$Nminus,
        original_other_group_data = data_list_current_timestamp$Nplus, # Pass original Nplus data
        pcor_matrix_target_group = network_pcor_raw[[i]]$Nminus, # For colnames reference
        group_name_to_permute = "Nminus",
        n_perm = 50, 
        alpha = 0.1,
        stabENG_general_params = stabENG_parameters,
        opt_lambda1_current = current_opt_lambda1,
        opt_lambda2_current = current_opt_lambda2,
        shared_labels_for_stabENG = shared_otu # Labels used in the initial stabENG call
      )
  } else {
      warning(paste0("Skipping Permutation Test for Nminus on Day ", i, " due to no data or no raw pcor matrix."))
      perm_thresholds_values$Nminus <- 0.05 # Default or handle as error
  }
  
  cat(sprintf("\nPermutation test thresholds for Day %s:\n  Nplus: %.4f\n  Nminus: %.4f\n",
              i, perm_thresholds_values$Nplus, perm_thresholds_values$Nminus))

  # Apply thresholds ONLY to pcor matrices
  network_pcor[[i]]$Nplus <- network_pcor_raw[[i]]$Nplus
  network_pcor[[i]]$Nminus <- network_pcor_raw[[i]]$Nminus
  if(!is.null(network_pcor[[i]]$Nplus) && !is.na(perm_thresholds_values$Nplus)) {
    network_pcor[[i]]$Nplus[abs(network_pcor[[i]]$Nplus) < perm_thresholds_values$Nplus] <- 0
  }
  if(!is.null(network_pcor[[i]]$Nminus) && !is.na(perm_thresholds_values$Nminus)) {
    network_pcor[[i]]$Nminus[abs(network_pcor[[i]]$Nminus) < perm_thresholds_values$Nminus] <- 0
  }
  network_list[[i]]$Nplus <- network_list_raw[[i]]$Nplus
  network_list[[i]]$Nminus <- network_list_raw[[i]]$Nminus
  if(!is.null(network_list[[i]]$Nplus) && !is.na(perm_thresholds_values$Nplus)) {
    network_list[[i]]$Nplus[abs(network_list[[i]]$Nplus) < perm_thresholds_values$Nplus] <- 0
  }
  if(!is.null(network_list[[i]]$Nminus) && !is.na(perm_thresholds_values$Nminus)) {
    network_list[[i]]$Nminus[abs(network_list[[i]]$Nminus) < perm_thresholds_values$Nminus] <- 0
  }
  
  cat('Filtered OTU number on Day', i, 'Nplus (pcor):', if(!is.null(network_pcor[[i]]$Nplus)) sum(network_pcor[[i]]$Nplus != 0)/2 else 0, # Counting non-zero edges in pcor
      'Nminus (pcor):', if(!is.null(network_pcor[[i]]$Nminus)) sum(network_pcor[[i]]$Nminus != 0)/2 else 0, '\n') # Counting non-zero edges in pcor
  cat('Raw OTU number on Day', i, 'Nplus (precision):', if(!is.null(network_list[[i]]$Nplus)) sum(network_list[[i]]$Nplus != 0)/2 else 0, # Counting non-zero edges in precision
      'Nminus (precision):', if(!is.null(network_list[[i]]$Nminus)) sum(network_list[[i]]$Nminus != 0)/2 else 0, '\n')


  # --- OTU Subsetting Logic ---
  # This logic now needs to decide which matrix (pcor or precision) drives the OTU subsetting.
  # Previously, it was based on network_list (precision matrix after thresholding).
  # If we want the simulation to be based on the sparse pcor structure,
  # then OTU subsetting should be based on network_pcor.
  # If the precision matrix from stabENG (network_list_raw) is considered the "true" structure for simulation,
  # then subsetting should be based on network_list_raw.

  # Option 1: OTU Subsetting based on the filtered pcor matrix (network_pcor)
  # This makes sense if the goal is to simulate based on the statistically filtered connections.
  
  # Refine based on filtered network_pcor for Nplus
    if (!is.null(network_pcor[[i]]$Nplus) && nrow(network_pcor[[i]]$Nplus) > 0 && ncol(network_pcor[[i]]$Nplus) > 0) {
        rows_to_keep_nplus_pcor <- apply(network_pcor[[i]]$Nplus, 1, function(x) any(x != 0))
        cols_to_keep_nplus_pcor <- apply(network_pcor[[i]]$Nplus, 2, function(x) any(x != 0))
        common_otus_pcor_nplus <- rownames(network_pcor[[i]]$Nplus)[rows_to_keep_nplus_pcor & cols_to_keep_nplus_pcor]
        if (length(common_otus_pcor_nplus) > 0) {
            # Keep the pcor matrix subsetted
            network_pcor[[i]]$Nplus <- network_pcor[[i]]$Nplus[common_otus_pcor_nplus, common_otus_pcor_nplus, drop = FALSE]
            # Also subset the corresponding raw precision matrix to these OTUs
            if (!is.null(network_list[[i]]$Nplus)) {
                 network_list[[i]]$Nplus <- network_list_raw[[i]]$Nplus[rownames(network_list_raw[[i]]$Nplus) %in% common_otus_pcor_nplus,
                                                                        colnames(network_list_raw[[i]]$Nplus) %in% common_otus_pcor_nplus, drop = FALSE]
            }
        } else {
            network_pcor[[i]]$Nplus <- matrix(0,0,0, dimnames=list(NULL,NULL))
            if (!is.null(network_list[[i]]$Nplus)) network_list[[i]]$Nplus <- matrix(0,0,0, dimnames=list(NULL,NULL))
            warning(paste0("No common OTUs found in Nplus pcor after Permutation Test filtering on Day ", i, ". Setting Nplus to empty matrix."))
        }
    } else {
        network_pcor[[i]]$Nplus <- matrix(0,0,0, dimnames=list(NULL,NULL))
        if (!is.null(network_list[[i]]$Nplus)) network_list[[i]]$Nplus <- matrix(0,0,0, dimnames=list(NULL,NULL))
    }

    # Refine based on filtered network_pcor for Nminus
    if (!is.null(network_pcor[[i]]$Nminus) && nrow(network_pcor[[i]]$Nminus) > 0 && ncol(network_pcor[[i]]$Nminus) > 0) {
        rows_to_keep_nminus_pcor <- apply(network_pcor[[i]]$Nminus, 1, function(x) any(x != 0))
        cols_to_keep_nminus_pcor <- apply(network_pcor[[i]]$Nminus, 2, function(x) any(x != 0))
        common_otus_pcor_nminus <- rownames(network_pcor[[i]]$Nminus)[rows_to_keep_nminus_pcor & cols_to_keep_nminus_pcor]
        if (length(common_otus_pcor_nminus) > 0) {
            network_pcor[[i]]$Nminus <- network_pcor[[i]]$Nminus[common_otus_pcor_nminus, common_otus_pcor_nminus, drop = FALSE]
            if (!is.null(network_list[[i]]$Nminus)) {
                network_list[[i]]$Nminus <- network_list_raw[[i]]$Nminus[rownames(network_list_raw[[i]]$Nminus) %in% common_otus_pcor_nminus,
                                                                         colnames(network_list_raw[[i]]$Nminus) %in% common_otus_pcor_nminus, drop = FALSE]
            }
        } else {
            network_pcor[[i]]$Nminus <- matrix(0,0,0, dimnames=list(NULL,NULL))
            if (!is.null(network_list[[i]]$Nminus)) network_list[[i]]$Nminus <- matrix(0,0,0, dimnames=list(NULL,NULL))
            warning(paste0("No common OTUs found in Nminus pcor after Permutation Test filtering on Day ", i, ". Setting Nminus to empty matrix."))
        }
    } else {
        network_pcor[[i]]$Nminus <- matrix(0,0,0, dimnames=list(NULL,NULL))
        if (!is.null(network_list[[i]]$Nminus)) network_list[[i]]$Nminus <- matrix(0,0,0, dimnames=list(NULL,NULL))
    }

    cat('OTU number after pcor-based subsetting on Day', i, 
        'Nplus (pcor):', if(!is.null(network_pcor[[i]]$Nplus)) nrow(network_pcor[[i]]$Nplus) else 0, 
        'Nminus (pcor):', if(!is.null(network_pcor[[i]]$Nminus)) nrow(network_pcor[[i]]$Nminus) else 0, '\n')
    cat('Corresponding OTU number on Day', i, 
        'Nplus (precision):', if(!is.null(network_list[[i]]$Nplus)) nrow(network_list[[i]]$Nplus) else 0, 
        'Nminus (precision):', if(!is.null(network_list[[i]]$Nminus)) nrow(network_list[[i]]$Nminus) else 0, '\n')


    otus_in_nplus_final <- if(!is.null(rownames(network_pcor[[i]]$Nplus))) rownames(network_pcor[[i]]$Nplus) else character(0)
    otus_in_nminus_final <- if(!is.null(rownames(network_pcor[[i]]$Nminus))) rownames(network_pcor[[i]]$Nminus) else character(0)
    
    shared_otu_perm_filtered <- character(0) 
    if (length(otus_in_nplus_final) > 0 && length(otus_in_nminus_final) > 0) {
        shared_otu_perm_filtered <- intersect(otus_in_nplus_final, otus_in_nminus_final)
    } else if (length(otus_in_nplus_final) > 0) {
        shared_otu_perm_filtered <- otus_in_nplus_final
    } else if (length(otus_in_nminus_final) > 0) {
        shared_otu_perm_filtered <- otus_in_nminus_final
    }
    cat('Shared OTU number for simulation basis (shared_otu_perm_filtered) on Day', i, ':', length(shared_otu_perm_filtered), '\n')

    # Now, ensure both network_pcor and network_list are subsetted to these shared_otu_perm_filtered
    if (length(shared_otu_perm_filtered) > 0) {
        if (!is.null(network_pcor[[i]]$Nplus) && nrow(network_pcor[[i]]$Nplus) > 0) {
             network_pcor[[i]]$Nplus <- network_pcor[[i]]$Nplus[rownames(network_pcor[[i]]$Nplus) %in% shared_otu_perm_filtered, 
                                                                colnames(network_pcor[[i]]$Nplus) %in% shared_otu_perm_filtered, drop = FALSE]
        }
        if (!is.null(network_list[[i]]$Nplus) && nrow(network_list[[i]]$Nplus) > 0) { # network_list was already subsetted if pcor had common OTUs
             network_list[[i]]$Nplus <- network_list[[i]]$Nplus[rownames(network_list[[i]]$Nplus) %in% shared_otu_perm_filtered, 
                                                                colnames(network_list[[i]]$Nplus) %in% shared_otu_perm_filtered, drop = FALSE]
        }

        if (!is.null(network_pcor[[i]]$Nminus) && nrow(network_pcor[[i]]$Nminus) > 0) {
             network_pcor[[i]]$Nminus <- network_pcor[[i]]$Nminus[rownames(network_pcor[[i]]$Nminus) %in% shared_otu_perm_filtered, 
                                                                  colnames(network_pcor[[i]]$Nminus) %in% shared_otu_perm_filtered, drop = FALSE]
        }
        if (!is.null(network_list[[i]]$Nminus) && nrow(network_list[[i]]$Nminus) > 0) {
             network_list[[i]]$Nminus <- network_list[[i]]$Nminus[rownames(network_list[[i]]$Nminus) %in% shared_otu_perm_filtered, 
                                                                  colnames(network_list[[i]]$Nminus) %in% shared_otu_perm_filtered, drop = FALSE]
        }
    } else {
        warning(paste0("No shared OTUs after pcor-based filtering for Day ", i, ". Simulation might be problematic."))
        if(!is.null(network_pcor[[i]]$Nplus)) network_pcor[[i]]$Nplus <- matrix(0,0,0, dimnames=list(NULL,NULL))
        if(!is.null(network_list[[i]]$Nplus)) network_list[[i]]$Nplus <- matrix(0,0,0, dimnames=list(NULL,NULL))
        if(!is.null(network_pcor[[i]]$Nminus)) network_pcor[[i]]$Nminus <- matrix(0,0,0, dimnames=list(NULL,NULL))
        if(!is.null(network_list[[i]]$Nminus)) network_list[[i]]$Nminus <- matrix(0,0,0, dimnames=list(NULL,NULL))
    }
    
  # --- End OTU Subsetting Logic ---

  # %% Plot network on Phylum level
  # Ensure Phylum_groups are defined based on the OTUs present in the filtered network_pcor
  # This requires otu_tax to be available and correctly subsetted if shared_otu_perm_filtered is used.
  # Assuming otu_tax is filtered globally to 'shared_otu' (which might be pre-permutation filtering)
  # For plotting, we need to use OTUs that are actually in the plotted matrix.
  
  plot_otus_nplus <- if(!is.null(network_pcor[[i]]$Nplus) && nrow(network_pcor[[i]]$Nplus) > 0) rownames(network_pcor[[i]]$Nplus) else character(0)
  plot_otus_nminus <- if(!is.null(network_pcor[[i]]$Nminus) && nrow(network_pcor[[i]]$Nminus) > 0) rownames(network_pcor[[i]]$Nminus) else character(0)

  if (length(plot_otus_nplus) > 0) {
    Phylum_groups_nplus <- as.factor(otu_tax[plot_otus_nplus,"Phylum"])
    png(filename=paste0(plot_path,"_network_Nplus_Phylum_PermTest_Filtered_vsized.png"))
    qgraph::qgraph(network_pcor[[i]]$Nplus, 
      layout = "circle",
      edge.color = ifelse(network_pcor[[i]]$Nplus > 0, "#2b732b", "red"),
      title = "PermTest Filtered Network Nplus by Phylum",
      vsize = 2.5,
      groups = Phylum_groups_nplus)
    dev.off()
  } else { cat("Skipping Nplus network plot for Day", i, "as no OTUs remain after filtering.\n") }
  
  if (length(plot_otus_nminus) > 0) {
    Phylum_groups_nminus <- as.factor(otu_tax[plot_otus_nminus,"Phylum"])
    png(filename=paste0(plot_path,"_network_Nminus_Phylum_PermTest_Filtered_vsized.png"))
    qgraph::qgraph(network_pcor[[i]]$Nminus, 
      layout = "circle",
      edge.color = ifelse(network_pcor[[i]]$Nminus > 0, "#2b732b", "red"),
      title = "PermTest Filtered Network Nminus by Phylum",
      vsize = 2.5,
      groups = Phylum_groups_nminus)
    dev.off()
  } else { cat("Skipping Nminus network plot for Day", i, "as no OTUs remain after filtering.\n") }
  
  # %%Visualize Edge weights (similar to before)
  # ... (density plot code remains the same, using the filtered network_pcor) ...
  cor_values_Nplus <- if(!is.null(network_pcor[[i]]$Nplus)) as.vector(network_pcor[[i]]$Nplus) else numeric(0)
  cor_values_Nminus <- if(!is.null(network_pcor[[i]]$Nminus)) as.vector(network_pcor[[i]]$Nminus) else numeric(0)
  cor_values_Nplus_nonzero <- cor_values_Nplus[cor_values_Nplus != 0]
  cor_values_Nminus_nonzero <- cor_values_Nminus[cor_values_Nminus != 0]

  if(length(cor_values_Nplus_nonzero) > 0 || length(cor_values_Nminus_nonzero) > 0) {
    ggplot(data.frame(Type = c(rep("N+", length(cor_values_Nplus_nonzero)), rep("N-", length(cor_values_Nminus_nonzero))), 
      Value = c(cor_values_Nplus_nonzero, cor_values_Nminus_nonzero)), 
      aes(x = Value, color = Type)) + 
      geom_density() + 
      labs(x = "Correlation Value (PermTest Filtered)", y = "Density")
    ggsave(filename=paste0(plot_path,"_correlation_distribution_PermTest_Filtered.png"))
  } else {
    cat("Skipping correlation distribution plot for Day", i, "as no non-zero correlations remain.\n")
  }
  
  #%% Simulation Part
  set.seed(10010)
  cat('Synthesize simulation data on day ',i,'\n')

  current_sim_otus <- if (length(shared_otu_perm_filtered) > 0) shared_otu_perm_filtered else shared_otu
  
  dat_sim_Nplus <- if (!is.null(data_list_current_timestamp$Nplus) && ncol(data_list_current_timestamp$Nplus) > 0) {
                       data_list_current_timestamp$Nplus[, colnames(data_list_current_timestamp$Nplus) %in% current_sim_otus, drop = FALSE]
                   } else { NULL }
  dat_sim_Nminus <- if (!is.null(data_list_current_timestamp$Nminus) && ncol(data_list_current_timestamp$Nminus) > 0) {
                        data_list_current_timestamp$Nminus[, colnames(data_list_current_timestamp$Nminus) %in% current_sim_otus, drop = FALSE]
                    } else { NULL }

  true_net_sim_Nplus <- if (!is.null(network_list[[i]]$Nplus) && ncol(network_list[[i]]$Nplus) > 0 && nrow(network_list[[i]]$Nplus) > 0) {
                            network_list[[i]]$Nplus[rownames(network_list[[i]]$Nplus) %in% current_sim_otus, 
                                                    colnames(network_list[[i]]$Nplus) %in% current_sim_otus, drop = FALSE]
                        } else { matrix(0, length(current_sim_otus), length(current_sim_otus), dimnames=list(current_sim_otus, current_sim_otus)) }
  true_net_sim_Nminus <- if (!is.null(network_list[[i]]$Nminus) && ncol(network_list[[i]]$Nminus) > 0 && nrow(network_list[[i]]$Nminus) > 0) {
                             network_list[[i]]$Nminus[rownames(network_list[[i]]$Nminus) %in% current_sim_otus, 
                                                      colnames(network_list[[i]]$Nminus) %in% current_sim_otus, drop = FALSE]
                         } else { matrix(0, length(current_sim_otus), length(current_sim_otus), dimnames=list(current_sim_otus, current_sim_otus)) }

  # Ensure proper dimensions
  if(length(current_sim_otus) > 0){
      if(nrow(true_net_sim_Nplus) == 0 || ncol(true_net_sim_Nplus) == 0) {
          true_net_sim_Nplus <- matrix(0, length(current_sim_otus), length(current_sim_otus), dimnames=list(current_sim_otus, current_sim_otus))
      }
      if(nrow(true_net_sim_Nminus) == 0 || ncol(true_net_sim_Nminus) == 0) {
          true_net_sim_Nminus <- matrix(0, length(current_sim_otus), length(current_sim_otus), dimnames=list(current_sim_otus, current_sim_otus))
      }
  }

  # Generate simulation data
  if (length(current_sim_otus) > 0 && !is.null(dat_sim_Nplus) && !is.null(dat_sim_Nminus) && 
      nrow(dat_sim_Nplus) > 0 && ncol(dat_sim_Nplus) > 0 && 
      nrow(dat_sim_Nminus) > 0 && ncol(dat_sim_Nminus) > 0) { 
      for (j in 1:n_sim) {
          Sim_list[[i]][[j]] <- list(
            Nplus = synthesize_scaled_data(dat_sim_Nplus, true_net_sim_Nplus),
            Nminus = synthesize_scaled_data(dat_sim_Nminus, true_net_sim_Nminus)
          )
      }
      rm(dat_sim_Nplus, dat_sim_Nminus);
  } else {
      warning(paste0("Skipping simulation data synthesis for Day ", i, " due to insufficient data."))
      Sim_list[[i]] <- vector("list", n_sim)
      for (j in 1:n_sim) { 
          Sim_list[[i]][[j]] <- list(Nplus = matrix(0,0,0), Nminus = matrix(0,0,0))
      }
  }
  
  cat('Calculate simulation data network on day ',i,' start.\n')

  Res_sim[[i]] <- vector("list", n_sim)
  sim_perm_thresholds_list <- vector("list", n_sim) # Store simulation-specific permutation thresholds

  for (j in 1:n_sim) {
      if ( (is.matrix(Sim_list[[i]][[j]]$Nplus) && nrow(Sim_list[[i]][[j]]$Nplus) > 0 && ncol(Sim_list[[i]][[j]]$Nplus) > 0) || 
           (is.matrix(Sim_list[[i]][[j]]$Nminus) && nrow(Sim_list[[i]][[j]]$Nminus) > 0 && ncol(Sim_list[[i]][[j]]$Nminus) > 0) ) {
          
          # Get labels for stabENG
          sim_data_labels <- character(0)
          if(ncol(Sim_list[[i]][[j]]$Nplus) > 0) sim_data_labels <- colnames(Sim_list[[i]][[j]]$Nplus)
          else if(ncol(Sim_list[[i]][[j]]$Nminus) > 0) sim_data_labels <- colnames(Sim_list[[i]][[j]]$Nminus)
          
          if(length(sim_data_labels) == 0 && length(current_sim_otus) > 0) sim_data_labels <- current_sim_otus

          if(length(sim_data_labels) > 0){
            # First, run stabENG to get raw simulation network
            stabENG_args_sim_data <- c(
              list(Y = Sim_list[[i]][[j]], labels = sim_data_labels),
              stabENG_parameters
            )
            Res_sim_raw <- do.call(stabENG, stabENG_args_sim_data)
            rm(stabENG_args_sim_data);
            
            # Store raw results
            Res_sim[[i]][[j]] <- Res_sim_raw
            
            # Remove diagonal elements from raw results
            if(!is.null(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus)) {
                diag(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus) <- 0
            }
            if(!is.null(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus)) {
                diag(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus) <- 0
            }
            
            # Get optimal lambdas from simulation
            sim_opt_lambda1 <- Res_sim_raw$opt.lambda1
            sim_opt_lambda2 <- Res_sim_raw$opt.lambda2
            
            cat(paste('Performing Freedman-Lane permutation test for simulation', j, 'on day', i, '\n'))
            
            # Perform permutation test on simulation data for Nplus
            sim_perm_threshold_nplus <- if (!is.null(Sim_list[[i]][[j]]$Nplus) && nrow(Sim_list[[i]][[j]]$Nplus) > 0 && 
                                             !is.null(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus)) {
                permtest_freedman_lane(
                  data_list_current_group = Sim_list[[i]][[j]]$Nplus,
                  original_other_group_data = Sim_list[[i]][[j]]$Nminus,
                  pcor_matrix_target_group = Res_sim[[i]][[j]]$opt.fit.pcor$Nplus,
                  group_name_to_permute = "Nplus",
                  n_perm = 25, # Reduced permutations for simulation efficiency
                  alpha = 0.1,
                  stabENG_general_params = stabENG_parameters,
                  opt_lambda1_current = sim_opt_lambda1,
                  opt_lambda2_current = sim_opt_lambda2,
                  shared_labels_for_stabENG = sim_data_labels
                )
            } else { NA }
            
            # Perform permutation test on simulation data for Nminus
            sim_perm_threshold_nminus <- if (!is.null(Sim_list[[i]][[j]]$Nminus) && nrow(Sim_list[[i]][[j]]$Nminus) > 0 && 
                                              !is.null(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus)) {
                permtest_freedman_lane(
                  data_list_current_group = Sim_list[[i]][[j]]$Nminus,
                  original_other_group_data = Sim_list[[i]][[j]]$Nplus,
                  pcor_matrix_target_group = Res_sim[[i]][[j]]$opt.fit.pcor$Nminus,
                  group_name_to_permute = "Nminus",
                  n_perm = 25, # Reduced permutations for simulation efficiency
                  alpha = 0.1,
                  stabENG_general_params = stabENG_parameters,
                  opt_lambda1_current = sim_opt_lambda1,
                  opt_lambda2_current = sim_opt_lambda2,
                  shared_labels_for_stabENG = sim_data_labels
                )
            } else { NA }
            
            # Store simulation-specific thresholds
            sim_perm_thresholds_list[[j]] <- list(Nplus = sim_perm_threshold_nplus, Nminus = sim_perm_threshold_nminus)
            
            cat(sprintf("Simulation %d permutation test thresholds for Day %s:\n  Nplus: %.4f\n  Nminus: %.4f\n",
                       j, i, ifelse(is.na(sim_perm_threshold_nplus), 0, sim_perm_threshold_nplus), 
                       ifelse(is.na(sim_perm_threshold_nminus), 0, sim_perm_threshold_nminus)))
            
            # Apply simulation-specific permutation thresholds to simulation results
            if(!is.null(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus) && !is.na(sim_perm_threshold_nplus)) {
                Res_sim[[i]][[j]]$opt.fit.pcor$Nplus[abs(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus) < sim_perm_threshold_nplus] <- 0
                cat(paste("Applied simulation-specific Nplus permutation threshold", round(sim_perm_threshold_nplus, 4), "to simulation", j, "\n"))
            }
            if(!is.null(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus) && !is.na(sim_perm_threshold_nminus)) {
                Res_sim[[i]][[j]]$opt.fit.pcor$Nminus[abs(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus) < sim_perm_threshold_nminus] <- 0
                cat(paste("Applied simulation-specific Nminus permutation threshold", round(sim_perm_threshold_nminus, 4), "to simulation", j, "\n"))
            }
            
            rm(Res_sim_raw); # Clean up
            
          } else {
             warning(paste0("Skipping stabENG for simulation ", j, " on Day ", i, " due to no labels for sim data."))
             Res_sim[[i]][[j]] <- list(opt.fit.pcor=list(Nplus=matrix(0,0,0), Nminus=matrix(0,0,0)))
             sim_perm_thresholds_list[[j]] <- list(Nplus = NA, Nminus = NA)
          }
      } else { 
          warning(paste0("Skipping stabENG for simulation ", j, " on Day ", i, " due to empty simulation data."))
          Res_sim[[i]][[j]] <- list(opt.fit.pcor=list(Nplus=matrix(0,0,0), Nminus=matrix(0,0,0)))
          sim_perm_thresholds_list[[j]] <- list(Nplus = NA, Nminus = NA)
      }
      
      cat('Simulation Progress on Day',i,':', 100*j/n_sim,'% done.\n')
      
      # Count edges after applying simulation-specific permutation thresholds
      nplus_edges_after_sim_perm <- if(!is.null(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus)) {
          sum(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus[upper.tri(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus)] != 0)
      } else { 0 }
      nminus_edges_after_sim_perm <- if(!is.null(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus)) {
          sum(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus[upper.tri(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus)] != 0)
      } else { 0 }
      
      cat(paste("Simulation", j, "edges after simulation-specific permutation filtering - Nplus:", nplus_edges_after_sim_perm, 
                "Nminus:", nminus_edges_after_sim_perm, "\n"))
      
      # Plot simulation network (using simulation-specific permutation filtered results)
      if (length(sim_data_labels) > 0) {
          Phylum_groups_sim <- as.factor(otu_tax[sim_data_labels,"Phylum"])
          
          if(!is.null(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus) && nrow(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus) > 0){
            png(filename=paste0(plot_path,"_simulation_network_",j,"_Nplus_Phylum_SimPermFiltered_vsized.png"))
            qgraph::qgraph(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus, 
              layout = "circle",
              edge.color = ifelse(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus > 0, "#2b732b", "red"),
              title = paste0("Sim ", j, " Network Nplus by Phylum (Day ", i, ") - Sim-Specific PermTest Filtered"),
              vsize = 2.5,
              groups = Phylum_groups_sim)
            dev.off()
          }
          
          if(!is.null(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus) && nrow(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus) > 0){
            png(filename=paste0(plot_path,"_simulation_network_",j,"_Nminus_Phylum_SimPermFiltered_vsized.png"))
            qgraph::qgraph(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus, 
              layout = "circle",
              edge.color = ifelse(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus > 0, "#2b732b", "red"),
              title = paste0("Sim ", j, " Network Nminus by Phylum (Day ", i, ") - Sim-Specific PermTest Filtered"),
              vsize = 2.5,
              groups = Phylum_groups_sim)
            dev.off()
          }
      }
  }
  gc()
  
  # True adjacency for confusion matrix should be based on the perm-filtered real data network
  true_adj_Nplus_for_confusion <- if (!is.null(network_pcor[[i]]$Nplus) && nrow(network_pcor[[i]]$Nplus) > 0 && ncol(network_pcor[[i]]$Nplus) > 0) {
                                      (network_pcor[[i]]$Nplus !=0)*1 # Already filtered to current_sim_otus effectively by shared_otu_perm_filtered logic
                                  } else { matrix(0, length(current_sim_otus), length(current_sim_otus), dimnames=list(current_sim_otus, current_sim_otus)) }
  true_adj_Nminus_for_confusion <- if (!is.null(network_pcor[[i]]$Nminus) && nrow(network_pcor[[i]]$Nminus) > 0 && ncol(network_pcor[[i]]$Nminus) > 0) {
                                       (network_pcor[[i]]$Nminus !=0)*1
                                   } else { matrix(0, length(current_sim_otus), length(current_sim_otus), dimnames=list(current_sim_otus, current_sim_otus)) }
  
  # Ensure row/col names for true_adj matrices if they became 0x0
   if(length(current_sim_otus) > 0){
      if(nrow(true_adj_Nplus_for_confusion) == 0 && ncol(true_adj_Nplus_for_confusion) == 0 && length(current_sim_otus) > 0) {
          dimnames(true_adj_Nplus_for_confusion) <- list(current_sim_otus, current_sim_otus)
      }
      if(nrow(true_adj_Nminus_for_confusion) == 0 && ncol(true_adj_Nminus_for_confusion) == 0 && length(current_sim_otus) > 0) {
          dimnames(true_adj_Nminus_for_confusion) <- list(current_sim_otus, current_sim_otus)
      }
  }

  confusion_results <- lapply(1:n_sim, function(j_sim) { # Renamed loop variable
      sim_pcor_nplus <- Res_sim[[i]][[j_sim]]$opt.fit.pcor$Nplus
      sim_pcor_nminus <- Res_sim[[i]][[j_sim]]$opt.fit.pcor$Nminus

      # Adjacency from simulation (ensure it's over current_sim_otus)
      adj_sim_nplus <- if (!is.null(sim_pcor_nplus) && nrow(sim_pcor_nplus) > 0 && ncol(sim_pcor_nplus) > 0) {
                            current_otus_in_sim_nplus <- colnames(sim_pcor_nplus)
                            shared_otus_for_comp_nplus <- intersect(current_sim_otus, current_otus_in_sim_nplus)
                            if(length(shared_otus_for_comp_nplus) > 0) {
                                (sim_pcor_nplus[shared_otus_for_comp_nplus, shared_otus_for_comp_nplus, drop=FALSE] !=0)*1
                            } else { matrix(0, length(current_sim_otus), length(current_sim_otus), dimnames=list(current_sim_otus, current_sim_otus))}
                       } else { matrix(0, length(current_sim_otus), length(current_sim_otus), dimnames=list(current_sim_otus, current_sim_otus)) }
      
      adj_sim_nminus <- if (!is.null(sim_pcor_nminus) && nrow(sim_pcor_nminus) > 0 && ncol(sim_pcor_nminus) > 0) {
                            current_otus_in_sim_nminus <- colnames(sim_pcor_nminus)
                            shared_otus_for_comp_nminus <- intersect(current_sim_otus, current_otus_in_sim_nminus)
                            if(length(shared_otus_for_comp_nminus) > 0) {
                               (sim_pcor_nminus[shared_otus_for_comp_nminus, shared_otus_for_comp_nminus, drop=FALSE] !=0)*1
                            } else { matrix(0, length(current_sim_otus), length(current_sim_otus), dimnames=list(current_sim_otus, current_sim_otus))}
                        } else { matrix(0, length(current_sim_otus), length(current_sim_otus), dimnames=list(current_sim_otus, current_sim_otus)) }

      # Align true_adj to the actual shared OTUs used for this specific simulation comparison
      true_adj_nplus_aligned <- if(ncol(adj_sim_nplus) > 0) true_adj_Nplus_for_confusion[colnames(adj_sim_nplus), colnames(adj_sim_nplus), drop=FALSE] else matrix(0,0,0)
      true_adj_nminus_aligned <- if(ncol(adj_sim_nminus) > 0) true_adj_Nminus_for_confusion[colnames(adj_sim_nminus), colnames(adj_sim_nminus), drop=FALSE] else matrix(0,0,0)


      Nplus_metrics <- if(length(current_sim_otus) > 0 && ncol(true_adj_nplus_aligned) > 0 && ncol(adj_sim_nplus) > 0 && 
                          all(dim(true_adj_nplus_aligned) == dim(adj_sim_nplus))) {
                           calculate_metrics(true_adj_nplus_aligned, adj_sim_nplus)
                       } else { list(TPR=NA, FPR=NA, Precision=NA, Recall=NA, F1=NA, AUC=NA, MCC=NA, TP=NA,FP=NA,TN=NA,FN=NA) } # Added TP/FP/TN/FN to NA list
      Nminus_metrics <- if(length(current_sim_otus) > 0 && ncol(true_adj_nminus_aligned) > 0 && ncol(adj_sim_nminus) > 0 &&
                           all(dim(true_adj_nminus_aligned) == dim(adj_sim_nminus))) {
                            calculate_metrics(true_adj_nminus_aligned, adj_sim_nminus)
                        } else { list(TPR=NA, FPR=NA, Precision=NA, Recall=NA, F1=NA, AUC=NA, MCC=NA, TP=NA,FP=NA,TN=NA,FN=NA) } # Added TP/FP/TN/FN to NA list
      return(list(Nplus = Nplus_metrics, Nminus = Nminus_metrics))
  })
  results_df <- do.call(rbind, lapply(confusion_results, function(x) {
    # Ensure that Nplus and Nminus components are data.frame compatible
    nplus_df <- if(is.list(x$Nplus)) as.data.frame(x$Nplus) else data.frame(TPR=NA, FPR=NA, Precision=NA, Recall=NA, F1=NA, AUC=NA, MCC=NA)
    nminus_df <- if(is.list(x$Nminus)) as.data.frame(x$Nminus) else data.frame(TPR=NA, FPR=NA, Precision=NA, Recall=NA, F1=NA, AUC=NA, MCC=NA)
    # Add prefixes
    names(nplus_df) <- paste0("Nplus.", names(nplus_df))
    names(nminus_df) <- paste0("Nminus.", names(nminus_df))
    cbind(nplus_df, nminus_df)
  }))
  rm(confusion_results); gc()

  if(nrow(results_df) > 0){
    results_df_long <- results_df %>%
      # dplyr::select(starts_with("Nplus.") | starts_with("Nminus.")) %>% # Already selected by cbind
      tidyr::pivot_longer(cols = everything(),
                  names_to = c("group", "metric"),
                  names_sep = "\\.",
                  values_to = "value")
    # Add matrix_id and times, ensuring correct lengths    num_metrics_per_sim_iteration <- ncol(results_df) # Nplus.TPR, Nplus.FPR etc. + Nminus.TPR etc.
    if (nrow(results_df_long) > 0 && num_metrics_per_sim_iteration > 0) {
        results_df_long <- results_df_long %>%
            dplyr::mutate(matrix_id = rep(1:n_sim, each = num_metrics_per_sim_iteration),
                          times = rep(i, nrow(results_df_long)))
    } else if (nrow(results_df_long) == 0 && n_sim > 0) { # Handle case where results_df_long is empty but shouldn't be
        warning(paste0("results_df_long is empty for Day ", i, " but n_sim > 0. Metrics might be all NA."))
        # Create an empty df with correct columns if needed for rbind
        results_df_long <- data.frame(group=character(), metric=character(), value=numeric(), matrix_id=integer(), times=character())
    }
    rm(results_df); gc()
    confusion_results_df <- rbind(confusion_results_df, results_df_long)
  } else {
    warning(paste0("No simulation metrics generated for Day ", i))
  }
  gc()
}

# Filter for TPR, FPR, and F1 metrics
if (is.data.frame(confusion_results_df) && nrow(confusion_results_df) > 0) {
    selected_metrics <- c("TPR", "FPR", "F1", "Precision", "Recall", "AUC", "MCC") # Expanded for better overview
    filtered_data <- confusion_results_df %>% 
        filter(metric %in% selected_metrics) %>%
        mutate(metric = factor(metric, levels = selected_metrics)) # Ensure order in plot

    # Plot
    p <- ggplot(filtered_data, aes(x = factor(as.integer(times)), y = value, fill = group)) +
      geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "black", position = position_dodge(0.8), na.rm = TRUE) +
      scale_fill_manual(values = c("Nplus" = "#00BFC4", "Nminus" = "#F8766D")) +
      facet_wrap(~metric, scales = "free_y", ncol = 4) + 
      labs(x = "Days", y = "Value", title = "Confusion Metrics (Permutation Test Filtered Networks)") +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "top",
        panel.grid.major = element_line(color = "gray90"),
        panel.border = element_rect(color = "black", fill = NA),
        strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels if many timepoints
      )
    ggsave(filename = file.path(plot_dir, "PermTest_Filtered_Confusion_boxplot_AllDays.png"), p, width=12, height=8) # Save in the last day's plot_dir or a general one
} else {
    cat("confusion_results_df is empty. Skipping final boxplot.\n")
}

saveRDS(network_pcor, file = file.path("DataImage", "network_pcor_list_PermTest.rds"))
saveRDS(network_list, file = file.path("DataImage", "network_precision_list_PermTest.rds"))
save.image(file.path("DataImage", paste0("big1226_Days_network_results_PermTest_Filtered_",format(Sys.Date(), "%Y%m%d"),".RData")))
cat("All done in: ", Sys.time() - start_time, "\n")