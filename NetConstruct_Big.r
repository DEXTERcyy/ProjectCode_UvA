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
cat("Packages loaded successfully\n")
#load("//home//14720078//ProjectCode//DataImage//big1226_Days_network_results_Big_Days_Filtered.RData")
start_time <- Sys.time() # Moved start_time here, after packages are loaded.
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
confusion_results_df <- list()
Sim_list <- list()
Res_sim <- list()
lambda1 <- list()
lambda2 <- list()
# pcor_screen <- function(
#   data_list, 
#   shared_otu, 
#   n_perm = 50, 
#   lambda1, 
#   lambda2, 
#   alpha = 0.1
# )
#   {
#     n_features <- length(shared_otu)
#     threshold <- list()
#     perm_data <- list()
#     # sorted pcor
#     null_pcors_Nplus <- numeric()
#     null_pcors_Nminus <- numeric()
    
#     # progress bar
#     pb <- txtProgressBar(min = 0, max = n_perm, style = 3)
    
#     for (p in 1:n_perm) {
#       perm_data$Nplus <- data_list$Nplus
#       perm_data$Nminus <- data_list$Nminus
#       for (j in 1:n_features) {
#         perm_data$Nplus[, j] <- sample(data_list$Nplus[, j])
#         perm_data$Nminus[, j] <- sample(data_list$Nminus[, j])
#       }
#       pcor_perm <- preprocess_and_estimate_network(perm_data, shared_otu, lambda1, lambda2)$pcor
#       diag(pcor_perm$Nplus) <- 0
#       diag(pcor_perm$Nminus) <- 0
#       null_pcors_Nplus <- c(null_pcors_Nplus, abs(pcor_perm$Nplus[upper.tri(pcor_perm$Nplus)]))
#       null_pcors_Nminus <- c(null_pcors_Nminus, abs(pcor_perm$Nminus[upper.tri(pcor_perm$Nminus)]))
#       setTxtProgressBar(pb, p)
#     }
    
#     close(pb)
    
#     threshold$Nplus <- quantile(null_pcors_Nplus, 1 - alpha)
#     threshold$Nminus <- quantile(null_pcors_Nminus, 1 - alpha)
    
#     return(threshold)
#   }
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
  rep.num = 15,
  nlambda1 = 20,
  lambda1.min = 0.01,
  lambda1.max = 1,
  nlambda2 = 20,
  lambda2.min = 0,
  lambda2.max = 0.1,
  lambda2.init = 0.01,
  ebic.gamma = 0.5,
  nCores = if (Sys.getenv("SLURM_CPUS_PER_TASK") != "") as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
)
# %% calculate network
timestamps <- timestamps[4]
for (i in timestamps){
  # Create the directory structure if it doesn't exist
  plot_dir <- file.path("Plots", "BigDataFilter")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  plot_path <- file.path(plot_dir, paste0("Day_", i))
  data_list_current_timestamp <- data_list_times[[i]] # Use a more descriptive name
  cat('Calculating network on day ',i,'\n')
  
  stabENG_args_real_data <- c(
    list(Y = data_list_current_timestamp, labels = shared_otu),
    stabENG_parameters
  )
  network_results <- do.call(stabENG, stabENG_args_real_data)
  cat('Network calculation on Day', i, 'completed.\n')
  
  network_list_raw[[i]]$Nplus <- network_results$opt.fit$Nplus # precision matrix estimates
  network_list_raw[[i]]$Nminus <- network_results$opt.fit$Nminus # precision matrix estimates
  diag(network_list_raw[[i]]$Nplus) = diag(network_list_raw[[i]]$Nminus) <- 0
  network_pcor_raw[[i]]$Nplus <- network_results$opt.fit.pcor$Nplus
  network_pcor_raw[[i]]$Nminus <- network_results$opt.fit.pcor$Nminus
  diag(network_pcor_raw[[i]]$Nplus) = diag(network_pcor_raw[[i]]$Nminus) <- 0
  lambda1[[i]] <- network_results$opt.lambda1
  lambda2[[i]] <- network_results$opt.lambda2
  
  rm(network_results, stabENG_args_real_data); gc() # Remove large stabENG output and args

  if (PCS_filter == TRUE) {
    cat('\nPerforming PCS for network on Day',i,'\n')
    pcs_thresholds_values <- list() # Initialize list for thresholds for current timestamp
    # PCS for Nplus group
    if (!is.null(data_list_current_timestamp$Nplus) && nrow(data_list_current_timestamp$Nplus) > 0) {
        pcs_thresholds_values$Nplus <- pcs_cv_threshold_stabENG_lasso_perm(
          data_list_unscaled_full = data_list_current_timestamp,
          initial_pcor_matrix_target = network_pcor_raw[[i]]$Nplus,
          group_name_target = "Nplus",
          other_group_name = "Nminus", # 指定另一个组的名称
          stabENG_params_list = stabENG_parameters, 
          otu_labels = shared_otu,
          fold = nrow(data_list_current_timestamp$Nplus), # 使用LOOCV
          plot_cv_curve = TRUE, # Changed to TRUE for Nplus as well for testing
          plot_path_prefix = paste0(plot_path, "_pcs_Nplus")
        )
    } else {
        warning(paste0("Skipping PCS for Nplus on Day ", i, " due to no data for Nplus group."))
        pcs_thresholds_values$Nplus <- 0.05 # Default tau if no data
    }

    # PCS for Nminus group
    if (!is.null(data_list_current_timestamp$Nminus) && nrow(data_list_current_timestamp$Nminus) > 0) {
        pcs_thresholds_values$Nminus <- pcs_cv_threshold_stabENG_lasso_perm(
          data_list_unscaled_full = data_list_current_timestamp, # 传递包含Nplus和Nminus的完整列表
          initial_pcor_matrix_target = network_pcor_raw[[i]]$Nminus,
          group_name_target = "Nminus",
          other_group_name = "Nplus",
          stabENG_params_list = stabENG_parameters, 
          otu_labels = shared_otu,
          fold = nrow(data_list_current_timestamp$Nminus),
          plot_cv_curve = TRUE, 
          plot_path_prefix = paste0(plot_path, "_pcs_Nminus")
        )
    } else {
        warning(paste0("Skipping PCS for Nminus on Day ", i, " due to no data for Nminus group."))
        pcs_thresholds_values$Nminus <- 0.05 # Default tau if no data
    }

    # filter edge sparsity for simulation
    network_list[[i]]$Nplus <- network_list_raw[[i]]$Nplus
    network_list[[i]]$Nminus <- network_list_raw[[i]]$Nminus
    network_list[[i]]$Nplus[abs(network_list[[i]]$Nplus) < pcs_thresholds_values$Nplus] <- 0
    network_list[[i]]$Nminus[abs(network_list[[i]]$Nminus) < pcs_thresholds_values$Nminus] <- 0
    cat('Filtered OTU number on Day', i, 'Nplus (precision):', nrow(network_list[[i]]$Nplus), 
        'Nminus (precision):', nrow(network_list[[i]]$Nminus), '\n')
    # check if Nplus and Nminus matrices have rownames and colnames
    if (is.null(rownames(network_list[[i]]$Nplus)) || is.null(colnames(network_list[[i]]$Nplus))) {
        warning(paste0("Nplus matrix on Day ", i, " has no rownames or colnames."))
        network_list[[i]]$Nplus <- matrix(0,0,0, dimnames=list(NULL,NULL))
    }
    # Refine network_list based on its own non-zero rows/cols first
    if (!is.null(network_list[[i]]$Nplus) && nrow(network_list[[i]]$Nplus) > 0 && ncol(network_list[[i]]$Nplus) > 0) {
        rows_to_keep_nplus <- apply(network_list[[i]]$Nplus, 1, function(x) any(x != 0))
        cols_to_keep_nplus <- apply(network_list[[i]]$Nplus, 2, function(x) any(x != 0))
        common_otus_self_nplus <- rownames(network_list[[i]]$Nplus)[rows_to_keep_nplus & cols_to_keep_nplus]
        if (length(common_otus_self_nplus) > 0) {
            network_list[[i]]$Nplus <- network_list[[i]]$Nplus[common_otus_self_nplus, common_otus_self_nplus, drop = FALSE]
        } else {
            network_list[[i]]$Nplus <- matrix(0,0,0, dimnames=list(NULL,NULL))
            warning(paste0("No common OTUs found in Nplus after PCS filtering on Day ", i, ". Setting Nplus to empty matrix."))
        }
    } else {
        network_list[[i]]$Nplus <- matrix(0,0,0, dimnames=list(NULL,NULL))
        warning(paste0("Nplus matrix is empty or invalid on Day ", i, ". Setting Nplus to empty matrix."))
    }

    if (!is.null(network_list[[i]]$Nminus) && nrow(network_list[[i]]$Nminus) > 0 && ncol(network_list[[i]]$Nminus) > 0) {
        rows_to_keep_nminus <- apply(network_list[[i]]$Nminus, 1, function(x) any(x != 0))
        cols_to_keep_nminus <- apply(network_list[[i]]$Nminus, 2, function(x) any(x != 0))
        common_otus_self_nminus <- rownames(network_list[[i]]$Nminus)[rows_to_keep_nminus & cols_to_keep_nminus]
        if (length(common_otus_self_nminus) > 0) {
            network_list[[i]]$Nminus <- network_list[[i]]$Nminus[common_otus_self_nminus, common_otus_self_nminus, drop = FALSE]
        } else {
            network_list[[i]]$Nminus <- matrix(0,0,0, dimnames=list(NULL,NULL))
            warning(paste0("No common OTUs found in Nminus after PCS filtering on Day ", i, ". Setting Nminus to empty matrix."))
        }
    } else {
        network_list[[i]]$Nminus <- matrix(0,0,0, dimnames=list(NULL,NULL))
        warning(paste0("Nminus matrix is empty or invalid on Day ", i, ". Setting Nminus to empty matrix."))
    }
    cat('Shared filtered OTU number on Day', i, 'Nplus (precision):', nrow(network_list[[i]]$Nplus), 
        'Nminus (precision):', nrow(network_list[[i]]$Nminus), '\n')

    # Determine shared OTUs after individual PCS filtering of precision matrices
    # Ensure that rownames exist and are valid before intersect
    otus_in_nplus_prec <- if(!is.null(rownames(network_list[[i]]$Nplus))) rownames(network_list[[i]]$Nplus) else character(0)
    otus_in_nminus_prec <- if(!is.null(rownames(network_list[[i]]$Nminus))) rownames(network_list[[i]]$Nminus) else character(0)
    
    if (length(otus_in_nplus_prec) == 0 && length(otus_in_nminus_prec) == 0) {
        shared_otu_PCS <- character(0) # No OTUs in either
    } else if (length(otus_in_nplus_prec) == 0) {
        shared_otu_PCS <- otus_in_nminus_prec # Use OTUs from Nminus if Nplus is empty
    } else if (length(otus_in_nminus_prec) == 0) {
        shared_otu_PCS <- otus_in_nplus_prec # Use OTUs from Nplus if Nminus is empty
    } else {
        shared_otu_PCS <- intersect(otus_in_nplus_prec, otus_in_nminus_prec)
    }
    cat('Shared OTU number for simulation basis (shared_otu_PCS) on Day', i, ':', length(shared_otu_PCS), '\n')

    # Further ensure network_list matrices are subsetted to these shared_otu_PCS for simulation "true" networks
    if (length(shared_otu_PCS) > 0) {
        if (nrow(network_list[[i]]$Nplus) > 0) { # Check if Nplus matrix is not already empty
             network_list[[i]]$Nplus <- network_list[[i]]$Nplus[rownames(network_list[[i]]$Nplus) %in% shared_otu_PCS, 
                                                                colnames(network_list[[i]]$Nplus) %in% shared_otu_PCS, drop = FALSE]
        }
        if (nrow(network_list[[i]]$Nminus) > 0) { # Check if Nminus matrix is not already empty
             network_list[[i]]$Nminus <- network_list[[i]]$Nminus[rownames(network_list[[i]]$Nminus) %in% shared_otu_PCS, 
                                                                  colnames(network_list[[i]]$Nminus) %in% shared_otu_PCS, drop = FALSE]
        }
    } else {
        # If shared_otu_PCS is empty, simulation might not be meaningful or possible with current setup
        warning(paste0("No shared OTUs after PCS filtering for Day ", i, ". Simulation might be problematic."))
        network_list[[i]]$Nplus <- matrix(0,0,0, dimnames=list(NULL,NULL))
        network_list[[i]]$Nminus <- matrix(0,0,0, dimnames=list(NULL,NULL))
    }
    
    # The network_pcor matrices should also be filtered consistently for plotting and "true" network in confusion matrix
    # Initialize network_pcor[[i]] if it's the first assignment for this i
    if (is.null(network_pcor[[i]])) network_pcor[[i]] <- list()

    if (!is.null(network_pcor_raw[[i]]$Nplus)) {
        network_pcor[[i]]$Nplus <- network_pcor_raw[[i]]$Nplus
        network_pcor[[i]]$Nplus[abs(network_pcor[[i]]$Nplus) < pcs_thresholds_values$Nplus] <- 0
        if (length(shared_otu_PCS) > 0 && nrow(network_pcor[[i]]$Nplus) > 0) {
            network_pcor[[i]]$Nplus <- network_pcor[[i]]$Nplus[rownames(network_pcor[[i]]$Nplus) %in% shared_otu_PCS, 
                                                               colnames(network_pcor[[i]]$Nplus) %in% shared_otu_PCS, drop = FALSE]
        } else if (length(shared_otu_PCS) == 0 && nrow(network_pcor[[i]]$Nplus) > 0) { # If shared_otu_PCS is empty, make pcor empty too
             network_pcor[[i]]$Nplus <- matrix(0,0,0, dimnames=list(NULL,NULL))
        } else if (nrow(network_pcor[[i]]$Nplus) == 0 && ncol(network_pcor[[i]]$Nplus) == 0) {
            # It's already an empty matrix, do nothing
        } else { # Default to empty if conditions are tricky
            network_pcor[[i]]$Nplus <- matrix(0,0,0, dimnames=list(NULL,NULL))
        }
    } else {
        network_pcor[[i]]$Nplus <- matrix(0,0,0, dimnames=list(NULL,NULL))
    }

    if (!is.null(network_pcor_raw[[i]]$Nminus)) {
        network_pcor[[i]]$Nminus <- network_pcor_raw[[i]]$Nminus
        network_pcor[[i]]$Nminus[abs(network_pcor[[i]]$Nminus) < pcs_thresholds_values$Nminus] <- 0
        if (length(shared_otu_PCS) > 0 && nrow(network_pcor[[i]]$Nminus) > 0) {
            network_pcor[[i]]$Nminus <- network_pcor[[i]]$Nminus[rownames(network_pcor[[i]]$Nminus) %in% shared_otu_PCS, 
                                                                 colnames(network_pcor[[i]]$Nminus) %in% shared_otu_PCS, drop = FALSE]
        } else if (length(shared_otu_PCS) == 0 && nrow(network_pcor[[i]]$Nminus) > 0) {
            network_pcor[[i]]$Nminus <- matrix(0,0,0, dimnames=list(NULL,NULL))
        } else if (nrow(network_pcor[[i]]$Nminus) == 0 && ncol(network_pcor[[i]]$Nminus) == 0) {
            # It's already an empty matrix, do nothing
        } else {
            network_pcor[[i]]$Nminus <- matrix(0,0,0, dimnames=list(NULL,NULL))
        }
    } else {
        network_pcor[[i]]$Nminus <- matrix(0,0,0, dimnames=list(NULL,NULL))
    }
    cat('Final OTU number for Nplus pcor on Day', i, ':', if(!is.null(network_pcor[[i]]$Nplus)) nrow(network_pcor[[i]]$Nplus) else 0, '\n')
    cat('Final OTU number for Nminus pcor on Day', i, ':', if(!is.null(network_pcor[[i]]$Nminus)) nrow(network_pcor[[i]]$Nminus) else 0, '\n')

  } else {
    cat('\nSkipping permutation test filtering on Day',i,'\n')
    # Use raw pcor
    network_list[[i]] <- network_list_raw[[i]]
    network_pcor[[i]] <- network_pcor_raw[[i]]
    PCS_threshold_Nplus <- 0
    PCS_threshold_Nminus <- 0
  }
  # %% Plot network on Phylum level
  # Phylum_groups <- as.factor(otu_tax[rownames(network_pcor[[i]]$Nplus),"Phylum"])
  # png(filename=paste0(plot_path,"_network_Nplus_Phylum_Stab_Filtered_vsized.png"))
  # qgraph::qgraph(network_pcor[[i]]$Nplus, 
  #   layout = "circle",
  #   edge.color = ifelse(network_pcor[[i]]$Nplus > 0, "blue", "red"),
  #   title = "Stab Network Nplus by Phylum",
  #   vsize = 2.5,
  #   groups = Phylum_groups)
  # dev.off()
  
  # png(filename=paste0(plot_path,"_network_Nminus_Phylum_Stab_Filtered_vsized.png"))
  # qgraph::qgraph(network_pcor[[i]]$Nminus, 
  #   layout = "circle",
  #   edge.color = ifelse(network_pcor[[i]]$Nminus > 0, "blue", "red"),
  #   title = "Stab Network Nminus by Phylum",
  #   vsize = 2.5,
  #   groups = Phylum_groups)
  # dev.off()
  # # %%Visualize Edge weights
  # cor_values_Nplus <- as.vector(network_pcor[[i]]$Nplus)
  # cor_values_Nminus <- as.vector(network_pcor[[i]]$Nminus)
  # cor_values_Nplus_nonzero <- cor_values_Nplus[cor_values_Nplus != 0]
  # cor_values_Nminus_nonzero <- cor_values_Nminus[cor_values_Nminus != 0]

  # ggplot(data.frame(Type = c(rep("N+", length(cor_values_Nplus_nonzero)), rep("N-", length(cor_values_Nminus_nonzero))), 
  #   Value = c(cor_values_Nplus_nonzero, cor_values_Nminus_nonzero)), 
  #   aes(x = Value, color = Type)) + 
  #   geom_density() + 
  #   labs(x = "Correlation Value", y = "Density")
  # ggsave(filename=paste0(plot_path,"_correlation_distribution.png"))
  gc()
}
