# #%% requirements
# # !ebic.gamma: if too sparse change to 0.6
cat(R.version.string, "\nR Home:", R.home(), "\nR lib:", .libPaths(), "\n")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# library(devtools, quietly = TRUE)
# library(BiocManager, quietly = TRUE)
# if (!requireNamespace("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
# if (!requireNamespace("SPRING", quietly = TRUE)) devtools::install_github("GraceYoon/SPRING")
# if (!requireNamespace("SpiecEasi", quietly = TRUE)) devtools::install_github("zdk123/SpiecEasi")
# if (!requireNamespace("stabJGL", quietly = TRUE)) devtools::install_github("camiling/stabJGL")
# if (!requireNamespace("NetCoMi", quietly = TRUE)) devtools::install_github("stefpeschel/NetCoMi", ref = "develop",repos = c("https://cloud.r-project.org/",BiocManager::repositories()))
# if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
# if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
# if (!requireNamespace("MLmetrics", quietly = TRUE)) install.packages("MLmetrics")
# if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
# if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
# if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr") # Likely installed with tidyverse
# if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")
# if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse") # Installs a suite of packages
# if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
# if (!requireNamespace("foreach", quietly = TRUE)) install.packages("foreach")
# if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
# if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel") # Usually base R, no need to install
# if (!requireNamespace("mixedCCA", quietly = TRUE)) install.packages("mixedCCA")
# if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
# if (!requireNamespace("qgraph", quietly = TRUE)) install.packages("qgraph")
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
source("//home//14720078//ProjectCode//Packages//stabENG.r")
source("//home//14720078//ProjectCode//Packages//MyPlot.r")
cat("Packages loaded successfully\n")
#load("//home//14720078//ProjectCode//DataImage//big1226_Days_network_results_Big_Days_Filtered.RData")
start_time <- Sys.time()
rawdata <- readRDS("//home//14720078//ProjectCode//data//DavarData1_substrate_phyloseq_1226_final_filtered.rds")
otu_raw <- otu_table(rawdata)
otu_RA <- transform_sample_counts(otu_raw, function(x) x / sum(x) )
otu_RA <- filter_taxa(otu_RA, function(x) mean(x) > 1e-3, TRUE)
shared_otu <- rownames(otu_RA)
otu_Ab <- as.data.frame(t(otu_raw[shared_otu,]))
sam_info <- as.data.frame(sample_data(rawdata))
otu_tax <- as.data.frame(tax_table(rawdata))
taxa_filtering <- TRUE
n_sim <- 15
perm_filter <- FALSE
sim_perm_filter <- FALSE
# %%split otu_Ab by condition group
# otu_Ab_Nplus <- otu_Ab[rownames(otu_Ab) %in% rownames(sam_info[sam_info$growthCondition=="plusN",]),]
# otu_Ab_Nminus <- otu_Ab[rownames(otu_Ab) %in% rownames(sam_info[sam_info$growthCondition=="minusN",]),]
# data_list <- list(Nplus = otu_Ab_Nplus, Nminus = otu_Ab_Nminus)
# by Group and Days
timestamps <- as.character(sort(as.integer(levels(sam_info$Days))))
#timestamps <- timestamps[1]
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
  if (taxa_filtering == TRUE) {
    # Sum the number of zeros in each column and sort by the number of zeros
    zero_counts_Nplus <- colSums(otu_Ab_Nplus_times[[i]] == 0)
    zero_counts_Nminus <- colSums(otu_Ab_Nminus_times[[i]] == 0)
    
    # Select the top 100 columns with the least zeros for each group
    sorted_Nplus <- names(sort(zero_counts_Nplus, decreasing = FALSE))[1:100]
    sorted_Nminus <- names(sort(zero_counts_Nminus, decreasing = FALSE))[1:100]
    
    # Use the intersecting column names for Nplus & Nminus
    common_taxa[[i]] <- intersect(sorted_Nplus, sorted_Nminus)
  }
  data_list_times[[i]] <- list(Nplus = otu_Ab_Nplus_times[[i]], 
                               Nminus = otu_Ab_Nminus_times[[i]])
}
# common_taxa <- Filter(Negate(is.null), common_taxa)
# common_taxa <- Filter(function(x) length(x) > 0, common_taxa)
if (taxa_filtering == TRUE) {
  common_taxa_intersection <- Reduce(intersect, common_taxa)
  # Filter the data_list_times to only include the common taxa
  for (i in timestamps) {
    if (length(common_taxa_intersection) > 0) {
      data_list_times[[i]]$Nplus <- data_list_times[[i]]$Nplus[, common_taxa_intersection, drop = FALSE]
      data_list_times[[i]]$Nminus <- data_list_times[[i]]$Nminus[, common_taxa_intersection, drop = FALSE]
    } else {
      stop("No intersection found in common_taxa.")
    }
  }
  shared_otu <- common_taxa_intersection
  # Filter the otu_tax to only include the common taxa
  otu_tax <- otu_tax[common_taxa_intersection,]
  if (!all(common_taxa_intersection %in% rownames(otu_tax))) {
    warning("Some taxa in common_taxa_intersection are not in otu_tax.")
  }
  if (!all(common_taxa_intersection %in% shared_otu)) {
    warning("Some taxa in common_taxa_intersection are not in shared_otu.")
  }
}
# %%
network_list_raw <- list()
network_list_filtered <- list()
network_list <- list()
network_pcor_raw <- list()
neytwork_pcor_filtered <- list()
network_pcor <- list()
confusion_results_df <- list()
Sim_list <- list()
Res_sim <- list()
lambda1 <- list()
lambda2 <- list()
permutation_threshold <- function(
  data_list, 
  shared_otu, 
  n_perm = 50, 
  lambda1, 
  lambda2, 
  alpha, 
  plot_threshold = FALSE, 
  plot_path = NULL
)
  {
    n_features <- length(shared_otu)
    threshold <- list()
    perm_data <- list()
    # sorted pcor
    null_pcors_Nplus <- numeric()
    null_pcors_Nminus <- numeric()
    
    # progress bar
    pb <- txtProgressBar(min = 0, max = n_perm, style = 3)
    
    for (p in 1:n_perm) {
      perm_data$Nplus <- data_list$Nplus
      perm_data$Nminus <- data_list$Nminus
      for (j in 1:n_features) {
        perm_data$Nplus[, j] <- sample(data_list$Nplus[, j])
        perm_data$Nminus[, j] <- sample(data_list$Nminus[, j])
      }
      pcor_perm <- preprocess_and_estimate_network(perm_data, shared_otu, lambda1, lambda2)$pcor
      diag(pcor_perm$Nplus) <- 0
      diag(pcor_perm$Nminus) <- 0
      null_pcors_Nplus <- c(null_pcors_Nplus, abs(pcor_perm$Nplus[upper.tri(pcor_perm$Nplus)]))
      null_pcors_Nminus <- c(null_pcors_Nminus, abs(pcor_perm$Nminus[upper.tri(pcor_perm$Nminus)]))
      setTxtProgressBar(pb, p)
    }
    
    close(pb)
    
    threshold$Nplus <- quantile(null_pcors_Nplus, 1 - alpha)
    threshold$Nminus <- quantile(null_pcors_Nminus, 1 - alpha)
    
    # visualization
    if (plot_threshold) {
      # Nplus plot
      df_perm_Nplus <- data.frame(
      value = c(abs(pcor_perm$Nplus[upper.tri(pcor_perm$Nplus)]), null_pcors_Nplus),
      type = c(rep("Observed", sum(upper.tri(pcor_perm$Nplus))), 
           rep("Permuted", length(null_pcors_Nplus)))
      )
      png(filename = paste0(plot_path, "_permutation_threshold_Nplus.png"))
      p_nplus <- ggplot(df_perm_Nplus, aes(x = .data$value, fill = .data$type)) +
      geom_density(alpha = 0.5) +
      geom_vline(xintercept = threshold$Nplus, linetype = "dashed", color = "red") +
      labs(title = "Permutation Test for Nplus Network Threshold",
         subtitle = paste0("Threshold = ", round(threshold$Nplus, 4),
                   " (Significance level: ", alpha * 100, "%)"),
         x = "Absolute value of partial correlation",
         y = "Density") +
      theme_minimal() +
      scale_fill_manual(values = c("Observed" = "blue", "Permuted" = "gray50"))
      print(p_nplus) # Ensure the plot is printed to the file
      dev.off()

      # Nminus plot
      df_perm_Nminus <- data.frame(
      value = c(abs(pcor_perm$Nminus[upper.tri(pcor_perm$Nminus)]), null_pcors_Nminus),
      type = c(rep("Observed", sum(upper.tri(pcor_perm$Nminus))), 
           rep("Permuted", length(null_pcors_Nminus)))
      )
      png(filename = paste0(plot_path, "_permutation_threshold_Nminus.png"))
      p_nminus <- ggplot(df_perm_Nminus, aes(x = .data$value, fill = .data$type)) +
      geom_density(alpha = 0.5) +
      geom_vline(xintercept = threshold$Nminus, linetype = "dashed", color = "red") +
      labs(title = "Permutation Test for Nminus Network Threshold",
         subtitle = paste0("Threshold = ", round(threshold$Nminus, 4),
                   " (Significance level: ", alpha * 100, "%)"),
         x = "Absolute value of partial correlation",
         y = "Density") +
      theme_minimal() +
      scale_fill_manual(values = c("Observed" = "blue", "Permuted" = "gray50"))
      print(p_nminus) # Ensure the plot is printed to the file
      dev.off()
    }
    
    cat(sprintf("\nThresholds based on %d permutations and %.1f%% significance level:\n Nplus: %.4f\n Nminus: %.4f\n", 
          n_perm, alpha*100, threshold$Nplus, threshold$Nminus))
    return(threshold)
  }
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
# %% calculate network
for (i in timestamps){
  # Create the directory structure if it doesn't exist
  plot_dir <- file.path("Plots", "BigDataDaysFilter")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  plot_path <- file.path(plot_dir, paste0("Day_", i))
  data_list <- data_list_times[[i]]
  cat('Calculating network on day ',i,'\n')
  network_results <- stabENG(data_list, labels = shared_otu, var.thresh = 0.1, rep.num = 25,
    nlambda1=20,lambda1.min=0.01,lambda1.max=1,nlambda2=20,lambda2.min=0,lambda2.max=0.1,
    lambda2.init=0.01,ebic.gamma=0.2, nCores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
  network_list_raw[[i]]$Nplus <- network_results$opt.fit$Nplus # precision matrix estimates
  network_list_raw[[i]]$Nminus <- network_results$opt.fit$Nminus # precision matrix estimates
  diag(network_list_raw[[i]]$Nplus) = diag(network_list_raw[[i]]$Nminus) <- 0
  network_pcor_raw[[i]]$Nplus <- network_results$opt.fit.pcor$Nplus
  network_pcor_raw[[i]]$Nminus <- network_results$opt.fit.pcor$Nminus
  diag(network_pcor_raw[[i]]$Nplus) = diag(network_pcor_raw[[i]]$Nminus) <- 0
  lambda1[[i]] <- network_results$opt.lambda1
  lambda2[[i]] <- network_results$opt.lambda2
  if (perm_filter == TRUE) {
    cat('\nPerforming permutation test for network on Day',i,'\n')
    perm_threshold <- permutation_threshold(
      data_list = data_list, 
      shared_otu = shared_otu, 
      n_perm = 50,
      alpha = 0.1,
      lambda1 = lambda1[[i]],
      lambda2 = lambda2[[i]],
      plot_threshold = TRUE,
      plot_path = plot_path
    )
    # filter edge sparsity for simulation
    network_list[[i]]$Nplus <- network_list_raw[[i]]$Nplus
    network_list[[i]]$Nminus <- network_list_raw[[i]]$Nminus
    network_list[[i]]$Nplus[abs(network_list[[i]]$Nplus) < perm_threshold$Nplus] <- 0
    network_list[[i]]$Nminus[abs(network_list[[i]]$Nminus) < perm_threshold$Nminus] <- 0
    # filter network_pcor
    network_pcor[[i]]$Nplus <- network_pcor_raw[[i]]$Nplus
    network_pcor[[i]]$Nminus <- network_pcor_raw[[i]]$Nminus
    network_pcor[[i]]$Nplus[abs(network_pcor[[i]]$Nplus) < perm_threshold$Nplus] <- 0
    network_pcor[[i]]$Nminus[abs(network_pcor[[i]]$Nminus) < perm_threshold$Nminus] <- 0
  }
  else {
    cat('\nSkipping permutation test filtering on Day',i,'\n')
    # Use raw pcor
    network_list[[i]] <- network_list_raw[[i]]
    network_pcor[[i]] <- network_pcor_raw[[i]]
    # Ensure diagonals are zero
    diag(network_list[[i]]$Nplus) = diag(network_list[[i]]$Nminus) <- 0
    diag(network_pcor[[i]]$Nplus) = diag(network_pcor[[i]]$Nminus) <- 0
    # Define thresholds as 0 if not calculated, for potential use later (e.g., sim_perm_filter)
    perm_threshold_Nplus <- 0
    perm_threshold_Nminus <- 0
  }
  # %% Plot network on Phylum level
  Phylum_groups <- as.factor(otu_tax[rownames(network_pcor[[i]]$Nplus),"Phylum"])
  png(filename=paste0(plot_path,"_network_Nplus_Phylum_Stab_Filtered_vsized.png"))
  qgraph::qgraph(network_pcor[[i]]$Nplus, 
    layout = "circle",
    edge.color = ifelse(network_pcor[[i]]$Nplus > 0, "blue", "red"),
    title = "Stab Network Nplus by Phylum",
    vsize = 2.5,
    groups = Phylum_groups)
  dev.off()
  
  png(filename=paste0(plot_path,"_network_Nminus_Phylum_Stab_Filtered_vsized.png"))
  qgraph::qgraph(network_pcor[[i]]$Nminus, 
    layout = "circle",
    edge.color = ifelse(network_pcor[[i]]$Nminus > 0, "blue", "red"),
    title = "Stab Network Nminus by Phylum",
    vsize = 2.5,
    groups = Phylum_groups)
  dev.off()
  # %%Visualize Edge weights
  cor_values_Nplus <- as.vector(network_pcor[[i]]$Nplus)
  cor_values_Nminus <- as.vector(network_pcor[[i]]$Nminus)
  cor_values_Nplus_filtered <- cor_values_Nplus[cor_values_Nplus != 0]
  cor_values_Nminus_filtered <- cor_values_Nminus[cor_values_Nminus != 0]

  ggplot(data.frame(Type = c(rep("N+", length(cor_values_Nplus_filtered)), rep("N-", length(cor_values_Nminus_filtered))), 
    Value = c(cor_values_Nplus_filtered, cor_values_Nminus_filtered)), 
    aes(x = Value, color = Type)) + 
    geom_density() + 
    labs(x = "Correlation Value", y = "Density")
  ggsave(filename=paste0(plot_path,"_correlation_distribution.png"))
  
  # Visualize network_pcor[[i]]$Nplus (Circular Layout)
  # if (!is.null(network_pcor[[i]]$Nplus) && nrow(network_pcor[[i]]$Nplus) > 0) {
  #   cir_plot(network_pcor[[i]]$Nplus, otu_tax, shared_otu, paste0(plot_path, "_network_Nplus_Circular"))
  # } else {
  #   warning("network_pcor[[i]]$Nplus is empty or invalid.")
  # }
  # cir_plot(network_pcor[[i]]$Nminus,otu_tax, shared_otu, paste0(plot_path,"_network_Nminus_Circular"))
  cat('Synthesize simulation count data on day ',i,'\n')
  for (j in 1:n_sim)
  {
    Sim_list[[i]][[j]] <- list(
    Nplus = synthesize_scaled_data(data_list$Nplus, network_list[[i]]$Nplus),
    Nminus = synthesize_scaled_data(data_list$Nminus, network_list[[i]]$Nminus))
  }

  cat('Calculate simulation data network on day ',i,' start.\n')
  for (j in 1:n_sim)
  {
    Res_sim[[i]][[j]] <- stabENG(Sim_list[[i]][[j]], labels = shared_otu, var.thresh = 0.1, rep.num = 25,
      nlambda1=20,lambda1.min=0.01,lambda1.max=1,nlambda2=20,lambda2.min=0,lambda2.max=0.1,
      lambda2.init=0.01,ebic.gamma=0.2, nCores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
    cat('Simulation Progress on Day',i,':', 100*j/n_sim,'% done.\n')
    # filter edge sparsity
    if (sim_perm_filter == TRUE)
      {
        perm_threshold_Nplus <- permutation_threshold(
          data_matrix = Sim_list[[i]][[j]]$Nplus, 
          pcor_matrix = Res_sim[[i]][[j]]$opt.fit.pcor$Nplus, 
          n_perm = 50,
          alpha = 0.05,
          plot_threshold = FALSE,
          plot_path = NULL
        )
        perm_threshold_Nminus <- permutation_threshold(
          data_matrix = Sim_list[[i]][[j]]$Nminus, 
          pcor_matrix = Res_sim[[i]][[j]]$opt.fit.pcor$Nminus, 
          n_perm = 50,
          alpha = 0.05,
          plot_threshold = FALSE,
          plot_path = NULL
        )
      }
    Res_sim[[i]][[j]]$opt.fit.pcor$Nplus[abs(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus) < perm_threshold_Nplus] <- 0
    Res_sim[[i]][[j]]$opt.fit.pcor$Nminus[abs(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus) < perm_threshold_Nminus] <- 0
    diag(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus) = diag(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus) <- 0
    # plot sim network
    png(filename=paste0(plot_path,"_simulation_network_",j,"_Nplus_Phylum_Stab_Filtered_vsized.png"))
    qgraph::qgraph(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus, 
      layout = "circle",
      edge.color = ifelse(Res_sim[[i]][[j]]$opt.fit.pcor$Nplus > 0, "blue", "red"),
      title = "Stab Network Nplus by Phylum",
      vsize = 2.5,
      groups = Phylum_groups)
    dev.off()
    
    png(filename=paste0(plot_path,"_simulation_network_",j,"_Nminus_Phylum_Stab_Filtered_vsized.png"))
    qgraph::qgraph(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus, 
      layout = "circle",
      edge.color = ifelse(Res_sim[[i]][[j]]$opt.fit.pcor$Nminus > 0, "blue", "red"),
      title = "Stab Network Nminus by Phylum",
      vsize = 2.5,
      groups = Phylum_groups)
    dev.off()
  }
  Sim_adj <- list()
  for (j in 1:n_sim)
  {
    Sim_adj[[i]][[j]] <- list(
      Nplus = (Res_sim[[i]][[j]]$opt.fit.pcor$Nplus !=0)*1,
      Nminus = (Res_sim[[i]][[j]]$opt.fit.pcor$Nminus !=0)*1
    )
  }
  #%% Confusion matrices
  calculate_mcc <- function(tp, tn, fp, fn)
    {
      numerator <- (tp * tn) - (fp * fn)
      denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  
      # If the denominator is 0, return 0 to avoid NaN errors
      if (denominator == 0) {
        return(0)
      } else {
        return(numerator / denominator)
      }
    }
  calculate_metrics <- function(true_adj, sim_adj)
    {
      # Flatten the matrices into vectors (upper triangular part, excluding diagonal)
      true_edges <- as.vector(true_adj[upper.tri(true_adj)])
      sim_edges <- as.vector(sim_adj[upper.tri(sim_adj)])
    
      # Confusion matrix
      cm <- confusionMatrix(as.factor(sim_edges), as.factor(true_edges), positive = "1")
      tn <- as.numeric(cm$table[1,1]) #true negatives
      fp <- as.numeric(cm$table[1,2]) #false positives
      fn <- as.numeric(cm$table[2,1]) #false negatives
      tp <- as.numeric(cm$table[2,2]) #true positives
    
      # Calculate TPR, FPR, Precision, Recall
      tpr <- tp / (tp + fn)  # Sensitivity / Recall
      fpr <- fp / (fp + tn)  # 1 - Specificity
      precision <- tp / (tp + fp)
      recall <- tpr
    
      # Calculate ROC and AUC
      roc_obj <- pROC::roc(true_edges, sim_edges)
      auc <- as.numeric(auc(roc_obj)) # plot here
    
      # F1 Score
      f1 <- F1_Score(sim_edges, true_edges)
    
      # MCC (Matthews correlation coefficient)
      mcc <- calculate_mcc(tp, tn, fp, fn)
    
      # Return metrics as a list
      return(list(TPR = tpr, FPR = fpr, Precision = precision, Recall = recall,
                  F1 = f1, AUC = auc, MCC = mcc))
    }
  cat('Calculate confusion matrices on day ',i,'\n')
  confusion_results <- lapply(1:n_sim, function(j)
    {
      true_adj_Nplus <- (network_pcor[[i]]$Nplus !=0)*1
      true_adj_Nminus <- (network_pcor[[i]]$Nminus !=0)*1
      Nplus_metrics <- calculate_metrics(true_adj_Nplus, Sim_adj[[i]][[j]]$Nplus)
      Nminus_metrics <- calculate_metrics(true_adj_Nminus, Sim_adj[[i]][[j]]$Nminus)
      return(list(Nplus = Nplus_metrics, Nminus = Nminus_metrics))
  })
  
  
  results_df <- do.call(rbind, lapply(confusion_results, as.data.frame))
  
  results_df_long <- results_df %>%
    dplyr::select(starts_with("Nplus.") | starts_with("Nminus.")) %>%
    tidyr::pivot_longer(cols = everything(), 
                names_to = c("group", "metric"), 
                names_sep = "\\.",
                values_to = "value") %>%
    dplyr::mutate(matrix_id = rep(1:n_sim, each = 14)) %>%
    dplyr::mutate(times = rep(i, n_sim*14))
  # append results_df_long to confusion_results_df
  confusion_results_df <- rbind(confusion_results_df, results_df_long)
  # # %% boxplot
  # metrics_to_plot <- c("TPR", "FPR", "Precision", "Recall", "F1", "AUC", "MCC")
  # filtered_df <- results_df_long %>%
  #   filter(metric %in% metrics_to_plot)
  
  # # Create boxplots
  # p <- ggplot(filtered_df, aes(x = group, y = value, color = group)) +
  #   geom_boxplot() +
  #   facet_wrap(~ metric, nrow = 1, scales = "free_y") + # facet_wrap for one row
  #   labs(title = "Distribution of Confusion Metrics",
  #        x = "Group",
  #        y = "Value",
  #        color = "Group") +
  #   theme_minimal()
  # ggsave(filename = paste0(plot_path,"_Filtered_Confusion_boxplot.png"), p)
}

# Filter for TPR, FPR, and F1 metrics
selected_metrics <- c("TPR", "FPR", "F1")
filtered_data <- confusion_results_df %>% filter(metric %in% selected_metrics)

# Plot
p <- ggplot(filtered_data, aes(x = factor(as.integer(times)), y = value, fill = group)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "black", position = position_dodge(0.8)) +
  scale_fill_manual(values = c("Nplus" = "#00BFC4", "Nminus" = "#F8766D")) +  # Custom colors
  facet_wrap(~metric, scales = "free_y") +  # Create separate plots for each metric, free y-scales!
  labs(x = "Times", y = "Value", title = "Confusion Metrics") +
  theme_minimal(base_size = 14) +  # Clean theme with adjusted font size
  theme(
    legend.position = "top",  # Move legend to top
    panel.grid.major = element_line(color = "gray90"),  # Add grid lines
    panel.border = element_rect(color = "black", fill = NA),  # Add border
    strip.text = element_text(face = "bold", size = 12)  # Make facet labels bold
  )
ggsave(filename = "//home//14720078//ProjectCode//Plots//BigDataDaysFilter//Filtered_Confusion_boxplot.png", p)
# save network_list and network_pcor as csv
write.csv(network_pcor, "//home//14720078//ProjectCode//DataImage//network_pcor.csv")
save.image("//home//14720078//ProjectCode//DataImage//big1226_Days_network_results_Big_Days_Filtered.RData")
cat("All done in: ", Sys.time() - start_time, "\n")