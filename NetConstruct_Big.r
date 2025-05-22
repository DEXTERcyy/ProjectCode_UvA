library(phyloseq)
library(SPRING)
library(SpiecEasi)
library(caret)
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
start_time <- Sys.time()
rawdata <- readRDS("//home//14720078//ProjectCode//data//DavarData1_substrate_phyloseq_1226_final_filtered.rds")
otu_raw <- otu_table(rawdata)
otu_RA <- transform_sample_counts(otu_raw, function(x) x / sum(x) )
otu_RA <- filter_taxa(otu_RA, function(x) mean(x) > 1e-3, TRUE)
shared_otu <- rownames(otu_RA)
otu_Ab <- as.data.frame(t(otu_raw[shared_otu,]))
sam_info <- as.data.frame(sample_data(rawdata))
otu_tax <- as.data.frame(tax_table(rawdata))

# %%split otu_Ab by condition group
otu_Ab_Nplus <- otu_Ab[rownames(otu_Ab) %in% rownames(sam_info[sam_info$growthCondition=="plusN",]),]
otu_Ab_Nminus <- otu_Ab[rownames(otu_Ab) %in% rownames(sam_info[sam_info$growthCondition=="minusN",]),]
data_list <- list(Nplus = otu_Ab_Nplus, Nminus = otu_Ab_Nminus)
# by Group and Days
timestamps <- as.character(sort(as.integer(levels(sam_info$Days))))
otu_Ab_Nplus_times <- list()
otu_Ab_Nminus_times <- list()
data_list_times <- list()
taxa_filtering <- TRUE
n_sim <- 5
for (i in timestamps)
  {
    otu_Ab_Nplus_times[[i]] <- otu_Ab[rownames(otu_Ab) %in% 
      rownames(sam_info[sam_info$growthCondition=="plusN" & sam_info$Days == i,]),]
    otu_Ab_Nminus_times[[i]] <- otu_Ab[rownames(otu_Ab) %in%
      rownames(sam_info[sam_info$growthCondition=="minusN" & sam_info$Days == i,]),]
    data_list_times[[i]] <- list(Nplus = otu_Ab_Nplus_times[[i]], Nminus = otu_Ab_Nminus_times[[i]])
    if (taxa_filtering == TRUE)
    {
      # Sum the number of zeros in each column and sort by the number of zeros
      zero_counts_Nplus <- colSums(otu_Ab_Nplus_times[[i]] == 0)
      zero_counts_Nminus <- colSums(otu_Ab_Nminus_times[[i]] == 0)
      
      # Select the top 100 columns with the least zeros for each group
      sorted_Nplus <- names(sort(zero_counts_Nplus, decreasing = FALSE))[1:100]
      sorted_Nminus <- names(sort(zero_counts_Nminus, decreasing = FALSE))[1:100]
      
      # Use the intersecting column names for Nplus & Nminus
      common_taxa[[i]] <- intersect(sorted_Nplus, sorted_Nminus)
    }
}
if (taxa_filtering == TRUE)
{
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
  # Filter the shared_otu to only include the common taxa
  shared_otu <- shared_otu[shared_otu %in% common_taxa_intersection]
  # Filter the otu_tax to only include the common taxa
  otu_tax <- otu_tax[rownames(otu_tax) %in% common_taxa_intersection, ]
  if (!all(common_taxa_intersection %in% rownames(otu_tax))) {
    warning("Some taxa in common_taxa_intersection are not in otu_tax.")
  }
  if (!all(common_taxa_intersection %in% shared_otu)) {
    warning("Some taxa in common_taxa_intersection are not in shared_otu.")
  }
}
# %%
network_list_raw <- list()
network_list <- list()
network_pcor_raw <- list()
network_pcor <- list()
confusion_results_df <- list()
for (i in timestamps)
{
  plot_path = paste0("//home//14720078//ProjectCode//Plots//BigDataDaysFilter//Day_",i)
  data_list <- data_list_times[[i]]
  cat('Calculating network on day ',i,'\n')
  network_results <- stabENG(data_list, labels = shared_otu, var.thresh = 0.1, rep.num = 25,
    nlambda1=20,lambda1.min=0.01,lambda1.max=1,nlambda2=20,lambda2.min=0,lambda2.max=0.1,
    lambda2.init=0.01,ebic.gamma=0.2, nCores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
  network_list_raw[[i]]$Nplus <- network_results$opt.fit$Nplus # precision matrix estimates
  network_list_raw[[i]]$Nminus <- network_results$opt.fit$Nminus # precision matrix estimates
  network_pcor_raw[[i]]$Nplus <- network_results$opt.fit.pcor$Nplus
  network_pcor_raw[[i]]$Nminus <- network_results$opt.fit.pcor$Nminus
  # filter edge sparsity for simulation
  network_list[[i]]$Nplus <- network_list_raw[[i]]$Nplus
  network_list[[i]]$Nplus[abs(network_list[[i]]$Nplus) < 0.01] <- 0
  network_list[[i]]$Nminus <- network_list_raw[[i]]$Nminus
  network_list[[i]]$Nminus[abs(network_list[[i]]$Nminus) < 0.1] <- 0
  diag(network_list[[i]]$Nplus) = diag(network_list[[i]]$Nminus) <- 0
  # filter network_pcor
  network_pcor[[i]]$Nplus <- network_pcor_raw[[i]]$Nplus
  network_pcor[[i]]$Nplus[abs(network_pcor[[i]]$Nplus) < 0.01] <- 0
  network_pcor[[i]]$Nminus <- network_pcor_raw[[i]]$Nminus
  network_pcor[[i]]$Nminus[abs(network_pcor[[i]]$Nminus) < 0.01] <- 0
  # network_pcor[[i]]$Nplus[abs(network_pcor[[i]]$Nplus) < 0.01] <- 0
  # network_pcor[[i]]$Nminus[abs(network_pcor[[i]]$Nminus) < 0.01] <- 0
  diag(network_pcor[[i]]$Nplus) = diag(network_pcor[[i]]$Nminus) <- 0
  # %% Plot network on Phylum level
  Phylum_groups <- as.factor(otu_tax[rownames(network_list[[i]]$Nplus),"Phylum"])
  png(filename=paste0(plot_path,"_network_Nplus_Phylum_Stab_Filtered_vsized.png"))
  qgraph::qgraph(network_list[[i]]$Nplus, 
    layout = "circle",
    edge.color = ifelse(network_list[[i]]$Nplus > 0, "blue", "red"),
    title = "Stab Network Nplus by Phylum",
    vsize = 2.5,
    groups = Phylum_groups)
  dev.off()
  
  png(filename=paste0(plot_path,"_network_Nminus_Phylum_Stab_Filtered_vsized.png"))
  qgraph::qgraph(network_list[[i]]$Nminus, 
    layout = "circle",
    edge.color = ifelse(network_list[[i]]$Nminus > 0, "blue", "red"),
    title = "Stab Network Nminus by Phylum",
    vsize = 2.5,
    groups = Phylum_groups)
  dev.off()
  # %%Visualize Edge weights (pCor or Prec?)
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
  
  # Visualize network_list[[i]]$Nplus (Circular Layout)
  cir_plot(network_list[[i]]$Nplus,otu_tax, shared_otu, paste0(plot_path,"_network_Nplus_Circular"))
  cir_plot(network_list[[i]]$Nminus,otu_tax, shared_otu, paste0(plot_path,"_network_Nminus_Circular"))
  #%%
  set.seed(10010)
  synthesize_scaled_data <- function(dat, net)
    {
      graph <- (net != 0)*1     # SpiecEasi::make_graph('cluster', dat$n_OTUs, dat$n_edge)
      attr(graph, "class") <- "graph"
      Prec <- SpiecEasi::graph2prec(graph)
      Cor <- cov2cor(SpiecEasi::prec2cov(Prec))
      X <- SpiecEasi::synth_comm_from_counts(dat, mar = 2, distr = 'zinegbin', Sigma = Cor, n = nrow(dat))
      return(X)
    }
  cat('Synthesize simulation data on day ',i,'\n')
  Sim_list <- list()
  for (j in 1:n_sim)
    {
      Sim_list[[i]][[j]] <- list(
        Nplus = synthesize_scaled_data(otu_Ab_Nplus, network_list[[i]]$Nplus),
        Nminus = synthesize_scaled_data(otu_Ab_Nminus, network_list[[i]]$Nminus)
      )
    }
  cat('Calculate simulation data network on day ',i,' start.\n')
}