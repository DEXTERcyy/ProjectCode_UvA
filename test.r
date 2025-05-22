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
load("//home//14720078//ProjectCode//DataImage//big1226_Days_network_results_Big_Days_Filtered.RData")
synthesize_scaled_data <- function(dat, net)
{
cat("  Start synthesize_scaled_data...\n") # Add print statement
graph <- (net != 0)*1     # SpiecEasi::make_graph('cluster', dat$n_OTUs, dat$n_edge)
attr(graph, "class") <- "graph"
cat("    Calling SpiecEasi::graph2prec...\n") # Add print statement
Prec <- SpiecEasi::graph2prec(graph)
cat("    Calling SpiecEasi::prec2cov...\n") # Add print statement
Cor <- cov2cor(SpiecEasi::prec2cov(Prec))
cat("    Calling SpiecEasi::synth_comm_from_counts...\n") # Add print statement
X <- SpiecEasi::synth_comm_from_counts(dat, mar = 2, distr = 'zinegbin', Sigma = Cor, n = nrow(dat))
cat("  Finished synthesize_scaled_data.\n") # Add print statement
return(X)
}
for (i in timestamps)
{
    data_list <- data_list_times[[i]]
    cat('Synthesize simulation count data on day ',i,'\n')
    for (j in 1:n_sim)
    {
        Sim_list[[i]][[j]] <- list(
        Nplus = synthesize_scaled_data(data_list$Nplus, network_list[[i]]$Nplus),
        Nminus = synthesize_scaled_data(data_list$Nminus, network_list[[i]]$Nminus))
    }
}