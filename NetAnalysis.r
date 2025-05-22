# %%
myNetAnalysis <- function(assoMat1, assoMat2, sparsMethod = "t-test")
  {
    set.seed(10010)
    dataType = assoType = "partialCorr"
    sampleSize = c(118,127) # 118 Nplus 117 Nminus
    counts1 <- counts2 <- NULL
    countsJointOrig <- NULL
    countsOrig1 <- countsOrig2 <- NULL
    groups <- NULL
    dissFunc = "signed"
    dissFuncPar = NULL
    simFunc = simFuncPar = NULL
    weighted = TRUE
    measure = "spieceasi"
    measurePar = NULL
    adjust = "none" #"lfdr"
    trueNullMethod = "convest"
    thresh = c(0.3,0.3)
    alpha = c(0.05,0.05) # p-value 0.1
    lfdrThresh = c(0.2,0.2)
    nboot = 1000L
    assoBoot = FALSE
    softThreshType = "signed"
    softThreshPower = c(NULL,NULL)
    softThreshCut = c(0.8,0.8)
    kNeighbor = 3L
    knnMutual = FALSE
    seed = 10010
    sparsReslt <- .sparsify(assoMat = assoMat1,
      countMat = counts1,
      sampleSize = sampleSize[1],
      measure = measure,
      measurePar = measurePar,
      assoType = assoType,
      sparsMethod = sparsMethod,
      thresh = thresh[1],
      alpha = alpha[1],
      adjust = adjust,
      lfdrThresh = lfdrThresh[1],
      trueNullMethod = trueNullMethod,
      nboot = nboot,
      assoBoot = assoBoot,
      cores = 18,
      softThreshType = softThreshType,
      softThreshPower = softThreshType[1],
      softThreshCut = softThreshCut[1],
      logFile = "log.txt",
      kNeighbor = kNeighbor,
      knnMutual = knnMutual,
      verbose = FALSE,
      seed = seed)

    assoEst1 <- assoMat1
    assoMat1 <- sparsReslt$assoNew
    power1 <- sparsReslt$power
    dissEst1 <- dissScale1 <- NULL

    dissMat1 <- .transToDiss(x = assoMat1,  dissFunc= dissFunc,
            dissFuncPar = dissFuncPar)

    if (sparsMethod == "softThreshold") {
      simMat1 <- sparsReslt$simMat
      adjaMat1 <- assoMat1
      assoBoot1 <- NULL

    } else {
      simMat1 <- .transToSim(x = dissMat1, simFunc = simFunc,
              simFuncPar = simFuncPar)

      adjaMat1 <- .transToAdja(x = simMat1, weighted = weighted)

      if (sparsMethod == "bootstrap" && !is.null(assoBoot) &&
        !is.list(assoBoot) && assoBoot == TRUE) {
        assoBoot1 <- sparsReslt$assoBoot
      } else {
      assoBoot1 <- NULL
      }
    }

    sparsReslt <- .sparsify(assoMat = assoMat2,
      countMat = counts2,
      sampleSize = sampleSize[2],
      measure = measure,
      measurePar = measurePar,
      assoType = assoType,
      sparsMethod = sparsMethod,
      thresh = thresh[2],
      alpha = alpha[2],
      adjust = adjust,
      lfdrThresh = lfdrThresh[2],
      trueNullMethod = trueNullMethod,
      nboot = nboot,
      assoBoot = assoBoot,
      cores = 18,
      softThreshType = softThreshType,
      softThreshPower = softThreshType[2],
      softThreshCut = softThreshCut[2],
      logFile = "log.txt",
      kNeighbor = kNeighbor,
      knnMutual = knnMutual,
      verbose = FALSE,
      seed = seed)

    assoEst2 <- assoMat2
    assoMat2 <- sparsReslt$assoNew
    power2 <- sparsReslt$power
    dissEst2 <- dissScale2 <- NULL

    dissMat2 <- .transToDiss(x = assoMat2, dissFunc = dissFunc,
              dissFuncPar = dissFuncPar)

    if (sparsMethod == "softThreshold") {
      simMat2 <- sparsReslt$simMat
      adjaMat2 <- assoMat2
      assoBoot2 <- NULL
    } else {
      simMat2 <- .transToSim(x = dissMat2, simFunc = simFunc,
              simFuncPar = simFuncPar)

      adjaMat2 <- .transToAdja(x = simMat2, weighted = weighted)

      if (sparsMethod == "bootstrap" && !is.null(assoBoot) &&
        !is.list(assoBoot) && assoBoot == TRUE) {
        assoBoot2 <- sparsReslt$assoBoot}
      else {
        assoBoot2 <- NULL
      }
    }
    # Create edge list
    g <- igraph::graph_from_adjacency_matrix(adjaMat1, weighted = weighted,
      mode = "undirected", diag = FALSE)

    if (is.null(igraph::E(g)$weight)) {
      isempty1 <- TRUE
      edgelist1 <- NULL
    } else {
      isempty1 <- FALSE

      edgelist1 <- data.frame(igraph::get.edgelist(g))
      colnames(edgelist1) <- c("v1", "v2")
      if (!is.null(assoMat1)) {
        edgelist1$asso <- sapply(1:nrow(edgelist1), function(i) {
          assoMat1[edgelist1[i, 1], edgelist1[i, 2]]
        })
      }

      edgelist1$diss <- sapply(1:nrow(edgelist1), function(i) {
        dissMat1[edgelist1[i, 1], edgelist1[i, 2]]
      })

      if (all(adjaMat1 %in% c(0,1))) {
        edgelist1$sim <-sapply(1:nrow(edgelist1), function(i) {
          simMat1[edgelist1[i, 1], edgelist1[i, 2]]
        })
      }

      edgelist1$adja <- sapply(1:nrow(edgelist1), function(i) {
        adjaMat1[edgelist1[i, 1], edgelist1[i, 2]]
      })
    }

    # Create edge list
    g <- igraph::graph_from_adjacency_matrix(adjaMat2, weighted = weighted, 
          mode = "undirected", diag = FALSE)
  
    if (is.null(igraph::E(g)$weight)) {
      isempty2 <- TRUE
      edgelist2 <- NULL
    } else {
      isempty2 <- FALSE

      edgelist2 <- data.frame(igraph::get.edgelist(g))
      colnames(edgelist2) <- c("v1", "v2")

      if (!is.null(assoMat2)) {
        edgelist2$asso <- sapply(1:nrow(edgelist2), function(i) {
          assoMat2[edgelist2[i, 1], edgelist2[i, 2]]
        })
      }

      edgelist2$diss <- sapply(1:nrow(edgelist2), function(i) {
        dissMat2[edgelist2[i, 1], edgelist2[i, 2]]
      })

      if (all(adjaMat2 %in% c(0, 1))) {
        edgelist2$sim <- sapply(1:nrow(edgelist2), function(i) {
          simMat2[edgelist2[i, 1], edgelist2[i, 2]]
        })
      }

      edgelist2$adja <- sapply(1:nrow(edgelist2), function(i) {
        adjaMat2[edgelist2[i, 1], edgelist2[i, 2]]
      })
    }

    if (isempty1) {
      message("\nNetwork 1 has no edges.")
    }
    if (isempty2) {
      message("Network 2 has no edges.")
    }

    #=============================================================================
    output <- list()

    output$edgelist1 <- edgelist1
    output$edgelist2 <- edgelist2
    output$assoMat1 <- assoMat1
    output$assoMat2 <- assoMat2
    output$dissMat1 <- dissMat1
    output$dissMat2 <- dissMat2
    output$simMat1 <- simMat1
    output$simMat2 <- simMat2
    output$adjaMat1 <- adjaMat1
    output$adjaMat2 <- adjaMat2

    output$assoEst1 <- assoEst1
    output$assoEst2 <- assoEst2
    output$dissEst1 <- dissEst1
    output$dissEst2 <- dissEst2
    output$dissScale1 <- dissScale1
    output$dissScale2 <- dissScale2

    output$assoBoot1 <- assoBoot1
    output$assoBoot2 <- assoBoot2

    output$countMat1 <- countsOrig1
    output$countMat2 <- countsOrig2
    if (!is.null(countsJointOrig)) output$countsJoint <- countsJointOrig
    output$normCounts1 <- counts1
    output$normCounts2 <- counts2
    output$groups <- groups
    output$sampleSize <- sampleSize
    output$softThreshPower <- list(power1 = power1, 
      power2 = power2) # calculated power
    output$assoType <- assoType
    output$parameters <- list(
      dataType = dataType,
      # group = group,
      # filtTax = filtTax,
      # filtTaxPar = filtTaxPar,
      # filtSamp = filtSamp,
      # filtSampPar = filtSampPar,
      # jointPrepro = jointPrepro,
      # zeroMethod = zeroMethod,
      # zeroPar = zeroPar,
      # needfrac = needfrac,
      # needint = needint,
      # normMethod = normMethod,
      # normPar = normPar,
      measure = measure,
      measurePar = measurePar,
      sparsMethod = sparsMethod,
      thresh = thresh,
      adjust = adjust,
      trueNullMethod = trueNullMethod,
      alpha = alpha,
      lfdrThresh = lfdrThresh,
      nboot = nboot,
      softThreshType = softThreshType,
      softThreshPower = softThreshPower,
      softThreshCut = softThreshCut,
      kNeighbor = kNeighbor,
      knnMutual = knnMutual,
      dissFunc = dissFunc,
      dissFuncPar = dissFuncPar,
      simFunc = simFunc,
      simFuncPar = simFuncPar,
      # scaleDiss = scaleDiss,
      weighted = weighted,
      sampleSize = sampleSize)
    output$call = match.call()
    output$twoNets <- TRUE
    class(output) <- "microNet"
    return(output)
}
# %%
library(NetCoMi)
load(file = "ProjectCode/DataImage/big1226_Days_network_results_Big_Days_Filtered.RData")
source(file="ProjectCode/Packages/NetCoMi/R/dot-sparsify.R")
source(file="ProjectCode/Packages/NetCoMi/R/transform.R")
# %%
NetCon <- myNetAnalysis(network_pcor$`4`$Nplus,network_pcor$`4`$Nminus)
NetRes <- netAnalyze(NetCon, clustMethod = "cluster_fast_greedy")
NetSum <- summary(NetRes)

# %% GCM heatmap plot
pdf(file="Plots/NetCoMi/GCM_heatmap_plot.pdf")
plotHeat(mat = NetRes$graphletLCC$gcm1,
  pmat = NetRes$graphletLCC$pAdjust1,
  type = "mixed",
  title = "GCM", 
  colorLim = c(-1, 1),
  mar = c(2, 0, 2, 0))
graphics::rect(xleft   = c( 0.5,  1.5, 4.5,  7.5),# Add rectangles highlighting the four types of orbits
        ybottom = c(11.5,  7.5, 4.5,  0.5),
        xright  = c( 1.5,  4.5, 7.5, 11.5),
        ytop    = c(10.5, 10.5, 7.5,  4.5),
        lwd = 2, xpd = NA)

text(6, -0.2, xpd = NA, 
"Significance codes:  ***: 0.001;  **: 0.01;  *: 0.05")
dev.off()
# %%Visualizing the network
pdf(file="Plots/NetCoMi/network_plot.pdf",width = 50, height = 30)
p <- plot(NetRes, 
  nodeColor = "cluster", 
  nodeSize = "eigenvector",
  title1 = "Network on OTU level", 
  showTitle = TRUE,
  cexTitle = 5)
legend(0.7, 1.1, cex = 5, title = "estimated association:",
legend = c("+","-"), lty = 2, lwd = 5, col = c("#009900","red"), 
bty = "n", horiz = TRUE)
dev.off()

# %% diff net
diff_season <- diffnet(NetCon,n1 = 29,n2 = 30,
  diffMethod = "fisherTest", 
  adjust = "none")

# %% net compare
tmp = summary(NetRes, 
  groupNames = c("Nplus", "Nminus"),
  showCentr = c("degree", "between", "closeness"), 
  numbNodes = 5)