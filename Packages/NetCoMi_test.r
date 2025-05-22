# %%
library(NetCoMi)
data("amgut1.filt")
data("amgut2.filt.phy")
set.seed(123456)
group <- sample(1:2, nrow(amgut1.filt), replace = TRUE)

amgut1 <- amgut1.filt[group == 1, ]
amgut2 <- amgut1.filt[group == 2, ]

net3 <- netConstruct(data = amgut1, 
                     data2 = amgut2,
                     measure = "pearson",
                     filtTax = "highestVar",
                     filtTaxPar = list(highestVar = 50),
                     filtSamp = "totalReads",
                     filtSampPar = list(totalReads = 1000),
                     zeroMethod = "multRepl", 
                     normMethod = "clr",
                     sparsMethod = "t-test")

props3 <- netAnalyze(net3, clustMethod = "cluster_fast_greedy")

plot(props3, sameLayout = TRUE)
netCompare()