## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  warning=FALSE, message=FALSE,
  collapse = TRUE,
  fig.align = "center",
  prompt = TRUE,
  comment = ""
)

## ----demo-plot, fig.width=8, fig.height=6, echo=FALSE-------------------------
knitr::include_graphics("exampleFig.png")

## ----eval=FALSE---------------------------------------------------------------
#  devtools::install_github("KiRinHong/miMediation")

## -----------------------------------------------------------------------------
library(miMediation)

## -----------------------------------------------------------------------------
# Load data
data(data.cecal)
# Take a look at the data
Trt <- data.cecal$treatment
table(Trt) # 0: control 1: antibotics
M <- data.cecal$mediators
head(M[,1:6])
Y <- data.cecal$outcome
summary(Y)
tree <- data.cecal$tree
tree

## ----tree-plot,fig.width=8, fig.height=8, fig.show='hold'---------------------
# set random seed here so that you can get the same result every time you run the code
set.seed(123)
cecal.rsltlst <- phyloMed(Trt, M, Y, tree, fdr.alpha = 0.1, n.perm = 1e5, graph = TRUE)

# take a look at phyloseq-class object
cecal.physeq <- cecal.rsltlst$clean.data
cecal.physeq
cecal.rslt <- cecal.rsltlst$rslt
# take a look at rslt (PhyloMed.P)
cecal.rslt$PhyloMed.P

