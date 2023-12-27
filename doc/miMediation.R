## ----include = FALSE----------------------------------------------------------
# pdf_document
# prettydoc::html_pretty:
#  theme: cayman
#  highlight: github
knitr::opts_chunk$set(
  warning=FALSE, message=FALSE,
  collapse = TRUE,
  fig.align = "center",
  prompt = TRUE,
  comment = ""
)

## ----demo-plot, fig.width=8, fig.height=4, echo=FALSE-------------------------
knitr::include_graphics("exampleFig.png")

## -----------------------------------------------------------------------------
library(miMediation)
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

## ----phy-tree-plot,fig.width=8, fig.height=8, fig.show='hold'-----------------
# set random seed here so that you can get the same result every time you run the code
set.seed(123)
cecal.rsltlst <- phyloMed(Trt, M, Y, tree = tree, fdr.alpha = 0.1, 
                          n.perm = 1e4, graph = "rectangular")

# take a look at phyloseq-class object
cecal.physeq <- cecal.rsltlst$clean.data
cecal.physeq
cecal.rslt <- cecal.rsltlst$rslt
# take a look at rslt (PhyloMed.P)
cecal.rslt$PhyloMed.P

## -----------------------------------------------------------------------------
# Load data
data(data.zeeviD)
# Take a look at the data
Trt <- data.zeeviD$treatment
table(Trt) # 0: control 1: treatment
M <- data.zeeviD$mediators
dim(M)
Y <- data.zeeviD$outcome
summary(Y)
tree <- data.zeeviD$tree
head(tree)

## ----tax-tree-plot,fig.width=8, fig.height=8, fig.show='hold',message=FALSE,warning=FALSE----
# run asymptotic result by default
demo.rsltlst <- phyloMed(Trt, M, Y, tree = tree, 
                         fdr.alpha = 0.1, graph = "circular")
# take a look at phyloseq-class object
demo.physeq <- demo.rsltlst$clean.data
demo.physeq
demo.rsltlst$rslt$PhyloMed.A

