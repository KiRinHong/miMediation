## ---- include = FALSE---------------------------------------------------------
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

