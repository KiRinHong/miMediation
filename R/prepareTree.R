#' Prepare phylogenetic tree to be used in \code{\link{phyloMed}} function
#' 
#' @description The \code{\link{phyloMed}} function requires a rooted and binary tree as input. The \code{prepareTree} 
#' is an utility function that allows users to preprocess their tree if it is unrooted and/or non-binary. The tree 
#' output from \code{prepareTree} function can be directly fed into \code{\link{phyloMed}} function.
#' 
#' @param tree A \code{phylo}-class object.
#' @param verbose  An optional logical value. If \code{TRUE}, tree manipulation process will be printed.
#' Default is \code{FALSE}.
#' @return A binary and rooted tree with nodes and tips conform with following \code{phylo} standards:
#' 
#' \itemize{
#' \item The internal nodes are numbered with value larger than the number of tips. 
#' \item The internal nodes are numbered sequentially, with values increasing away from the root.
#' }
#' 
#' @author Qilin Hong \email{qhong8@@wisc.edu}
#' @examples 
#' # Load real data
#' data(data.cecal)
#' tree = prepareTree(data.cecal$tree)
#' 
#' @importFrom TreeTools Renumber
#' @importFrom ape root multi2di
#' @export

prepareTree <- function(tree, verbose = FALSE){
  if(class(tree) != "phylo") stop("Input tree is not a phylo class!")
  tree$edge = tree$edge[order(tree$edge[,2]),] # order the tips
  if(.is_binary(tree)){
    if(verbose) cat("The phylogeny tree is already binary!\n")
    if(.is_rooted(tree, .ntaxa(tree))){
      tree.new = tree
    }else{
      outgroup = .pick_new_outgroup(tree)
      if(verbose) cat("Root the tree!\n")
      tree.new = root(tree, outgroup = outgroup, resolve.root = TRUE)
    }
  }else{
    if(.is_rooted(tree, .ntaxa(tree))){
      if(verbose) cat("Resolve the multichotomies in the order they appear on the tree!\n")
      tree.new = multi2di(tree, random = FALSE)
    }else{
      outgroup = .pick_new_outgroup(tree)
      if(verbose) cat("Root the tree!\n")
      tree.new = root(tree, outgroup = outgroup, resolve.root = TRUE)
      if(verbose) cat("Resolve the multichotomies in the order they appear on the tree!\n")
      tree.new = multi2di(tree.new, random = FALSE)
    }
  }
  tree.new = Renumber(tree.new) # conform to certain principle
  return(tree.new)
}

# calculate the number of taxa
.ntaxa <- function(tree){
  length(tree$tip.label)
}
# check whether the tree is not
.is_rooted <- function(tree, K){
  if(!is.null(tree$root.edge)) return(TRUE)
  if(tabulate(tree$edge[,1])[K+1]>2) FALSE else TRUE
}
# check whether the tree is binary
.is_binary <- function(tree){
  .ntaxa(tree)-tree$Nnode+.is_rooted(tree,.ntaxa(tree)) == 2
}
.pick_new_outgroup <- function(tree.unrooted){
  treeDT = cbind(cbind(tree.unrooted$edge, tree.unrooted$edge.length)[1:.ntaxa(tree.unrooted),], tree.unrooted$tip.label)
  colnames(treeDT) = c("from", "to", "length", "id")
  # Take the longest terminal branch as outgroup
  new.outgroup = treeDT[which.max(treeDT[, "length"]),"id"]
  return(new.outgroup)
}

