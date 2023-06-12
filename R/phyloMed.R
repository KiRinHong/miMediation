#' Phylogeny-based test of mediation effect in microbiome (PhyloMed)
#' 
#' @description \code{phyloMed} enables us to test the mediation effect in high-dimensional microbial composition.
#' The method leverages the hierarchical phylogeny relationship among different microbial taxa to decompose the
#' complex mediation model on the full microbial composition into multiple simple independent local mediation models
#' on subcompositions. The \code{phyloMed} function (a) performs the mediation test for the subcomposition at each internal node
#' of the phylogenetic tree and pinpoint the mediating nodes with significant test p-values; and (b) combine all subcomposition 
#' p-values to assess the overall mediation effect of the entire microbial community.
#' 
#' @details \code{phyloMed} could leverage phylogeny or taxonomy relationship among taxa. 
#' If the \code{tree} is a unrooted and/or non-binary phylogenetic tree, \code{phyloMed} will preprocess the tree: 
#' (a) root the tree with the longest tip branch as outgroup if it is unrooted; and/or 
#' (b) resolve multichotomies into dichotomies based on the order they appear if it is non-binary. 
#' If the \code{tree} is a taxonomy table, \code{phyloMed} will group taxa based on different levels of taxonomic ranks.
#' \code{phyloMed} uses the treatment-mediator association test p-value and mediation-outcome association test p-value
#' to construct the subcomposition mediation test statistic at each local model (Hong et al., Manuscript). The two p-values can come from
#' either the asymptotic test or the permutation test. Asymptotic test is faster but less accurate when the study sample size is small.
#' By default (\code{n.perm=NULL}), only asymptotic test will be performed. Otherwise, if \code{n.perm} is set to a positive number,
#' results from two versions of PhyloMed will be output, one based on the asymptotic p-value and the other based on the permutation 
#' p-value. Graph only highlights the mediating nodes identified from permutation version when both versions are performed. 
#' By default (\code{graph=NULL}), graph will not be plotted. 
#' 
#' @param treatment A numeric vector of the treatment. 
#' @param mediators A named numeric matrix containing microbiome abundance. Each row is a subject and each column is a taxon. 
#' Column name contains the taxon name. 
#' @param outcome A numeric vector of continuous or binary outcome.
#' @param confounders An optional numeric vector or matrix containing confounders that may affect the 
#' treatment, mediators and outcome. 
#' Each row is a subject and each column is a specific confounder, e.g., age or sex. 
#' Default is \code{NULL}.
#' @param interaction An optional logical value. If \code{TRUE}, the interaction term between treatment and mediator 
#' will be taken into account. 
#' Default is \code{FALSE}.
#' @param tree A phylogenetic tree (\code{phylo-class} object) or a taxonomy table (\code{matrix-class} object). 
#' The tip labels in the phylogenetic tree or the row names in the taxonomy table should overlap with the column names in the \code{mediators} matrix. 
#' The column names in the taxonomy table should start from the higher level to lower level, e.g., from kingdom to genus. See Details.
#' @param pi.method An optional character string denotes the method to used in estimate proportion of null. 
#' Default method is \code{"product"}, an alternative method is \code{"maxp"}. Can be abbreviated.
#' @param fdr.alpha An optional numeric value for the desired FDR significance level in identifying mediating nodes on the tree. 
#' Default is \code{0.05}.
#' @param n.perm An optional numeric value for the maximum number of permutations. 
#' Default is \code{NULL}. See Details.
#' @param verbose An optional logical value. If \code{TRUE}, information of the test on each node will be printed.
#' Default is \code{FALSE}.
#' @param graph An optional character string denotes the layout of the graph, which contains a phylogenetic tree or taxonomic tree with
#' identified mediating nodes highlighted. Can be \code{"circular"} or \code{"rectangular"}. Can be abbreviated.
#' Default is \code{NULL}. See Details. 
#' 
#' @return 
#' A \code{phyloseq}-class object named \code{clean.data} and a list named \code{rslt}.
#' 
#' \code{clean.data} contains the following components:
#' \item{\code{sample_data}}{Input treatment, outcome and confounders.}
#' \item{\code{otu_table}}{The abundance data for the taxa that are present on the tips of the \code{phy_tree} or on the rows of the \code{tax_table}.}
#' \item{\code{tax_table}}{The taxonomy table with rows exactly match the taxa in the \code{otu_table}. \code{NULL} if input \code{tree} is a phylogeny tree.}
#' \item{\code{phy_tree}}{The binary and rooted phylogenetic tree with tips exactly match the taxa in the \code{otu_table}: (a) The internal nodes are numbered 
#' with value larger than the number of tips; (b) The internal nodes are numbered sequentially, with values increasing away from the root. 
#' \code{NULL} if input \code{tree} is a taxonomy table.}
#' 
#' If \code{n.perm} is not \code{NULL}, the function will return two lists in \code{rslt} named \code{PhyloMed.A} and \code{PhyloMed.P}, respectively.
#' Otherwise, only one list named \code{PhyloMed.A} will be returned.
#' 
#' Each list contains the following components:
#' \item{\code{node.pval}}{A numeric vector of subcomposition mediation p-values for all internal nodes.}
#' \item{\code{sig.clade}}{A list of significant nodes with their descendants.}
#' \item{\code{null.prop}}{A vector of the estimated proportion of different types of null hypotheses across all local mediation tests.}
#' \item{\code{global.pval}}{A global test p-value using harmonic mean.}
#' 
#' If \code{graph} is not \code{NULL}, the phylogenetic or taxonomic tree will be plotted. The layout depends on the input of \code{graph}. 
#' The size of the circle at each internal node is proportional to 
#' \eqn{-\log_{10}}(subcomposition p-value), the larger circle indicates a smaller p-value. The significant nodes are highlighted by blue rectangle. 
#' 
#' @author Qilin Hong \email{qhong8@@wisc.edu}
#' @references 
#' Hong, Q., Chen, G., & Tang, Z. Z. (2023). PhyloMed: a phylogeny-based test of mediation effect in microbiome. Genome Biology, 24(1), 1-21.
#' 
#' @keywords PhyloMed
#' @examples
#' # Load real data
#' data(data.cecal)
#' # Run test with phylogeny tree
#' Trt = data.cecal$treatment
#' M = data.cecal$mediators
#' Y = data.cecal$outcome
#' tree = data.cecal$tree
#' rslt.phylomed = phyloMed(Trt, M, Y, tree = tree, graph = "rectangular")
#' # Run test with taxonomy table
#' Trt = data.zeeviD$treatment
#' M = data.zeeviD$mediators
#' Y = data.zeeviD$outcome
#' tree = data.zeeviD$tree
#' rslt.phylomed = phyloMed(Trt, M, Y, tree = tree, graph = "circular")
#' 
#' @importFrom fdrtool gcmlcm
#' @importFrom SKAT SKAT_Null_Model SKAT
#' @importFrom phyloseq otu_table tax_table sample_data sample_names phyloseq
#' @importFrom ape keep.tip root multi2di
#' @importFrom MASS ginv
#' @importFrom TreeTools Renumber
#' @import ggplot2 ggtree harmonicmeanp data.table data.tree
#' @export

phyloMed <- function(treatment, mediators, outcome, confounders = NULL, interaction = FALSE, 
                     tree, pi.method = "product", fdr.alpha = 0.05, 
                     n.perm = NULL, verbose = FALSE, graph = NULL){
  #treatment=Trt; mediators=M; outcome=Y; pi.method = "product"; n.perm = NULL; confounders = NULL;interaction = FALSE; fdr.alpha = 0.1; graph = "circular"; verbose=TRUE
  if(sum(is.na(cbind(treatment, mediators, outcome, confounders)))>0) stop("Input data contain NAs!")
  if(is.data.frame(treatment)) treatment = unlist(treatment)
  if(is.data.frame(mediators)) mediators = as.matrix(mediators)
  if(is.data.frame(outcome)) outcome = unlist(outcome)
  if(any(round(mediators) != mediators))
    stop("Mediators contain relative abundance (proportion), please use the abundance (count) table!")
  
  n.sample = length(treatment)
  if(any(nrow(mediators)!=n.sample, length(outcome)!=n.sample)) stop("Input data must be of same length!")
  
  if(!is.numeric(treatment)) stop("Treatment is not a numeric vector!")
  if(!is.numeric(outcome)) stop("Outcome is not a numeric vector!")
  if(!any(class(tree) %in% "phylo")){
    if(any(is.matrix(tree), is.data.frame(tree))){
      if(is.data.frame(tree)) tax.tab = as.matrix(tree)
      if(is.matrix(tree)) tax.tab = tree
      cat("No phylogenetic tree available, construct taxonomic tree!")
      if(dim(tax.tab)[1] < dim(tax.tab)[2]){
        tax.tab = t(tax.tab)
        cat("Reorgnize the taxonomy table with taxonomic ranks as columns!")
      }
      if(any(is.na(tax.tab))) stop("Taxonomy table contains NAs!")
      empty.sum = apply(tax.tab, 1, function(x) sum(grepl("^\\s*$", x)))
      if(any(empty.sum > 0)) stop("Taxonomy table contains empty strings!")
      if(length(unique(tax.tab[,1])) != 1) stop("The first column of taxonomy table is not unique! Please specify unique Kingdom!")
      tax.tab.unique = unique(tax.tab)
      unique.len = apply(tax.tab.unique, 2, function(x) length(unique(x)))
      if(any(unique.len != sort(unique.len))) stop("Please reorganize taxonomy table from highest rank to lowest rank!")
      for (i in ncol(tax.tab.unique):2) {
        tmp.id = which(table(tax.tab.unique[,i]) > 1)
        if(length(tmp.id) == 0){
          next
        }else{
          tmp.tab = tax.tab.unique[tax.tab.unique[,i] %in% names(table(tax.tab.unique[,i]))[tmp.id],]
          if(nrow(unique(tmp.tab[,1:i])) == length(tmp.id)){
            next
          }else{
            print(unique(tmp.tab[,1:i]))
            stop(sprintf("Taxonomic rank %s contains conflicts, i.e., the higher rank classification is different for the same lower rank classification! Please rename the conflict ones!", 
                         colnames(tax.tab.unique)[i])) 
          }
        }
      }
      input.type = "table"
      tab = tax.tab
      ori.tab.colnames = colnames(tab)
      colnames(tab) = paste0("Rank", 1:ncol(tab))
      if(all(colnames(mediators) %in% rownames(tab))){
        if(ncol(mediators) == nrow(tab)){
          mediators = mediators[,rownames(tab)]
        }else{
          cat("Prune the taxonomic table based on the column names of mediators!")
          tab = tab[colnames(mediators),]
          mediators = mediators[,rownames(tab)]
        }
      }else if(all(rownames(tab) %in% colnames(mediators))){
        cat("Subset the mediators based on the taxonomy table!")
        mediators = mediators[,rownames(tab)]
      }else if(length(intersect(rownames(tab), colnames(mediators))) > 0){
        taxa.cross = intersect(rownames(tab), colnames(mediators))
        cat(sprintf("Prune the taxonomy table based on %g overlapped taxa!", length(taxa.cross)))
        tab = tab[taxa.cross,]
        cat(sprintf("Subset the mediators based on %g overlapped taxa!", length(taxa.cross)))
        mediators = mediators[,rownames(tab)]
      }else{
        stop("The column names of mediators do not match the row names of taxonomic table!")
      }
    }else{
      stop("Input data should be either a phylogeny tree or taxonomy table!")
    }
  }else{
    input.type = "tree"
    cat("Run phyloMed based on phylogenetic tree!")
    tree = .prepareTree(tree, verbose = verbose)
    if(all(colnames(mediators) %in% tree$tip.label)){
      if(ncol(mediators) == .ntaxa(tree)){
        mediators = mediators[,tree$tip.label]
      }else{
        cat("Prune the phylogenetic tree based on the column names of mediators!")
        tree.trim = keep.tip(tree, colnames(mediators))
        tree = .prepareTree(tree.trim, verbose = verbose)
        mediators = mediators[,tree$tip.label]
      }
    }else if(all(tree$tip.label %in% colnames(mediators))){
      cat("Subset the mediators based on the tip labels on the tree!")
      mediators = mediators[,tree$tip.label]
    }else if(length(intersect(tree$tip.label, colnames(mediators))) > 0){
      taxa.cross = intersect(tree$tip.label, colnames(mediators))
      cat(sprintf("Prune the phylogenetic tree based on %g overlapped taxa!", length(taxa.cross)))
      tree.trim = keep.tip(tree, taxa.cross)
      tree = .prepareTree(tree.trim, verbose = verbose)
      cat(sprintf("Subset the mediators based on %g overlapped taxa!", length(taxa.cross)))
      mediators = mediators[,tree$tip.label]
    }else{
      stop("The column names of mediators do not match the tip labels!")
    }
  }
  
  method = match.arg(tolower(pi.method), choices = c("product", "maxp"))
  if(!is.null(graph)) layout.method = match.arg(tolower(graph), choices = c("circular", "rectangular"))
  if(verbose) cat(sprintf("Use %s's method to obtain probability of null hypotheses estimates\n", toupper(method)))
  
  if(is.null(confounders)){ 
    conf.ori = matrix(1, nrow = n.sample, ncol = 1)
  }else if(all(is.vector(confounders), unique(confounders) == 1)){
    conf.ori = matrix(confounders, nrow = n.sample, ncol = 1) # handle the case if user include the intercept into covariate
  }else if(all(!is.vector(confounders), unique(confounders[,1]) == 1)){
    conf.ori = confounders # handle the case if user include the intercept into covariate
  }else{
    conf.ori = cbind(1, confounders)
  }
  
  Trt.ori = treatment
  outcome.ori = outcome
  M = mediators
  if(input.type == "table"){
    K = 0
    n.rank = ncol(tab)
    n.level = n.rank - 1
    for(k in 1:n.level){
      K = K + sum(table(tab[,n.rank-k]) > 1)
    }
    #### here, visit every internal node of the tree and generate p-values for testing alpha and beta respectively
    B.max = n.perm
    if(!is.null(n.perm)){
      R.sel = .choose_r(fdr.alpha/K, 0.05)
      if(R.sel > B.max) cat(sprintf("The maximal permutation times (%g) is smaller than the chosen maximal number of successes (%g). 
    The permutation p-value may not be precise enough. Please increase the maxmimal permutation times!", B.max, R.sel))
      if(verbose){
        cat(sprintf("Adapative permutation procedure perfroms %g times at most\n", B.max))
        cat(sprintf("Adapative permutation procedure requires %g exceeds\n", R.sel))
      }
    }
    W.data = data.table(data.frame(tab, t(M)))
    n.rank = ncol(tab)
    otucols = names(W.data)[-(1:n.rank)]
    n.level = n.rank-1
    
    subtree = chi.stat.alpha = z.stat.alpha = pval.alpha.asym = pval.alpha.perm = 
      chi.stat.beta = z.stat.beta = pval.beta.asym = pval.beta.perm = NULL

    for(k in 1:n.level){
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")
      
      tmp = table(tab[,n.rank-k])
      level.uni = sort(names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)    
      
      tt = W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
      W.count = data.matrix(tt[, otucols, with=FALSE])      
      
      current.parent = ori.tab.colnames[n.rank-k]
      for(m in 1:m.level){
        current.child = level.uni[m]
        current.node = paste0(current.parent, ".", current.child)
        if(verbose) cat(sprintf("=====Processing internal node #%s=====\n", current.node))
        Mc = t(W.count[which(W.tax == level.uni[m]), , drop=FALSE])
        remove.index = which(colSums(Mc) == 0)
        if(length(remove.index) == ncol(Mc)){
          if(verbose) cat("Some children have all observations being 0, skip node\n")
          next
        }else{
          if(length(remove.index)>0){
            Mc = Mc[, -remove.index, drop=FALSE] 
          }
          if(ncol(Mc) == 1){
            if(verbose) cat("Only exists one child, skip node\n")
            next
          }else{
            Mc2 = Mc + 0.5
            idx = which(rowSums(Mc) != 0)
            ref.id = which.max(colMeans(Mc2/rowSums(Mc2)))
            if(length(idx) != 0){
              Trt = Trt.ori[idx]; conf = conf.ori[idx,,drop=FALSE]; outcome = outcome.ori[idx]
              G = log(Mc2[idx,-ref.id]/as.numeric(Mc2[idx,ref.id]))
            }else{
              Trt = Trt.ori; conf = conf.ori; outcome = outcome.ori
              G = log(Mc2[,-ref.id]/as.numeric(Mc2[,ref.id]))
            }
            
            if(interaction){
              G2 = cbind(G, Trt*G)
            }
            
            # rank < 3
            condition1 = FALSE
            if(interaction){
              condition1 = (qr(cbind(Trt,G2))$rank < (1+ncol(G2)))
              if(all(condition1, verbose)) cat("Matrix (T, G, Trt*G) is not full rank, skip node\n")
            }
            
            # remove the subjects with all subcomponents being zero may result in collinearity
            # the outcome may be unique after removing the subjects when is binary
            condition2 = (qr(cbind(conf,Trt))$rank < (ncol(conf)+1)) || (length(unique(outcome)) == 1)
            if(all(condition2, verbose)) cat("Matrix (Covariates, Trt) is not full rank or outcome is unique, skip node\n")
            
            # only have two effective observations
            condition3 = (length(outcome) == 2)
            if(all(condition3, verbose)) cat("The effective sample size equals to 2, skip node\n")
            
            if(any(condition1,condition2,condition3)){
              chi.stat.alpha = c(chi.stat.alpha, NA)
              z.stat.alpha = c(z.stat.alpha, NA)
              pval.alpha.asym = c(pval.alpha.asym, NA)
              pval.alpha.perm = c(pval.alpha.perm, NA)
              chi.stat.beta = c(chi.stat.beta, NA)
              z.stat.beta = c(z.stat.beta, NA)
              pval.beta.asym = c(pval.beta.asym, NA)
              pval.beta.perm = c(pval.beta.perm, NA)
              next
            }
            subtree = c(subtree, current.node)
            
            Rconf = diag(nrow(conf)) - conf %*% solve(t(conf) %*% conf) %*% t(conf)
            
            if(is.matrix(G)){
              mod.full = lm(G ~ 0+conf+Trt)
              est = mod.full$coef["Trt",]
              TestAlpha = .test_alpha_vc(G, Trt, conf)
              stat = TestAlpha$stat
              pval = TestAlpha$pval
              
              chi.stat.alpha = c(chi.stat.alpha, stat)
              z.stat.alpha = c(z.stat.alpha, sqrt(stat)*sign(sum(est)))
              pval.alpha.asym = c(pval.alpha.asym, pval)
            }else{
              mod.full = summary(glm(G~0+conf+Trt,family = quasi()))
              est = mod.full$coefficients["Trt","Estimate"]
              mod.reduce = summary(glm(G~0+conf,family = quasi()))
              mod.reduce.disp = mod.reduce$dispersion
              mod.reduce.resid = mod.reduce$deviance.resid
              TestAlpha = .test_alpha(G, Trt, conf, mod.reduce.resid, mod.reduce.disp)
              stat = TestAlpha$stat
              pval = TestAlpha$pval
              
              chi.stat.alpha = c(chi.stat.alpha, stat)
              z.stat.alpha = c(z.stat.alpha, sqrt(stat)*sign(est))
              pval.alpha.asym = c(pval.alpha.asym, pval)
            }
            
            if(length(unique(outcome)) > 2) {
              # continuous traits
              if(interaction){
                obj = SKAT_Null_Model(outcome~0+conf+Trt, out_type = "C")
                TestBeta = .test_beta(outcome, G2, Trt, conf, obj = obj, test.type = "vc")
              }else{
                if(is.matrix(G)){
                  obj = SKAT_Null_Model(outcome~0+conf+Trt, out_type = "C")
                  TestBeta = .test_beta_vc(outcome, G, Trt, conf, obj = obj) 
                }else{
                  mod = summary(lm(outcome~cbind(conf[,-1], Trt)))
                  mod.s2 = mod$sigma^2
                  mod.resid = mod$residuals
                  TestBeta = .test_beta(outcome, G, Trt, conf, resid.obs = mod.resid, s2.obs = mod.s2, test.type = "mv")
                }
              }
            }else{
              # binary traits
              if(interaction){
                obj = SKAT_Null_Model(outcome~0+conf+Trt, out_type = "D")
                TestBeta = .test_beta(outcome, G2, Trt, conf, obj = obj, test.type = "vc")
              }else{
                if(is.matrix(G)){
                  obj = SKAT_Null_Model(outcome~0+conf+Trt, out_type = "D")
                  TestBeta = .test_beta_vc(outcome, G, Trt, conf, obj = obj)
                }else{
                  mod = glm(outcome~cbind(conf[,-1], Trt), family = "binomial")
                  mod.est = mod$coefficients
                  TestBeta = .test_beta(outcome, G, Trt, conf, est.obs = mod.est, test.type = "mv")
                }
              }
            }
            
            chi.stat.beta = c(chi.stat.beta, TestBeta$stat)
            z.stat.beta = c(z.stat.beta, sqrt(TestBeta$stat)*sign(sum(TestBeta$est)))
            pval.beta.asym = c(pval.beta.asym, TestBeta$pval)
            
            if(!is.null(n.perm)){
              #### permutation for alpha
              if(is.infinite(stat)){
                pval.alpha.perm = c(pval.alpha.perm, 1/(B.max+1))
                if(verbose){
                  cat(sprintf("# of permutations for alpha: %g, Test statistic = Inf\n", B.max))
                  cat("Use Pseudo-ECDF approximation p-value\n")
                }
              }else{
                m = 1
                Nexc = 0
                alpha.stat.perm = numeric(B.max)
                while (Nexc < R.sel & m < B.max) {
                  TestAlpha_perm = .test_alpha_vc(G, Rconf[sample(nrow(Rconf)),] %*% Trt, conf)
                  stat_perm = TestAlpha_perm$stat
                  if(stat_perm >= stat) Nexc = Nexc + 1
                  alpha.stat.perm[m] = stat_perm
                  m = m + 1
                }
                if(m < B.max){
                  pval.alpha.perm = c(pval.alpha.perm, Nexc/(m-1))
                  if(verbose){
                    cat(sprintf("# of permutations for alpha: %g\n", m-1))
                    cat("Use ECDF approximation p-value\n")
                  }
                }else{
                  if(Nexc <= 10){
                    tmp.alpha.perm = tryCatch(.gpd_approx(alpha.stat.perm, 250, stat), error=function(err) NA)
                    pval.alpha.perm = c(pval.alpha.perm, tmp.alpha.perm)
                    if(verbose){
                      cat(sprintf("# of permutations for alpha: %g\n", B.max))
                      cat("Use GPD approximation p-value\n")
                    }
                    if(is.na(tmp.alpha.perm)){
                      pval.alpha.perm = c(pval.alpha.perm, (Nexc+1)/(B.max+1))
                      if(verbose) cat("Fail to fit GPD, use Pseudo-ECDF approximation p-value\n")
                    }
                  }else{
                    pval.alpha.perm = c(pval.alpha.perm, Nexc/B.max)
                    if(verbose){
                      cat(sprintf("# of permutations for alpha: %g\n", B.max))
                      cat("Use ECDF approximation p-value\n")
                    }
                  }
                }
              }
              
              #### permutation for beta
              if(is.infinite(TestBeta$stat)){
                pval.beta.perm = c(pval.beta.perm, 1/(B.max+1))
                if(verbose){
                  cat(sprintf("# of permutations for beta: %g, Test statistic = Inf\n", B.max))
                  cat("Use Pseudo-ECDF approximation p-value\n")
                }
              }else{
                m = 1
                Nexc = 0
                beta.stat.perm = numeric(B.max)
                while (Nexc < R.sel & m < B.max) {
                  if(interaction){
                    tmp_beta = .test_beta(outcome, Rconf[sample(nrow(Rconf)),] %*% G2, Trt, conf, obj = obj, test.type = "vc") 
                  }else{
                    if(length(unique(outcome)) > 2){
                      if(is.matrix(G)){
                        tmp_beta = .test_beta_vc(outcome, Rconf[sample(nrow(Rconf)),] %*% G, Trt, conf, obj = obj) 
                      }else{
                        tmp_beta = .test_beta(outcome, Rconf[sample(nrow(Rconf)),] %*% G, Trt, conf, resid.obs = mod.resid, s2.obs = mod.s2, test.type = "mv")
                      }
                    }else{
                      if(is.matrix(G)){
                        tmp_beta = .test_beta_vc(outcome, Rconf[sample(nrow(Rconf)),] %*% G, Trt, conf, obj = obj)
                      }else{
                        tmp_beta = .test_beta(outcome, Rconf[sample(nrow(Rconf)),] %*% G, Trt, conf, est.obs = mod.est, test.type = "mv")
                      }
                    }
                  }
                  if(tmp_beta$stat >= TestBeta$stat) Nexc = Nexc + 1
                  beta.stat.perm[m] = tmp_beta$stat
                  m = m + 1
                }
                if(m < B.max){
                  pval.beta.perm = c(pval.beta.perm, Nexc/(m-1))
                  if(verbose){
                    cat(sprintf("# of permutations for beta: %g\n", m))
                    cat("Use ECDF approximation p-value\n")
                  }
                }else{
                  if(Nexc <= 10){
                    tmp.beta.perm = tryCatch(.gpd_approx(beta.stat.perm, 250, TestBeta$stat), error=function(err) NA)
                    pval.beta.perm = c(pval.beta.perm, tmp.beta.perm)
                    if(verbose){
                      cat(sprintf("# of permutations for beta: %g\n", B.max))
                      cat("Use GPD approximation p-value\n")
                    }
                    if(is.na(tmp.beta.perm)){
                      pval.beta.perm = c(pval.beta.perm, (Nexc+1)/(B.max+1))
                      if(verbose) cat("Fail to fit GPD, use Pseudo-ECDF approximation p-value\n")
                    }
                  }else{
                    pval.beta.perm = c(pval.beta.perm, Nexc/B.max)
                    if(verbose){
                      cat(sprintf("# of permutations for beta: %g\n", B.max))
                      cat("Use ECDF approximation p-value\n")
                    }
                  }
                }
              }
            }
          }
        }
      }
    }# lineage loop
    
  }
  
  if(input.type == "tree"){
    K = .ntaxa(tree)
    #### here, visit every internal node of the tree and generate p-values for testing alpha and beta respectively
    B.max = n.perm
    if(!is.null(n.perm)){
      R.sel = .choose_r(fdr.alpha/K, 0.05)
      if(R.sel > B.max) cat(sprintf("The maximal permutation times (%g) is smaller than the chosen maximal number of successes (%g). 
    The permutation p-value may not be precise enough. Please increase the maxmimal permutation times!", B.max, R.sel))
      if(verbose){
        cat(sprintf("Adapative permutation procedure perfroms %g times at most\n", B.max))
        cat(sprintf("Adapative permutation procedure requires %g exceeds\n", R.sel))
      }
    }
    treestructure = .phylostructure(tree)
    ## here perform tree-based method
    chi.stat.alpha = chi.stat.beta = z.stat.alpha = z.stat.beta = pval.alpha.asym = pval.beta.asym = pval.alpha.perm = pval.beta.perm = numeric(tree$Nnode)
    for (i in (K + 1):(K + tree$Nnode)) {
      if(verbose) cat(sprintf("=====Processing internal node #%g=====\n", i))
      child.left = treestructure$phylochildren[i,1]
      child.right = treestructure$phylochildren[i,2]
      
      Mc.left = rowSums(M[,treestructure$descendant[child.left,], drop = FALSE])
      Mc.right = rowSums(M[,treestructure$descendant[child.right,], drop = FALSE])
      
      if(mean(Mc.left) < mean(Mc.right)){
        Mc = cbind(Mc.right, Mc.left)
      }else{
        Mc = cbind(Mc.left, Mc.right)
      }
      
      Mc2 = Mc + 0.5
      idx = which(rowSums(Mc) != 0)
      if(length(idx) != 0){
        Trt = Trt.ori[idx]; conf = conf.ori[idx,,drop=FALSE]; outcome = outcome.ori[idx]
        G = log(Mc2[idx,-2]/as.numeric(Mc2[idx,2]))
      }else{
        Trt = Trt.ori; conf = conf.ori; outcome = outcome.ori
        G = log(Mc2[,-2]/as.numeric(Mc2[,2]))
      }
      
      if(interaction){
        G2 = cbind(G, Trt*G)
      }
      
      # For some internal node, one child have all obs being 0. We need to skip such internal node, and set the results to NA
      condition1 = any(colSums(Mc) == 0)
      if(all(condition1,verbose)) cat(sprintf("Some children have all observations being 0, skip internal node #%g\n", i))
      
      # rank < 3
      condition2 = FALSE
      if(interaction){
        condition2 = (qr(cbind(Trt,G2))$rank < 3)
        if(all(condition2, verbose)) cat(sprintf("Matrix (T, G, Trt*G) is not full rank, skip internal node #%g\n", i))
      }
      
      # remove the subjects with all subcomponents being zero may result in collinearity
      # the outcome may be unique after removing the subjects when is binary
      condition3 = (qr(cbind(conf,Trt))$rank < (ncol(conf)+1)) || (length(unique(outcome)) == 1)
      if(all(condition3, verbose)) cat(sprintf("Matrix (Covariates, Trt) is not full rank or outcome is unique, skip internal node #%g\n", i))
      
      # only have two effective observations
      condition4 = (length(outcome) == 2)
      if(all(condition4, verbose)) cat(sprintf("The effective sample size equals to 2, skip internal node #%g\n", i))
      
      if(any(condition1,condition2,condition3,condition4)){
        chi.stat.alpha[i-K] = NA
        z.stat.alpha[i-K] = NA
        pval.alpha.asym[i-K] = NA
        pval.alpha.perm[i-K] = NA
        chi.stat.beta[i-K] = NA
        z.stat.beta[i-K] = NA
        pval.beta.asym[i-K] = NA
        pval.beta.perm[i-K] = NA
        next
      }
      
      # compute residual forming matrix Rconf (Smith's method)
      Rconf = diag(nrow(conf)) - conf %*% solve(t(conf) %*% conf) %*% t(conf)
      
      mod.full = summary(glm(G~0+conf+Trt,family = quasi()))
      est = mod.full$coefficients["Trt","Estimate"]
      mod.reduce = summary(glm(G~0+conf,family = quasi()))
      mod.reduce.disp = mod.reduce$dispersion
      mod.reduce.resid = mod.reduce$deviance.resid
      TestAlpha = .test_alpha(G, Trt, conf, mod.reduce.resid, mod.reduce.disp)
      stat = TestAlpha$stat
      pval = TestAlpha$pval
      
      chi.stat.alpha[i-K] = stat
      z.stat.alpha[i-K] = sqrt(stat)*sign(est)
      pval.alpha.asym[i-K] = pval
      
      if(length(unique(outcome)) > 2) {
        # continuous traits
        if(interaction){
          obj = SKAT_Null_Model(outcome~0+conf+Trt, out_type = "C")
          TestBeta = .test_beta(outcome, G2, Trt, conf, obj = obj, test.type = "vc") # est[1] ~ mediator est[2] ~ exposure * mediator
        }else{
          mod = summary(lm(outcome~cbind(conf[,-1], Trt)))
          mod.s2 = mod$sigma^2
          mod.resid = mod$residuals
          TestBeta = .test_beta(outcome, G, Trt, conf, resid.obs = mod.resid, s2.obs = mod.s2, test.type = "mv")
        }
      } else {
        # binary traits
        if(interaction){
          obj = SKAT_Null_Model(outcome~0+conf+Trt, out_type = "D")
          TestBeta = .test_beta(outcome, G2, Trt, conf, obj = obj, test.type = "vc") # est[1] ~ mediator est[2] ~ exposure * mediator
        }else{
          mod = glm(outcome~cbind(conf[,-1], Trt), family = "binomial")
          mod.est = mod$coefficients
          TestBeta = .test_beta(outcome, G, Trt, conf, est.obs = mod.est, test.type = "mv")
        }
      }
      
      chi.stat.beta[i-K] = TestBeta$stat
      z.stat.beta[i-K] = sqrt(TestBeta$stat)*sign(sum(TestBeta$est))
      pval.beta.asym[i-K] = TestBeta$pval
      
      if(!is.null(n.perm)){
        #### permutation for alpha
        if(is.infinite(stat)){
          pval.alpha.perm[i-K] = 1/(B.max+1)
          if(verbose){
            cat(sprintf("# of permutations for alpha: %g, Test statistic = Inf\n", B.max))
            cat("Use Pseudo-ECDF approximation p-value\n")
          }
        }else{
          m = 1
          Nexc = 0
          alpha.stat.perm = numeric(B.max)
          while (Nexc < R.sel & m < B.max) {
            TestAlpha_perm = .test_alpha(G, Rconf[sample(nrow(Rconf)),] %*% Trt, conf, mod.reduce.resid, mod.reduce.disp)
            stat_perm = TestAlpha_perm$stat
            if(stat_perm >= stat) Nexc = Nexc + 1
            alpha.stat.perm[m] = stat_perm
            m = m + 1
          }
          if(m < B.max){
            pval.alpha.perm[i-K] = Nexc/(m-1)
            if(verbose){
              cat(sprintf("# of permutations for alpha: %g\n", m-1))
              cat("Use ECDF approximation p-value\n")
            }
          }else{
            if(Nexc <= 10){
              pval.alpha.perm[i-K] = tryCatch(.gpd_approx(alpha.stat.perm, 250, stat), error=function(err) NA)
              if(verbose){
                cat(sprintf("# of permutations for alpha: %g\n", B.max))
                cat("Use GPD approximation p-value\n")
              }
              if(is.na(pval.alpha.perm[i-K])){
                pval.alpha.perm[i-K] = (Nexc+1)/(B.max+1)
                if(verbose) cat("Fail to fit GPD, use Pseudo-ECDF approximation p-value\n")
              }
            }else{
              pval.alpha.perm[i-K] = Nexc/B.max
              if(verbose){
                cat(sprintf("# of permutations for alpha: %g\n", B.max))
                cat("Use ECDF approximation p-value\n")
              }
            }
          }
        }
        
        #### permutation for beta
        if(is.infinite(TestBeta$stat)){
          pval.beta.perm[i-K] = 1/(B.max+1)
          if(verbose){
            cat(sprintf("# of permutations for beta: %g, Test statistic = Inf\n", B.max))
            cat("Use Pseudo-ECDF approximation p-value\n")
          }
        }else{
          m = 1
          Nexc = 0
          beta.stat.perm = numeric(B.max)
          while (Nexc < R.sel & m < B.max) {
            if(interaction){
              # est[1] ~ mediator est[2] ~ exposure * mediator
              tmp_beta = .test_beta(outcome, Rconf[sample(nrow(Rconf)),] %*% G2, Trt, conf, obj = obj, test.type = "vc") 
            }else{
              if(length(unique(outcome)) > 2){
                tmp_beta = .test_beta(outcome, Rconf[sample(nrow(Rconf)),] %*% G, Trt, conf, resid.obs = mod.resid, s2.obs = mod.s2, test.type = "mv")
              }else{
                tmp_beta = .test_beta(outcome, Rconf[sample(nrow(Rconf)),] %*% G, Trt, conf, est.obs = mod.est, test.type = "mv")
              }
            }
            if(tmp_beta$stat >= TestBeta$stat) Nexc = Nexc + 1
            beta.stat.perm[m] = tmp_beta$stat
            m = m + 1
          }
          if(m < B.max){
            pval.beta.perm[i-K] = Nexc/(m-1)
            if(verbose){
              cat(sprintf("# of permutations for beta: %g\n", m))
              cat("Use ECDF approximation p-value\n")
            }
          }else{
            if(Nexc <= 10){
              pval.beta.perm[i-K] = tryCatch(.gpd_approx(beta.stat.perm, 250, TestBeta$stat), error=function(err) NA)
              if(verbose){
                cat(sprintf("# of permutations for beta: %g\n", B.max))
                cat("Use GPD approximation p-value\n")
              }
              if(is.na(pval.beta.perm[i-K])){
                pval.beta.perm[i-K] = (Nexc+1)/(B.max+1)
                if(verbose) cat("Fail to fit GPD, use Pseudo-ECDF approximation p-value\n")
              }
            }else{
              pval.beta.perm[i-K] = Nexc/B.max
              if(verbose){
                cat(sprintf("# of permutations for beta: %g\n", B.max))
                cat("Use ECDF approximation p-value\n")
              }
            }
          }
        }
        
      }
    }
  }
  
  
  if(method == "product"){
    alpha1.asym = .pi0_JC(na.omit(z.stat.alpha))
    alpha2.asym = .pi0_JC(na.omit(z.stat.beta))
    tmp.asym = .nullEstimation_prod(pval.alpha.asym, pval.beta.asym, alpha1.asym, alpha2.asym)
  }else if(method == "maxp"){
    tmp.asym = .nullEstimation_minus(pval.alpha.asym, pval.beta.asym, alpha1.asym, alpha2.asym, z.stat.alpha, z.stat.beta)
  }

  rawp.asym = tmp.asym$rawp
  if(length(which(rawp.asym==0))>0)  rawp.asym[rawp.asym == 0] = runif(length(which(rawp.asym==0)), min = 0, max = 1e-7)
  null.prop.est.asym = c(tmp.asym$alpha00, tmp.asym$alpha10, tmp.asym$alpha01)
  names(null.prop.est.asym) = c("H00","H10","H01")
  
  p.asym.adj = p.adjust(rawp.asym, method = "BH")
  sig.nodeID.asym = which(p.asym.adj < fdr.alpha)
  
  if(input.type == "table"){
    if(length(sig.nodeID.asym) > 0){
      sig.clade.asym = lapply(subtree[sig.nodeID.asym], 
                              function(x){
                                tmp = unlist(strsplit(x, split = "\\.", ))
                                return(as.vector(tab[which(tab[,which(ori.tab.colnames == tmp[1])] == tmp[2]), ncol(tab)]))
                              })
      names(sig.clade.asym) = subtree[sig.nodeID.asym]
    }else{
      sig.clade.asym = NULL
    }
    names(rawp.asym) = subtree
    df = as.data.frame(tab)
    for (i in 1:ncol(df)) {
      df[,i] = paste0(ori.tab.colnames[i],".",df[,i])
    }
    df$pathString = do.call(paste, c(df, sep = "/"))
    tax.tree = as.Node(df, pathDelimiter = "/")
    tree.vis = as.phylo.Node(tax.tree)
    rawp = rep(NA, length(tree.vis$node.label))
    names(rawp) = tree.vis$node.label
    node.name = subtree
    # all(node.name %in% tree.vis$node.label)
    if(all(!is.null(graph), is.null(n.perm))){
      rawp[node.name] = rawp.asym
      tree.vis$node.label = rawp
      tree.vis$tip.label = gsub(paste0("^",ori.tab.colnames[length(ori.tab.colnames)],"\\."), "", tree.vis$tip.label)
      
      g.asym = ggtree(tree.vis, layout=layout.method) +
        geom_point2(aes(subset=!isTip), shape=21, size=-log10(as.numeric(rawp))*3, fill = "red") +
        geom_tiplab(size=2) + 
        theme_tree(plot.margin=margin(5,5,5,5))
      if(length(sig.nodeID.asym) > 0) {
        nodelab = tree.vis$node.label[names(sig.clade.asym)]
        nodeids = nodeid(tree.vis, nodelab)
        cladedat = data.frame(id=nodeids, class=names(sig.clade.asym))
        g.asym = g.asym +
          geom_hilight(data=cladedat, mapping=aes(node=id), alpha=.6, fill="steelblue") +
          geom_cladelab(data=cladedat, mapping=aes(node=id, label=class), vjust=-.5, hjust=-.5, fontsize=3, fontface=4)
      }
      print(g.asym)
    }
  }
  
  if(input.type == "tree"){
    taxa.names = tree$tip.label
    if(length(sig.nodeID.asym) > 0){
      sig.clade.asym = lapply(sig.nodeID.asym+K, function(x) taxa.names[which(treestructure$descendant[x,])])
      names(sig.clade.asym) = paste0("Node", sig.nodeID.asym+K)
    }else{
      sig.clade.asym = NULL
    }
    names(rawp.asym) = paste0("Node", K + 1:length(rawp.asym))
    if(all(!is.null(graph), is.null(n.perm))){
      tree.vis = tree
      tree.vis$node.label = rawp.asym
      g.asym = ggtree(tree.vis, layout = layout.method, branch.length = "none") + 
        geom_point2(aes(subset=!isTip), shape=21, size=-log10(as.numeric(rawp.asym))*3, fill = "red") +
        geom_tiplab(size=2) + 
        theme_tree(plot.margin=margin(5,5,5,5))
      
      if(length(sig.nodeID.asym) > 0) {
        labels = rep(NA, 2*K-1); labels[sig.nodeID.asym+K] = names(sig.clade.asym)
        g.asym = g.asym + geom_text2(aes(subset=!isTip, label=labels), vjust=-.5, hjust=-.5, angle = 0, size=3)
        for (i in 1:length(sig.nodeID.asym)) {
          g.asym = g.asym + geom_hilight(node = sig.nodeID.asym[i]+K, fill = "steelblue", alpha = .6)
        }
      }
      print(g.asym)
    }
  }
  
  rawp.asym.rm = na.omit(rawp.asym)
  L = length(rawp.asym.rm)
  globalp.asym = p.hmp(rawp.asym.rm, w = rep(1/L, L), L = L)
  names(globalp.asym) = "HMP"
  rslt = list(PhyloMed.A = list(node.pval = rawp.asym, sig.clade = sig.clade.asym, 
                                null.prop = null.prop.est.asym, global.pval = globalp.asym))
  
  if(!is.null(n.perm)){
    if(method == "product"){
      alpha1.perm = .pi0_JC(na.omit(abs(qnorm(pval.alpha.perm/2, lower.tail = FALSE))*sign(z.stat.alpha)))
      alpha2.perm = .pi0_JC(na.omit(abs(qnorm(pval.beta.perm/2, lower.tail = FALSE))*sign(z.stat.beta)))
      tmp.perm = .nullEstimation_prod(pval.alpha.perm, pval.beta.perm, alpha1.perm, alpha2.perm)
    }else if(method == "maxp"){
      tmp.perm = .nullEstimation_minus(pval.alpha.perm, pval.beta.perm, alpha1.perm, alpha2.perm,
                                       abs(qnorm(pval.alpha.perm/2, lower.tail = FALSE))*sign(z.stat.alpha),
                                       abs(qnorm(pval.beta.perm/2, lower.tail = FALSE))*sign(z.stat.beta))
    }

    rawp.perm = tmp.perm$rawp
    if(length(which(rawp.perm==0))>0)  rawp.perm[rawp.perm == 0] = runif(length(which(rawp.perm==0)), min = 0, max = 1e-7)
    null.prop.est.perm = c(tmp.perm$alpha00, tmp.perm$alpha10, tmp.perm$alpha01)
    names(null.prop.est.perm) = c("H00","H10","H01")
    p.perm.adj = p.adjust(rawp.perm, method = "BH")
    sig.nodeID.perm = which(p.perm.adj < fdr.alpha)
    
    if(input.type == "table"){
      if(length(sig.nodeID.perm) > 0){
        sig.clade.perm = lapply(subtree[sig.nodeID.perm], 
                                function(x){
                                  tmp = unlist(strsplit(x, split = "\\.", ))
                                  return(as.vector(tab[which(tab[,which(ori.tab.colnames == tmp[1])] == tmp[2]), ncol(tab)]))
                                })
        names(sig.clade.perm) = subtree[sig.nodeID.perm]
      }else{
        sig.clade.perm = NULL
      }
      names(rawp.perm) = subtree
      if(!is.null(graph)){
        rawp[node.name] = rawp.perm
        tree.vis$node.label = rawp
        
        g.perm = ggtree(tree.vis, layout=layout.method) +
          geom_point2(aes(subset=!isTip), shape=21, size=-log10(as.numeric(rawp))*3, fill = "red") +
          geom_tiplab(size=2) + 
          theme_tree(plot.margin=margin(5,5,5,5))
        if(length(sig.nodeID.perm) > 0) {
          nodelab = tree.vis$node.label[names(sig.clade.perm)]
          nodeids = nodeid(tree.vis, nodelab)
          cladedat = data.frame(id=nodeids, class=names(sig.nodeID.perm))
          g.perm = g.perm +
            geom_hilight(data=cladedat, mapping=aes(node=id), alpha=.6, fill="steelblue") +
            geom_cladelab(data=cladedat, mapping=aes(node=id, label=class), vjust=-.2, hjust=-.2, fontsize=3, fontface=4)
        }
        print(g.perm)
      }
    }
    
    if(input.type == "tree"){
      if(length(sig.nodeID.perm) > 0){
        sig.clade.perm = lapply(sig.nodeID.perm+K, function(x) taxa.names[which(treestructure$descendant[x,])])
        names(sig.clade.perm) = paste0("Node", sig.nodeID.perm+K)
      }else{
        sig.clade.perm = NULL
      }
      names(rawp.perm) = paste0("Node", K + 1:length(rawp.perm))
      
      if(!is.null(graph)){
        tree.vis = tree
        tree.vis$node.label = rawp.perm
        g.perm = ggtree(tree.vis, layout = layout.method, branch.length = "none") + 
          geom_point2(aes(subset=!isTip), shape=21, size=-log10(as.numeric(rawp.asym))*3, fill = "red") +
          geom_tiplab(size=3) + 
          theme_tree(plot.margin=margin(5,5,5,5))
        if(length(sig.nodeID.perm) > 0) {
          labels = rep(NA, 2*K-1); labels[sig.nodeID.perm+K] = names(sig.clade.perm)
          g.perm = g.perm + geom_text2(aes(subset=!isTip, label=labels), vjust=-.5, hjust=-.5, angle = 0, size=3)
          for (i in 1:length(sig.nodeID.perm)) {
            g.perm = g.perm + geom_hilight(node = sig.nodeID.perm[i]+K, fill = "steelblue", alpha = .6)
          }
        }
        print(g.perm)
      }
    }
    
    rawp.perm.rm = na.omit(rawp.perm)
    globalp.perm = p.hmp(rawp.perm.rm, w = rep(1/L, L), L = L)
    names(globalp.perm) = "HMP"
    rslt = list(PhyloMed.A = list(node.pval = rawp.asym, sig.clade = sig.clade.asym, 
                                  null.prop = null.prop.est.asym, global.pval = globalp.asym),
                PhyloMed.P = list(node.pval = rawp.perm, sig.clade = sig.clade.perm, 
                                  null.prop = null.prop.est.perm, global.pval = globalp.perm))
  }
  
  OTU = otu_table(M, taxa_are_rows = FALSE)
  if(input.type == "table"){
    TAX = tax_table(tab)
    physeq.tree = NULL
  }
  if(input.type == "tree"){
    TAX = NULL
    physeq.tree = tree
  }

  meta.df = as.data.frame(cbind(Trt.ori, outcome.ori, conf.ori[,-1]))
  names(meta.df)[1:2] = c("treatment", "outcome")
  rownames(meta.df) = rownames(OTU)
  meta.data = sample_data(meta.df)
  # sample_names(meta.data) = rownames(OTU)
  physeq = phyloseq(OTU, TAX, meta.data, physeq.tree)
  rslt.comb = list(clean.data = physeq, rslt = rslt)
  return(rslt.comb)
}


#### INTERNAL FUNCTIONS ####

# generate the tree structure, adapted from PhyloScan
.phylostructure <- function (tree) {
  K = .ntaxa(tree)
  phyloparent = numeric(tree$Nnode + K)
  phylochildren = matrix(0, tree$Nnode + K, 2)
  for (i in 1:nrow(tree$edge)) {
    i1 = tree$edge[i, 1]
    i2 = tree$edge[i, 2]
    phyloparent[i2] = i1
    if (phylochildren[i1, 1] == 0) 
      phylochildren[i1, 1] = i2
    else phylochildren[i1, 2] = i2
  }
  descendant = matrix(FALSE, tree$Nnode + K, K)
  for (i in 1:K) descendant[i, i] = TRUE
  processed = logical(tree$Nnode + K)
  processed[1:K] = TRUE
  while (!all(processed)) {
    for (i in (K + 1):(K + tree$Nnode)) {
      if (all(processed[phylochildren[i, ]])) {
        descendant[i, descendant[phylochildren[i, 1], ]] = TRUE
        descendant[i, descendant[phylochildren[i, 2], ]] = TRUE
        processed[i] = TRUE
      }
    }
  }
  list(phylotree = tree, phylochildren = phylochildren, 
       phyloparent = phyloparent, descendant = descendant)
}
# use Nexc most exterme test statistic to approximate GPD 
.gpd_approx <- function(yperm, nexc, obsts){
  yperm = na.omit(yperm)
  y = sort(yperm, decreasing = TRUE)
  tmp = .gpd_params_est(nexc, y)
  a_hat = tmp[1]
  k_hat = tmp[2]
  t = tmp[3]
  ### goodness-fit test
  nexc_re = .gpd_goft(nexc, y[1:nexc]-t)
  if(nexc_re[2] == nexc){
    z0 = obsts - t
    p = .gpd_pval(a_hat, k_hat, z0, length(yperm), nexc)
    return(p)
  }else{
    tmp = .gpd_params_est(nexc_re[2], y)
    a_hat = tmp[1]
    k_hat = tmp[2]
    t = tmp[3]
    z0 = obsts - t
    p = .gpd_pval(a_hat, k_hat, z0, length(yperm), nexc_re[2])
    return(p)
  }
}

# Parameter estimators for the generalized Pareto distribution
.gpd_params_est <- function(nexc, y){
  t = (y[nexc] + y[nexc+1])/2
  z = y[1:nexc] - t
  z2 = z^2
  m = mean(z)
  m2 = mean(z2)
  a_hat = m*m2/(m2-m^2)/2
  k_hat = (m^2/(m2-m^2)-1)/2
  return(c(a_hat, k_hat, t))
}

# p-value for the generalized Pareto distribution approximation
.gpd_pval <- function(a_hat, k_hat, z0, nperm, nexc){
  p = nexc/nperm*((1-k_hat*z0/a_hat)**(1/k_hat))
  return(p)
}

# iteratively reduce Nexc by 10 until a goodness-of-fit satisfy
.gpd_goft <- function(nexc, y){
  # y: the sorted test statistic - threshold t
  nexc = length(y)
  p = .gp_test(y)
  
  nexc.c = seq(0,nexc-10,10)
  z = y
  i = 0
  
  re = c()
  for(i in nexc.c) {
    z = y[1:(nexc-i)]
    p = .gp_test(z)
    re = rbind(re, c(i, p))
    if(!is.na(p) & p > 0.05) break
    i = i + 10
  }
  
  if(nrow(re) >=2) {
    nexc.c2 = seq(re[(nrow(re)-1),1]+1, re[nrow(re),1],1)
    re = c()
    for(i in nexc.c2) {
      z = y[1:(nexc-i)]
      p = .gp_test(z)
      re = rbind(re, c(i, p))
      if(!is.na(p) & p > 0.05) break
      i = i + 10
    }
  }
  
  p = re[nrow(re),2]
  len = nexc-re[nrow(re),1]
  
  return(c(p, len))
}

# Bootstrap goodness-of-fit test for the Generalized Pareto distribution
# test of fit for the GPD with unknown parameter
# refer to "goft"
.gp_test <- function(x, B = 2999){
  x = x[!is.na(x)]
  x = as.vector(x)
  n = length(x)   #  sample size without NA values
  samplerange = max(x) - min(x)
  gammap = .amle_method(x, k = ceiling(.2 * n))[1]    
  gamman = .combined_method(x)[1]
  r1 = .R1(x)     # observed value of R^-
  r2 = .R2(x)     # observed value of R^+
  # use "combined" for fitting GPD with negative shape parameter
  # use "amle" for fitting GPD with non-negative shape parameter
  p.value1 = sum(replicate(B, .R1(.rgp(n, shape = gamman))) < r1) / B  # bootstrap p-value for H_0^- 
  p.value2 = sum(replicate(B, .R2(.rgp(n, shape = gammap))) < r2) / B  # bootstrap p-value for H_0^+ 
  p.value = max(p.value1, p.value2)    # p-value of the intersection-union test
  return(p.value)
}

# Asymptotic maximum likelihood estimators 
.amle_method <- function(x, k){
  x = sort(x)
  n = length(x)
  nk = n - k
  x1 = x[(nk+1):n]
  w = log(x1)
  g = - (w[1] - sum(w) / k)    
  sigma = g * exp(w[1] + g * log(k / n))
  return(c(g, sigma))
}

# Combined estimators 
.combined_method <- function(x){
  m = mean(x)
  maxi = max(x)
  g = m / (m - maxi) 
  sigma = - g * maxi     
  return(c(g, sigma))
}

# Test statistic for H_0^-
.R1 <- function(x){
  gamma_neg = .combined_method(x)[1]
  Fn = ecdf(x)
  x1 = x[x != max(x)]
  z1 = (1 - Fn(x1))^( - gamma_neg) 
  return(abs(cor(x1, z1)))
}

# Test statistic for H_0^+
.R2  <- function(x){
  n = length(x)
  Fn = ecdf(x)
  gamma_positive = .amle_method(x, ceiling(.2 * n))[1]
  x1 = x[x != max(x)]
  y1 = (1 - Fn(x1))^( - gamma_positive) 
  x.star = log(x1)
  y.star = log( y1 -1 )
  if (gamma_positive <= 0.5)	return(cor(x1, y1))  
  if (gamma_positive > 0.5)  return((cor(x.star, y.star)))
}

# Simulation of random numbers from the gPd
.rgp  <- function (n, shape){
  if (shape != 0) 
    return((1 / shape) * (runif(n)^(-shape) - 1))
  else return(rexp(n, 1))
}
# Test alpha in mediator model
.test_alpha <- function(G, Trt, covariates, resid, disp){
  # Trt = tmp; covariates = conf; resid = mod.reduce.resid; disp = mod.reduce.disp
  Trt = matrix(Trt, nrow = length(Trt), ncol = 1)
  X1 = model.matrix(G~0+covariates)
  D.resid2 = diag(resid^2)/disp
  A11 = t(X1) %*% X1; B11 = t(X1) %*% D.resid2 %*% X1
  A12 = t(X1) %*% Trt; B12 = t(X1) %*% D.resid2 %*% Trt
  A21 = t(Trt) %*% X1; B21 = t(Trt) %*% D.resid2 %*% X1
  A22 = t(Trt) %*% Trt; B22 = t(Trt) %*% D.resid2 %*% Trt
  # continuous traits
  W = B22 - A21 %*% ginv(A11) %*% B12 - B21 %*% ginv(A11) %*% A12 + 
    A21 %*% ginv(A11) %*% B11 %*% ginv(A11) %*% A12
  U = as.vector(t(Trt) %*% resid)/disp
  V = W/disp
  stat = as.numeric(U %*% ginv(V) %*% U)
  pval = 1 - pchisq(stat, 1) 
  return(list(stat=stat, pval=pval))
}
# Test beta in outcome model
.test_beta <- function(outcome, G, Trt, conf, resid.obs=NULL, s2.obs=NULL, est.obs=NULL, obj=NULL, test.type="mv"){
  m = ncol(G)
  n = nrow(G)
  if(is.vector(G)){
    m = 1
    n = length(G)
  }
  ## get U and V
  if(test.type=="mv"){
    tmp = .cal_sumstats(outcome, G, cbind(conf[,-1], Trt), resid.obs, s2.obs, est.obs)
    U = tmp$U
    V = tmp$V
    
    stat = as.numeric(U %*% solve(V) %*% U)
    pval = 1 - pchisq(stat, m) 
  }
  
  if(test.type=="vc"){
    pval = SKAT(G, obj, is_check_genotype=FALSE, kernel="linear")$p.value
    stat = qchisq(1-pval, df=1)
  }
  ## estimate parameters
  if(length(unique(outcome)) > 2 ){
    # continuous trait
    est = summary(lm(outcome ~ G + cbind(conf[,-1], Trt)))$coefficients[1+1:m,1] 
  }else{
    # binary trait
    est = summary(glm(outcome ~ G + cbind(conf[,-1], Trt), family = "binomial"))$coefficients[1+1:m,1] 
  }
  
  return(list(stat=stat, est=est, pval=pval))
}

# modify Lan's code to Bowen: cal_score_sumstats.R 07/09/2019
# continuous: resid/s2 binary: est
.cal_sumstats <- function(Y, G, covariates, resid, s2, est){
  m = ncol(G)
  n = nrow(G)
  if(is.vector(G)){
    m = 1
    n = length(G)
  }
  if(!is.matrix(G)){
    G = matrix(G, nrow = n, ncol = m)
  }
  
  X1 <- model.matrix(Y~covariates)
  if(length(unique(Y)) > 2 ) {
    # continuous traits
    W = t(G) %*% G-
      (t(G) %*% X1) %*%
      solve(t(X1) %*% X1) %*% (t(X1) %*% G)
    
    obs.stat = list(U = as.vector(t(G) %*% resid) / s2, V = W/s2)
  } else {
    # binary traits
    mod = glm(Y~covariates, family = "binomial")
    sum.U = numeric(m)
    sum1 = matrix(0, m, m)
    sum2 = matrix(0, m, ncol(X1))
    sum3 = matrix(0, ncol(X1), ncol(X1))
    sum4 = matrix(0, ncol(X1), m)
    for(i in 1:n){
      lambdaX = as.numeric(est %*% X1[i, ])
      b1 = exp(lambdaX) / (1 + exp(lambdaX))
      b2 = exp(lambdaX) / ((1 + exp(lambdaX)))^2
      
      U.part1 = (Y[i] - b1) * G[i, ]
      U.part1[is.na(U.part1)] = 0
      V.part1 = b2 * G[i, ] %*% t(G[i, ])
      V.part1[is.na(V.part1)] = 0
      V.part2 = b2 * G[i, ] %*% t(X1[i, ])
      V.part2[is.na(V.part2)] = 0
      V.part3 = b2 * X1[i, ] %*% t(X1[i, ])
      V.part3[is.na(V.part3)] = 0
      V.part4 = b2 * X1[i, ] %*% t(G[i, ])
      V.part4[is.na(V.part4)] = 0
      
      sum.U = sum.U + U.part1
      sum1 = sum1 + V.part1
      sum2 = sum2 + V.part2
      sum3 = sum3 + V.part3
      sum4 = sum4 + V.part4
    }
    
    obs.stat = list(U = sum.U, V = sum1 - sum2 %*% solve(sum3) %*% sum4)
  }
  
  return(obs.stat)
}


.pi0_JC <- function(z){
  xi = c(0:100)/100
  tmax = sqrt(log(length(z)))
  tt = seq(0, tmax, 0.05)
  
  epsest = NULL
  
  for (j in 1:length(tt)) {
    t = tt[j]
    f = t*xi
    f = exp(f^2/2)
    w = (1 - abs(xi))
    co = 0*xi
    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z))
    }
    epshat = sum(w*f*co)/sum(w)
    epsest = c(epsest,epshat)
  }
  tmp = min(epsest)
  if(tmp > 1) tmp = 1
  return(tmp)
}

.nullEstimation_prod <- function (pval.alpha, pval.beta, alpha1, alpha2) {
  # alpha00: a=0 and b=0
  # alpha01: a=0 and b!=0
  # alpha10: a!=0 and b=0
  input.pvals = cbind(pval.alpha, pval.beta)
  idx.na = which(!complete.cases(input.pvals))
  input.pvals = input.pvals[complete.cases(input.pvals), ]
  
  alpha00 = alpha1 * alpha2
  alpha01 = alpha1 * (1-alpha2)
  alpha10 = (1-alpha1) * alpha2
  w = alpha00+alpha01+alpha10
  alpha00 = alpha00/w
  alpha01 = alpha01/w
  alpha10 = alpha10/w
  
  tmp = pmax(input.pvals[,1], input.pvals[,2])
  nmed = length(tmp)
  cdf12 = input.pvals
  if(length(which(input.pvals == 1)) > 0)
    input.pvals[input.pvals == 1] = input.pvals[input.pvals == 1] - runif(length(which(input.pvals == 1)), min = 0, max = 1e-7)
  input.pvals = input.pvals + runif(length(tmp)*2, min = 0, max = 1e-8)
  xx1 = c(0, input.pvals[order(input.pvals[, 1]), 1])
  yy1 = c(0, seq(1, nmed, by = 1)/nmed)
  gfit1 = gcmlcm(xx1, yy1, type = "lcm")
  xknots1 = gfit1$x.knots[-1]
  Fknots1 = cumsum(diff(gfit1$x.knots) * gfit1$slope.knots)
  xx2 = c(0, input.pvals[order(input.pvals[, 2]), 2])
  yy2 = c(0, seq(1, nmed, by = 1)/nmed)
  gfit2 = gcmlcm(xx2, yy2, type = "lcm")
  xknots2 = gfit2$x.knots[-1]
  Fknots2 = cumsum(diff(gfit2$x.knots) * gfit2$slope.knots)
  if (alpha1 != 1) 
    Fknots1 = (Fknots1 - alpha1 * xknots1)/(1 - alpha1)
  else Fknots1 = rep(0, length(xknots1))
  if (alpha2 != 1) 
    Fknots2 = (Fknots2 - alpha2 * xknots2)/(1 - alpha2)
  else Fknots2 = rep(0, length(xknots2))
  orderq1 = orderq2 = gcdf1 = gcdf2 = tmp
  for (i in 1:length(xknots1)) {
    if (i == 1) {
      gcdf1[orderq1 <= xknots1[i]] = (Fknots1[i]/xknots1[i]) * orderq1[orderq1 <= xknots1[i]]
    }
    else {
      if (sum(orderq1 > xknots1[i - 1] & orderq1 <= xknots1[i]) > 0) {
        temp = orderq1[orderq1 > xknots1[i - 1] & orderq1 <= xknots1[i]]
        gcdf1[orderq1 > xknots1[i - 1] & orderq1 <= 
                xknots1[i]] = Fknots1[i - 1] + (Fknots1[i] - 
                                                  Fknots1[i - 1])/(xknots1[i] - xknots1[i - 
                                                                                          1]) * (temp - xknots1[i - 1])
      }
    }
  }
  for (i in 1:length(xknots2)) {
    if (i == 1) {
      gcdf2[orderq2 <= xknots2[i]] = (Fknots2[i]/xknots2[i]) * orderq2[orderq2 <= xknots2[i]]
    }
    else {
      if (sum(orderq2 > xknots2[i - 1] & orderq2 <= xknots2[i]) > 0) {
        temp = orderq2[orderq2 > xknots2[i - 1] & orderq2 <= xknots2[i]]
        gcdf2[orderq2 > xknots2[i - 1] & orderq2 <= 
                xknots2[i]] = Fknots2[i - 1] + (Fknots2[i] - 
                                                  Fknots2[i - 1])/(xknots2[i] - xknots2[i - 
                                                                                          1]) * (temp - xknots2[i - 1])
      }
    }
  }
  gcdf1 = ifelse(gcdf1 > 1, 1, gcdf1)
  gcdf2 = ifelse(gcdf2 > 1, 1, gcdf2)
  cdf12[, 1] = gcdf1
  cdf12[, 2] = gcdf2
  rawp = (tmp * cdf12[, 2] * alpha01) + (tmp * cdf12[, 1] * alpha10) + (tmp^2 * alpha00)
  
  rawp.wNA = numeric(length(pval.alpha))
  if(length(idx.na) > 0){
    rawp.wNA[idx.na] = NA
    rawp.wNA[-idx.na] = rawp
  }else{
    rawp.wNA = rawp
  }
  
  rslt = list(alpha10 = alpha10, alpha01 = alpha01, alpha00 = alpha00, alpha1 = alpha1, alpha2 = alpha2, rawp = rawp.wNA)
  return(rslt)
}

.nullEstimation_minus <- function (pval.alpha, pval.beta, alpha1, alpha2, z1, z2) {
  # alpha00: a=0 and b=0
  # alpha01: a=0 and b!=0
  # alpha10: a!=0 and b=0
  
  input.pvals = cbind(pval.alpha, pval.beta)
  input.z = cbind(z1, z2)
  idx.na = which(!complete.cases(input.pvals))
  input.pvals = input.pvals[complete.cases(input.pvals), ]
  input.z = input.z[complete.cases(input.z), ]
  z.max = z.min = numeric(nrow(input.z))
  for (i in 1:nrow(input.z)) {
    z.max[i] = ifelse(abs(input.z[i,1]) >= abs(input.z[i,2]), input.z[i,1], input.z[i,2])
    z.min[i] = ifelse(abs(input.z[i,1]) <= abs(input.z[i,2]), input.z[i,1], input.z[i,2])
  }
  
  w = .pi0_JC(z.min)
  alphastar = min(w, 1)
  alpha00 = min(alpha1+alpha2-alphastar, 1)
  alpha10 = max(alphastar-alpha1, 0)
  alpha01 = max(alphastar-alpha2, 0)
  alpha11 = max(1-alphastar, 0)
  w = alpha00+alpha01+alpha10
  alpha00.r = alpha00/w
  alpha01.r = alpha01/w
  alpha10.r = alpha10/w
  alpha1.r = alpha00.r + alpha01.r
  alpha2.r = alpha00.r + alpha10.r
  
  # alpha1 = min(mean(input.pvals[,1]>=lambda)/(1-lambda), 1)
  # alpha2 = min(mean(input.pvals[,2]>=lambda)/(1-lambda), 1)
  # alpha00 = min(mean(input.pvals[,2]>=lambda & input.pvals[,1]>=lambda)/(1-lambda)^2, 1)
  # alpha10 = max(alpha2-alpha00, 0)
  # alpha01 = max(alpha1-alpha00, 0)
  
  tmp = pmax(input.pvals[,1], input.pvals[,2])
  nmed = length(tmp)
  cdf12 = input.pvals
  if(length(which(input.pvals == 1)) > 0)
    input.pvals[input.pvals == 1] = input.pvals[input.pvals == 1] - runif(length(input.pvals == 1), min = 0, max = 1e-7)
  input.pvals = input.pvals + runif(length(tmp)*2, min = 0, max = 1e-8)
  xx1 = c(0, input.pvals[order(input.pvals[, 1]), 1])
  yy1 = c(0, seq(1, nmed, by = 1)/nmed)
  gfit1 = gcmlcm(xx1, yy1, type = "lcm")
  xknots1 = gfit1$x.knots[-1]
  Fknots1 = cumsum(diff(gfit1$x.knots) * gfit1$slope.knots)
  xx2 = c(0, input.pvals[order(input.pvals[, 2]), 2])
  yy2 = c(0, seq(1, nmed, by = 1)/nmed)
  gfit2 = gcmlcm(xx2, yy2, type = "lcm")
  xknots2 = gfit2$x.knots[-1]
  Fknots2 = cumsum(diff(gfit2$x.knots) * gfit2$slope.knots)
  if (alpha1.r != 1) 
    Fknots1 = (Fknots1 - alpha1.r * xknots1)/(1 - alpha1.r)
  else Fknots1 = rep(0, length(xknots1))
  if (alpha2.r != 1) 
    Fknots2 = (Fknots2 - alpha2.r * xknots2)/(1 - alpha2.r)
  else Fknots2 = rep(0, length(xknots2))
  orderq1 = orderq2 = gcdf1 = gcdf2 = tmp
  for (i in 1:length(xknots1)) {
    if (i == 1) {
      gcdf1[orderq1 <= xknots1[i]] = (Fknots1[i]/xknots1[i]) * orderq1[orderq1 <= xknots1[i]]
    }
    else {
      if (sum(orderq1 > xknots1[i - 1] & orderq1 <= xknots1[i]) > 0) {
        temp = orderq1[orderq1 > xknots1[i - 1] & orderq1 <= xknots1[i]]
        gcdf1[orderq1 > xknots1[i - 1] & orderq1 <= 
                xknots1[i]] = Fknots1[i - 1] + (Fknots1[i] - 
                                                  Fknots1[i - 1])/(xknots1[i] - xknots1[i - 
                                                                                          1]) * (temp - xknots1[i - 1])
      }
    }
  }
  for (i in 1:length(xknots2)) {
    if (i == 1) {
      gcdf2[orderq2 <= xknots2[i]] = (Fknots2[i]/xknots2[i]) * orderq2[orderq2 <= xknots2[i]]
    }
    else {
      if (sum(orderq2 > xknots2[i - 1] & orderq2 <= xknots2[i]) > 0) {
        temp = orderq2[orderq2 > xknots2[i - 1] & orderq2 <= xknots2[i]]
        gcdf2[orderq2 > xknots2[i - 1] & orderq2 <= 
                xknots2[i]] = Fknots2[i - 1] + (Fknots2[i] - 
                                                  Fknots2[i - 1])/(xknots2[i] - xknots2[i - 
                                                                                          1]) * (temp - xknots2[i - 1])
      }
    }
  }
  gcdf1 = ifelse(gcdf1 > 1, 1, gcdf1)
  gcdf2 = ifelse(gcdf2 > 1, 1, gcdf2)
  cdf12[, 1] = gcdf1
  cdf12[, 2] = gcdf2
  rawp = (tmp * cdf12[, 2] * alpha01.r) + (tmp * cdf12[, 1] * alpha10.r) + (tmp^2 * alpha00.r)

  rawp.wNA = numeric(length(pval.alpha))
  if(length(idx.na) > 0){
    rawp.wNA[idx.na] = NA
    rawp.wNA[-idx.na] = rawp
  }else{
    rawp.wNA = rawp
  }
  rslt = list(alpha10 = alpha10.r, alpha01 = alpha01.r, alpha00 = alpha00.r, alpha1 = alpha1, alpha2 = alpha2, rawp = rawp.wNA)
  return(rslt)
}

# # determines the maximal number of permutations 
# .choose_b <- function(alpha, c) {
#   error <- alpha * c
#   B <- alpha*(1 - alpha) / (error^2)
#   return(B)
# }

# determines the number of test statistics that should 
# be sampled in adaptive permutation 
# 68% CI, 1 SE
.choose_r <- function(alpha, c) {
  error <- alpha * c
  R <- 0
  foundR <- FALSE
  while(!foundR) {
    R <- R + 1
    brange <- qnbinom(c(0.1586553, 0.8413447), R, alpha)
    pvalRange <- R / (R + brange)
    diff <- max(abs(pvalRange - alpha))
    if(diff < error) {
      foundR <- TRUE
    }
  }
  return(R)
}

# variance component test for alpha
# Ref: YT Huang and X Lin, BMC Bioinformatics 2013, 14:210
# Ref: YT Huang, Biometrics 2019, https://doi.org/10.1111/biom.13073
.test_alpha_vc <- function(G, Trt, covariates){
  G.res = resid(lm(G~0+covariates))
  Trt.res = resid(lm(Trt~0+covariates))
  Y1 = t(G.res)
  X = Trt.res
  n = dim(Y1)[2]
  p = dim(Y1)[1]
  Y = as.numeric(as.matrix(Y1))
  
  ## Calculate working covariance
  covY.ind = matrix(0, p, p)
  diag(covY.ind) = 1
  v.work = covY.ind
  
  ## Calculate observed Q
  Y1.bar = apply(Y1, 1, mean)
  v.work.inv = chol2inv(chol(v.work))	
  Q.half = 0
  for (i in 1:n) Q.half = Q.half+(X[i]*t(Y1[,i]-Y1.bar)%*%v.work.inv)
  Q = sum(Q.half^2)	
  
  ## Calculate Q under the null by permutation
  Q.perm = NULL
  Qj.perm = NULL
  set.seed(126)
  for (pp in 1:1000){
    x.perm = sample(X)
    Q.temp.perm = t(as.numeric(as.matrix(Y1-Y1.bar))/rep(1, n))*rep(x.perm, each=p)
    Q.half.perm = apply(matrix(Q.temp.perm, nrow=p), 1, sum)
    Q.perm = c(Q.perm, sum(Q.half.perm^2))
  }
  Q.perm.sub = Q.perm
  
  ## Calculate p-value scaled chi-square p-value
  stat = Q*2*mean(Q.perm.sub)/var(Q.perm.sub)
  pval = pchisq(stat, df = 2*mean(Q.perm.sub)^2/var(Q.perm.sub), lower.tail=F)
  return(list(stat = stat, pval = pval))
}

.test_beta_vc <- function(outcome, G, Trt, conf, obj){
  m = ncol(G)
  pval = SKAT(G, obj, is_check_genotype=FALSE, kernel="linear")$p.value
  stat = qchisq(1-pval, df=1)
  if(length(unique(outcome)) > 2 ){
    # continuous trait
    est = summary(lm(outcome ~ G + cbind(conf[,-1], Trt)))$coefficients[1+1:m,1] 
  }else{
    # binary trait
    est = summary(glm(outcome ~ G + cbind(conf[,-1], Trt), family = "binomial"))$coefficients[1+1:m,1] 
  }
  return(list(stat=stat, est=est, pval=pval))
}



.prepareTree <- function(tree, verbose = FALSE){
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
