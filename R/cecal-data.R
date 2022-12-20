#' Real data
#' 
#' @description This data is derived from a real microbiome dataset (Cho et al., 2012), which includes 48 samples (38 in antibiotics vs. 10 in controls).
#' @docType data
#' @usage data(data.cecal)
#' 
#' @format \code{data.cecal} contains the following components:
#' \describe{
#'   \item{treatment}{treatment indicator: Antibiotics group and control group are coded as 1 and 0, respectively.}
#'   \item{mediators}{an abundance matrix with the top 100 most abundant taxa that have at least 20\% non-zero observations.}
#'   \item{outcome}{body fat percentage.}
#'   \item{tree}{a phylogenetic tree.}
#' }
#' @references 
#' Cho, I. et al. (2012). 
#' Antibiotics in early life alter the murine colonic microbiome and adiposity. 
#' \emph{Nature} 488:621-626.
#' @source \url{https://doi.org/10.1038/nature11400}
#' @keywords dataset
#' 
#' @examples 
#' data(data.cecal)
"data.cecal"
