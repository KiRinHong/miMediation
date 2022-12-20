#' Real data
#' 
#' @description This data is derived from a real microbiome dataset (Zeevi, D. et al., 2015), which includes microbiome samples from 200 healthy subjects.
#' @docType data
#' @usage data(data.zeeviD)
#' 
#' @format \code{data.zeeviD} contains the following components:
#' \describe{
#'   \item{treatment}{treatment indicator.}
#'   \item{mediators}{an abundance matrix with the top 100 most abundant taxa.}
#'   \item{outcome}{continuous outcome.}
#'   \item{tree}{a taxonomy table.}
#' }
#' @references 
#' Zeevi, D. et al. (2015). 
#' Personalized nutrition by prediction of glycemic responses. 
#' \emph{Cell} 163:1079-1094.
#' @source \url{https://doi.org/10.1016/j.cell.2015.11.001}
#' @keywords dataset
#' 
#' @examples 
#' data(data.zeeviD)
"data.zeeviD"
