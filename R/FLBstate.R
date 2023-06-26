#' FLBstate
#'
#' Data for the functioning of the B cell state prediction.
#'
#' @name FLBstate
#' @docType package
NULL

#' FL B cell Signature Matrix
#'
#' A data frame representing linear defining gene signatures for each of the three B cell states
#'
#' @docType data
#' @name signature_mat
#' @usage data(signature_mat)
#' @format A data frame with 5 columns and 150 rows
NULL

#' FL Bulk Expression Example Matrix
#'
#' A data frame representing gene expression of the 150 signature genes across 91 samples. The values are log2TPM+1 values.
#'
#' @docType data
#' @name test_exp
#' @usage data(test_exp)
#' @format A data frame with 91 columns and 150 rows
NULL
