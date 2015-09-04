#' seqMetaPipeline: Analysis Pipeline for SeqMeta
#'
#' Functions to make seqMeta analyses and reports consistent.
#'
#'
#' To learn more about seqMetaPipeline, start with the vignettes:
#'
#' @name seqMetaPipeline
#' @docType package
#' @import seqMeta plyr gap
#' @importFrom survival coxph
NULL



#' Phenotype data for 5700 subjects
#'
#' A dataset containing the phenotype data for 5700 subjects
#'
#' @format A data frame with 5718 rows and 5 variables:
#' \describe{
#'   \item{id}{subject identifier}
#'   \item{y}{phenotype measure}
#'   \item{age}{age of the subject}
#'   \item{pc1}{principle component 1}
#'   \item{pc2}{principle component 2}
#' }
#' @source Internal randomized dataset
"pheno1"