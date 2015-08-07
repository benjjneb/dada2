############################################################################
#' derep - A class representing dereplicated sequences
#' 
#' A \code{\link{list}} with the following three members.
#' \itemize{
#'  \item{$uniques: Named integer vector. Named by the unique sequence, valued by abundance.}
#'  \item{$quals: Numeric matrix of average quality scores by position for each unique. Uniques are rows, positions are cols.}
#'  \item{$map: Integer vector of length the number of reads, and value the index (in $uniques) of the unique to which that read was assigned.}
#' }
#' This can be created from a FastQ sequence file using
#' \code{\link{derepFastq}}
#' 
#' @seealso \code{\link{derepFastq}}
#' 
#' @name derep-class
#' @rdname derep-class
setClass("derep", contains = "list")
############################################################################
#' The object class returned by \code{\link{dada}}
#'
#' A multi-item List with the following named values...
#' \itemize{
#'  \item{$genotypes: }{Integer vector, named by sequence valued by abundance, of the denoised genotypes.}
#'  \item{$clustering: }{An informative data.frame containing information on each cluster.}
#'  \item{$quality: }{The average quality scores for each cluster (row) by position (col).}
#'  \item{$map: }{Integer vector that maps the unique (index) to the cluster/genotype (value).}
#'  \item{$birth_subs: }{A data.frame containing the substitutions at the birth of each new cluster.}
#'  \item{$trans: }{The matrix of transitions by type (row), eg. A2A, A2C..., and quality score (col)
#'          observed in the final output of the dada algorithm.}
#'  \item{$err_in: }{The err matrix used for this invocation of dada.}
#'  \item{$err_out: }{The err matrix estimated from the output of dada. NULL if err_function not provided.}
#'  \item{$opts: }{A list of the dada_opts used for this invocation of dada.}
#'  \item{$uniques: }{The uniques vector(s) used for this invocation of dada.}
#'  \item{$call: }{The function call used for this invocation of dada.}
#' }
#' 
#' @seealso \code{\link{dada}}
#' 
#' @name dada-class
#' @rdname dada-class
setClass("dada", contains = "list")