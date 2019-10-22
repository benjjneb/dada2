############################################################################
#' A class representing dereplicated sequences
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
#'  \item{$denoised: }{Integer vector, named by sequence valued by abundance, of the denoised sequences.}
#'  \item{$clustering: }{An informative data.frame containing information on each cluster.}
#'  \item{$sequence: }{A character vector of each denoised sequence. Identical to names($denoised).}
#'  \item{$quality: }{The average quality scores for each cluster (row) by position (col).}
#'  \item{$map: }{Integer vector that maps the unique (index of derep$unique) to the denoised sequence (index of dada$denoised).}
#'  \item{$birth_subs: }{A data.frame containing the substitutions at the birth of each new cluster.}
#'  \item{$trans: }{The matrix of transitions by type (row), eg. A2A, A2C..., and quality score (col)
#'          observed in the final output of the dada algorithm.}
#'  \item{$err_in: }{The err matrix used for this invocation of dada.}
#'  \item{$err_out: }{The err matrix estimated from the output of dada. NULL if err_function not provided.}
#'  \item{$opts: }{A list of the dada_opts used for this invocation of dada.}
#' }
#' 
#' @seealso \code{\link{dada}}
#' 
#' @name dada-class
#' @rdname dada-class
setClass("dada", contains = "list")
############################################################################
#' The named integer vector format used to represent collections of unique DNA sequences.
#'
#' The uniques vector is an \code{integer} vector that is named by the unique sequence, and 
#' valued by the abundance of that sequence. This format is commonly used within the 
#' \code{\link{dada2-package}}, for function inputs and outputs. The \code{\link{getUniques}}
#' function coerces a variety of input objects into the uniques-vector format, including
#' \code{\link{dada-class}} and \code{\link{derep-class}} objects.
#' 
#' @seealso \code{\link{getUniques}}
#' 
#' @name uniques-vector
#' @rdname uniques-vector
NULL