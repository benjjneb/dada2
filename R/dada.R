
#` The main interface to DADA from R.
#` Type and sanity checks arguments (TODO).
#` 
#' @export
#'
dada <- function(uniques,
                 err = matrix(c(0.991, 0.003, 0.003, 0.003, 0.003, 0.991, 0.003, 0.003, 0.003, 0.003, 0.991, 0.003, 0.003, 0.003, 0.003, 0.991), nrow=4, byrow=T),
                 score = matrix(c(5, -4, -4, -4, -4, 5, -4, -4, -4, -4, 5, -4, -4, -4, -4, 5), nrow=4, byrow=T),
                 gap_penalty = -8,
                 self_consist = FALSE) {
  # ADD TYPE VALIDATION HERE
  prev <- NULL
  cur <- dada_uniques(names(uniques), unname(uniques), err, score, gap_penalty)
  while(self_consist && !identical(cur, prev)) {
    prev <- cur
    err <- cur$trans + 1   # ADD ONE PSEUDOCOUNT TO EACH TRANSITION
    err <- t(apply(err, 1, function(x) x/sum(x)))  # apply returns a transposed result
    cur <- dada_uniques(names(uniques), unname(uniques), err, score, gap_penalty)
  }
  cur
}
