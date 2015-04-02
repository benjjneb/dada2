#' @export
threepiece <- function(parms, qave) {
  if(!(is.numeric(parms) && length(parms) == 4)) {
    stop("Invalid parms argument to threepiece.")
  }
  e0 <- parms[1]
  knot1 <- parms[2]
  knot2 <- parms[3]
  e40 <- parms[4]
  pred <- rep(e0, length(qave))
  pred[qave>knot1] <- e0 + (qave[qave>knot1]-knot1)*(e40-e0)/(knot2-knot1)
  pred[qave>knot2] <- e40
  pred <- 10^pred
  return(pred)
}

#' @export
make_log10_lsq_obj <- function(err_model, qave, obs) {
  objfun <- function(parms) {
    pred <- err_model(parms, qave)
    return(sum((log10(pred)-log10(obs))^2))
  }
  return(objfun)
}
# might need to step the log

tp_init_parms <- c(-1, 20, 35, -4)

#' @export
showErrors <- function(dq, nti, ntj, erri=TRUE, erro=TRUE) {
  ACGT <- c("A", "C", "G", "T")
  tij <- 4*(which(ACGT==nti)-1) + which(ACGT==ntj)
  nij <- paste0(nti,"2",ntj)
  subdf <- as.data.frame(t(dq$subqual))
  subdf$qave <- as.numeric(rownames(subdf))
  subdf[,nij] <- subdf[,nij]/(subdf[,paste0(nti,"2A")]+subdf[,paste0(nti,"2C")]+subdf[,paste0(nti,"2G")]+subdf[,paste0(nti,"2T")])

  erridf <- dq$err_in
  if(is.list(erridf)) erridf <- erridf[[1]]
  erridf <- as.data.frame(t(erridf))
  colnames(erridf) <- paste0(rep(ACGT, each=4), "2", ACGT)
  rownames(erridf) <- seq(dq$opts$QMIN, dq$opts$QMAX, length.out=nrow(erridf))
  erridf$qave <- as.numeric(rownames(erridf))
  
  errodf <- as.data.frame(t(dq$err_out))
  errodf$qave <- as.numeric(rownames(errodf))
  
  p <- ggplot(data=subdf, aes_string(x="qave", y=paste0("log10(",nij,")")))
  p <- p + geom_point()
  if(erri) p <- p + geom_line(data=erridf)
  if(erro) p <- p + geom_line(data=errodf, linetype="dashed")
  return(p)
}
