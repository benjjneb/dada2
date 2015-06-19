#' Plot Substitution Pairs from DADA Result
#' 
#' This is similar to original DADA article, Figure 6.
#' 
#' @param dadaOut (Required). The object returned by \code{\link{dada}}.
#' 
#' @param facetByGrp (Optional). Logical(1).
#'  Whether to plot all substitution groups together in one panel
#'  or separately on a grid of panels with a linear model fit.
#' 
#' @return A \code{\link{ggplot}2} object that will be rendered
#'  to default device if \code{\link{print}ed},
#'  or can be stored and further modified.
#'  See \code{\link{ggsave}} for additional options.
#' 
#' @importFrom reshape2 melt
#' @importFrom reshape2 dcast
#' @import data.table
#' @import ggplot2
#' 
#' @export
#' 
#' @examples 
#' # Examples here.
#' testFile = system.file("extdata", "test-nonunique.fastq.gz", package="dadac")
#' test1 = derepFastq(testFile, verbose = TRUE)
#' test1$quals[test1$quals > 40] <- 40
#' res1 <- dada(uniques = test1$uniques, quals = test1$quals,
#'              err = dadac:::inflateErr(tperr1, 2), 
#'              OMEGA_A = 1e-40, 
#'              USE_QUALS = TRUE, 
#'              err_function = makeErrfun(threepiece, dadac:::tp_init_parms), 
#'              self_consist = TRUE) 
#' plot_substitutions(res1)
plot_substitutions = function(dadaOut, facetByGrp = TRUE){
  transdt = data.table(melt(dadaOut$trans))
  setnames(transdt, c("Substitution", "Quality", "Count"))
  transdt[, Sub1 := substr(Substitution, 1, 1)]
  transdt[, Sub2 := substr(Substitution, 3, 3)]
  # For a plot like from the DADA-1 paper, Figure 6,
  # Map color to complementary pairs of errors
  # red = (A→G,T→C) cyan=(C→T,G→A) green=(A→T,T→A) black=(C→A,G→T) blue=(A→C,T→G) purple=(C→G,G→C).
  CompSubGroups = c(A2G = "A2GT2C", T2C = "A2GT2C",
                    C2T = "C2TG2A", G2A = "C2TG2A",
                    A2T = "A2TT2A", T2A = "A2TT2A",
                    C2A = "C2AG2T", G2T = "C2AG2T",
                    A2C = "A2CT2G", T2G = "A2CT2G",
                    C2G = "C2GG2C", G2C = "C2GG2C") 
  transdt[, SubGrp := CompSubGroups[as.character(Substitution)]]
  # Redefine the Forward/Reverse pairings
  transdt[(substr(SubGrp, 1, 3) == Substitution), Direction := "Forward"]
  transdt[(substr(SubGrp, 4, 6) == Substitution), Direction := "Reverse"]
  # Define plot as in DADA1 article Figure 6
  # (but better because it is not matlab)
  tCast = dcast.data.table(data = transdt[(Direction != "NoChange")][Count > 0],
                           formula = SubGrp + Quality ~ Direction,
                           value.var = "Count",
                           drop = TRUE)
  # Remove missing
  tCast <- tCast[!is.na(Forward)][!is.na(Reverse)]
  # Define ggplot2 plot
  p1 = ggplot(tCast, aes(Forward, Reverse, color = SubGrp, size = Quality)) + 
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(intercept=0, slope=1, alpha = 0.5, linetype = 2) +
    geom_point(alpha = 0.7) + 
    ggtitle("Substitution Pairs Count Comparison")
  if(facetByGrp){
    # If TRUE, facet and add simple linear model fit
    p1 <- p1 + 
      stat_smooth(method = "lm") +
      facet_wrap(~SubGrp, drop = TRUE)
  }
  return(p1)
}
