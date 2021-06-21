#' @title Creatinine Normalisation
#' @description Creatinine Normalisation (CN) is a useful method much like region of interest normalisation that can normalise spectra based on the total area of the creatinine signal at the chemical shift 3.05ppm.
#' @param X A numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @param ppm A numerical array holding the chemical shift values of the X matrix
#' @param cre3 A concatenated numerical value of the lower and upper ppm values where the creatinine peak at 3.05 starts and ends
#' @param cre4 A concatenated numerical value of the lower and upper ppm values where the creatinine peak at 4.05 starts and ends
#' @return A list of the normalised X matrix, a numerical array of the corresponding dilution factors, and ratio of the creatinine signals (ideally should equal 0.667)
#' @details #Limitations:
#' * While useful, there are complicating factors associated with CN's use such as inter- and intra-individual creatinine signal fluctuation which hinders accurate normalisation.
#' * CN relies on the idea that creatinine concentrations do not change in urine samples, the only changes to creatinine signal are created by effects of dilution but this is incorrect.
#' * Creatinine concentrations fluctuate depending on a person's physiology and behaviour so use CN consciously.
#' @author \email{kylebario1@@gmail.com}
#' @family {Attribute-Based}
#' @examples
#' \dontrun{
#'   Xn <- creNorm(X)
#'   creDilf <- Xn[[2]]
#'   Xn <- Xn[[1]]}
#' @export
#' @importFrom metabom8 get_idx

creNorm <- function(X, ppm, cre3 = c(3.043, 3.055), cre4 = c(4.05,4.07)){
  i3 <- get_idx(cre3, ppm)
  i4 <- get_idx(cre4, ppm)
  cre_a3 <- apply(X[,i3], 1, sum)
  cre_a4 <- apply(X[,i4], 1, sum)
  cre_r <- cre_a4/cre_a3
  Xn <- sapply(1:nrow(X), function(i){
    X[i,]/cre_a3[i]
  })
  return(list(Xn, cre_a3, cre_r))
}
