#' @title Creatinine Normalisation
#' @description Creatinine Normalisation (CN) is a useful method much like region of interest normalisation that can normalise spectra based on the total area of the creatinine signal at the chemical shift 3.05ppm.
#' @details
#' ### How It Works:
#' 1. For each experimental spectra, the integral (total area) of the region defined by the bounds in the argument `cre3` is calculated.
#' 2. Each intensity in that spectra is then divided by the integral. This is iterated over for all spectra.
#' ### Advantages:
#' * CN is computationally light and is not as difficult as other normalisations like [pqNorm()]
#' * Can be useful for people who have not done any physical exercise recently or of similar age, gender and muscle mass.
#' ### Limitations:
#' * Complications arise from inter- and intra-individual creatinine signal fluctuation caused by exercise (which causes muscle breakdown and creatinine production) when using CN.
#' * CN relies on the idea that creatinine concentrations do not change in urine samples, the only changes to creatinine signal are created by effects of dilution but this is incorrect.
#' * Creatinine concentrations fluctuate depending on a person's physiology and behaviour so use CN consciously. Consider using in conjunction with a reference-based normalisation like [pqNorm()].
#' @family {Attribute-Based}
#' @param X A numerical matrix with rows being the spectra and columns being the chemical shift variables.
#' @param ppm A numerical array holding the chemical shift values of the X matrix.
#' @param cre3 A concatenated numerical value of the lower and upper ppm values where the creatinine peak at 3.05 starts and ends.
#' @param cre4 A concatenated numerical value of the lower and upper ppm values where the creatinine peak at 4.05 starts and ends.
#' @return A list of:
#' 1. The normalised X matrix,
#' 2. A numerical array of the corresponding dilution factors, and
#' 3. The ratio of the creatinine signals (ideally should equal 0.667).
#' * Following the example below will extract the results quickly and easily.
#' @author \email{kylebario1@@gmail.com}
#' @seealso More on the methodology of CN and issue with using it are outlined here: \url{https://doi.org/10.1021/ac051632c}
#' @examples
#' Xcn <- creNorm(X)
#' Xcre <- Xcn[[1]]
#' creDilf <- Xcn[[2]]
#' @export
#' @importFrom metabom8 get_idx

creNorm <- function(X, ppm, cre3 = c(3.043, 3.055), cre4 = c(4.05,4.07)){
  i3 <- get_idx(cre3, ppm)
  i4 <- get_idx(cre4, ppm)
  cre_a3 <- apply(X[,i3], 1, sum)
  cre_a4 <- apply(X[,i4], 1, sum)
  ratio <- cre_a4/cre_a3
  Xn <- sapply(1:nrow(X), function(i){
    X[i,]/cre_a3[i]
  })
  return(list(Xn, cre_a3, ratio))
}
