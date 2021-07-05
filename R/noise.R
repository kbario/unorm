#' @title Maximum Noise Estimation
#' @description This function estimates the maximum level of noise for each spectra in an X matrix. The output of this function can be used as the noi variable in both [pqNorm()] and [hmNorm()]
#' @details This function uses a combination of the functions mean and sd to find the standard deviation of the noise region as well as it's mean. The standard deviation is multiplied by five to capture any extreme outliers and added to the mean noise level to produce an estimation of maximum noise. The output of this function is intended for use in both the [pqNorm()] and [hmNorm()] functions so they can remove the noise from the spectra and produce clear results.
#' @family {preproc}
#' @param Xo The numerical matrix containing the NMR data you wish to estimate the noise of. The rows must contain information of one whole spectrum and the columns contain the specific chemical shift variables.
#' @param ppm An array of chemical shift variables. ppm should be column matched to the X matrix you using.
#' @param shift The concatenated ppm values that define the lower and upper bounds of the noise region. Default is c(9.5,11)
#' @param sd_mult The value that the standard deviation will be multiplied by. Default = 5. This does not need to be changed.
#' @return This function returns an estimation of the maximum noise level for each spectra.
#' @seealso Noise estimation methodology was adapted from histogram matching methods paper which can be found here: \url{http://dx.doi.org/10.1007/s11306-018-1400-6}
#' @author \email{kylebario1@@gmail.com}
#' @examples
#' noi <- noise(X, c(9.5,11))
#' @importFrom metabom8 get_idx
#' @importFrom stats sd
#' @export

noise <- function(Xo, ppm, shift = c(9.5,11), sd_mult = 5){
  rm <- apply(X, 1, function(i){
    rm <- sd_mult*sd(i[get_idx(shift, ppm)])+(mean((i[get_idx(shift, ppm)]), trim = 0.05))
    return(rm)
  })
  return(rm)
}
