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
#' cre <- creNorm(X, ppm)
#' Xn <- cre$Xn
#' dilf <- cre$dilf
#' creRatio <- cre$ratio
#' @export
#' @importFrom metabom8 get_idx

creNorm <- function(X, ppm = NULL, cre3 = c(3, 3.1), cre4 = c(4, 4.1), err = 5){
  if (is.null(dim(X))){
    if (is.null(length(X))){
      stop("Please provide a valid X variable. X is neither a matrix or an array")
    }
    if (is.null(ppm)){
      stop("Please provide a X-matched ppm. None was provided and ppm cannot be determined from a single spectrum")
    }
    cat('\033[0;34mCalculating Dilfs... \033[0m')
    i3 <- shift_pickr(X, ppm, cre3, 0.005)
    i4 <- shift_pickr(X, ppm, cre4, 0.005)
    a3 <- sum(X[i3])
    a4 <- sum(X[i4])
    r <- a4/a3
    cat('\033[1;32mDone.\n\033[0m')
    cat('\033[0;34mChecking creatinine peak ratio... \033[0m')
    er <- ((2/3)/100)*err
    lo <- (2/3)-er
    up <- (2/3)+er
    if(r<=up & r>=lo){
      cat('\033[1;32mRatio is within limit.\n\033[0m')
    } else{
      cat('\033[1;31mThe provided spectra is outside of the error limit.\n\033[1;33mRefer to Df_cre for more information.\n\033[0m')
    }
    e <- as.array(r<=up & r>=lo)
    df <- data.frame(a3, r, e)
    colnames(df) <- c('dilf', 'ratio', paste('ratio within a', err, '% error margin'))
    cat('\033[0;34mNormalising X... \033[0m')
    Xn <- X/a3
    rownames(Xn) <- rownames(X)
    cat('\033[1;32mDone.\n\033[0m')
    assign("X_cre", Xn, envir = .GlobalEnv)
    assign("dilf_cre", a3, envir = .GlobalEnv)
    assign("Df_cre", df, envir = .GlobalEnv)
  }
  if (!is.null(dim(X))){
    cat('\033[0;34mCalculating Dilfs... \033[0m')
    if (is.null(ppm)){
      p <- as.numeric(colnames(X))
    } else {
      p <- ppm
    }
    i3 <- shift_pickr(X, p, cre3, 0.005)
    i4 <- shift_pickr(X, p, cre4, 0.005)
    a3 <- sapply(1:nrow(X),function(i){
      j <- i3[i,]
      a <- sum(X[i,j])
      return(a)
      })
    a4 <- sapply(1:nrow(X),function(i){
      j <- i4[i,]
      a <- sum(X[i,j])
      return(a)
    })
    r <- a3/a4
    cat('\033[1;32mDone.\n\033[0m')
    cat('\033[0;34mChecking creatinine peak ratios... \033[0m')
    er <- ((2/3)/100)*err
    lo <- (2/3)-er
    up <- (2/3)+er
    if(all(r<=up & r>=lo)){
      cat('\033[1;32mAll within limits.\n\033[0m')
    } else{
      cat('\033[1;31mspec', which(r>=up | r<=lo), 'are outside error limits.\n\033[1;33mRefer to Df_cre for more information.\n\033[0m')
    }
    e <- as.array(r<=up & r>=lo)
    df <- data.frame(a3, r, e)
    colnames(df) <- c('dilf', 'ratio', paste('ratio within a', err, '% error margin'))
    cat('\033[0;34mNormalising X... \033[0m')
    Xn <- t(sapply(1:nrow(X), function(i){
      X[i,]/a3[i]
    }))
    rownames(Xn) <- rownames(X)
    cat('\033[1;32mDone.\n\033[0m')
    assign("X_cre", Xn, envir = .GlobalEnv)
    assign("dilf_cre", a3, envir = .GlobalEnv)
    assign("Df_cre", df, envir = .GlobalEnv)
  }
}
