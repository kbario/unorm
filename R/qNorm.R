#' @title Quantile Normalisation
#' @description A method of NMR spectral normalisation
#' @details ### How It Works:
#' ### Advantages:
#' ### Limitations:
#' @family {}
#' @param
#' @param
#' @return
#' @seealso
#' @author
#' @examples
#' @importFrom
#' @export
#'

qNorm <- function(X){
  Xs <- t(apply(X, 1, sort))
  Xo <- t(apply(X, 1, order))

}

xs <- Xs[1,]
bins <- seq(min(xs), max(xs), by = 500)

xsbin <- cut(xs,breaks = bins)
cs <- table(xsbin)

cum <- sapply(1:length(cs), function(x){
  pr <- (sum(cs[1:x]))/(sum(cs))
  return(pr)
})

# smaller
xs <- Xs[1,]
cs <- table(cut(xs,breaks = seq(min(xs), max(xs), by = 500)))
cdef <- sapply(1:length(cs), function(x){
  pr <- (sum(cs[1:x]))/(sum(cs))
  return(pr)
})


cd <- function(X, num = 2){sapply(1:num, function(x){
  maxx <- max(Xs)
  minn <- min(Xs)
  xs <- Xs[x,]
  cs <- table(cut(xs,breaks = seq(minn, maxx, by = 500)))
  cdef <- sapply(1:length(cs), function(x){
    pr <- (sum(cs[1:x]))/(sum(cs))
    return(pr)
    })
  return(cdef)
})}

cdef <- cd(Xs, 2)

plot(cdef[,1], cdef[,2])


















