rm(list = ls())
load("/Users/kylebario/OneDrive - Murdoch University/R/rproj/smol_bins/osmoPipelineSmolscans.rda")


# xfactor norm
Xxf <- xfNorm(X1, osmoNew$Result)
xfDilf = Xxf[[2]]
Xxf = t(Xxf[[1]])

Xxf_area <- apply(Xxf, 1, sum)
plot(Xxf_area)

matspec(Xxf, ppm1)






# pqNorm
Xpqn <- pqNorm(X1)
pqnDilf <- Xpqn[[2]]

### things don't seem to be working
### lets try something else

pqNorm <- function(X){
  Xm <- apply(X, 2, median)
  dilf <- sapply(1:nrow(X), function(x){
    d <- density((X[x,]/Xm), from = -10, to = 10)
    dilf <- d$x[which.max(d$y)]
    return(dilf)
  })
  Xn <- sapply(1:length(dilf), function(x){
    X[x,]/dilf[x]
  })
  return(list(Xn, dilf))
}


Xm <- apply(X1, 2, median)

dilf <- sapply(1:nrow(X1), function(x){
  d <- density((X1[x,]/Xm), from = -10, to = 10)})
  dilf <- d$x[which.max(d$y)]
  return(dilf)
})

density(X1[2,]/Xm, from = -10, to = 10)
plot(density(X1[2,]/Xm, from = -10, to = 10))
d <- density(X1[2,]/Xm, from = -10, to = 10)
which_d <- d$x[which.max(d$y)]

plot(density(X1[2,]/Xm, from = -10, to = 10))
abline(v = which_d)


plot(X1[2,]/Xm)
X1_2_hist <- hist(X1[2,]/Xm)

max(X1_2_hist$counts)
X1_2_hist$breaks
max(X1_2_hist$density)


X1_2_hist$breaks[(which.max(X1_2_hist$counts)-1)]
X1_2_hist$breaks[(which.max(X1_2_hist$counts))]
X1_2_hist$breaks[(which.max(X1_2_hist$counts)+1)]

length(X1_2_hist$breaks[(which.max(X1_2_hist$counts)-1)]:X1_2_hist$breaks[(which.max(X1_2_hist$counts)-1)])

vector <- c(min(X1_2_hist$breaks), seq(from = X1_2_hist$breaks[(which.max(X1_2_hist$counts)-1)], to = X1_2_hist$breaks[(which.max(X1_2_hist$counts)+1)], length.out = 100), max(X1_2_hist$breaks))


X1_2_hist2 = hist(X1[2,]/Xm, breaks = vector)

X1_2_hist2$breaks[which.max(X1_2_hist2$counts)]
X1_2_hist2$breaks[which.max(X1_2_hist2$counts)-1]
X1_2_hist2$breaks[which.max(X1_2_hist2$counts)+1]

vector <- c(min(X1_2_hist2$breaks), seq(from = X1_2_hist2$breaks[which.max(X1_2_hist2$counts)-1], to = X1_2_hist2$breaks[which.max(X1_2_hist2$counts)+1], length.out = 100), max(X1_2_hist2$breaks))

X1_2_hist3 = hist(X1[2,]/Xm, breaks = vector)
X1_2_hist3$mids[which.max(X1_2_hist3$counts)]

X1_2_hist3$breaks[which.max(X1_2_hist3$counts)-1]
X1_2_hist3$breaks[which.max(X1_2_hist3$counts)+1]


# trying to find the area of the histogram that holds over 95% of the counts
sort(X1_2_hist3$density, decreasing = T)


# instead of finding the dist that has at least 95% of the counts,
# define points between 0 and 5 and see if it contains over 95% of the stuff
# creat a histogram
hgram <- hist((X1[6,]/Xm))
#create a vector with breaks 501 breaks (100 breaks between for each 1) between 0 and 5
victor <- c(min(hgram$breaks), seq(from = -0.005, to = 5.005, length.out = 502), max(hgram$breaks))
hgram <- hist((X1[6,]/Xm), breaks = victor)

plot(hgram, xlim = c(0,5))
hgram$mids

sum(h$counts)
sum(h$counts[c(2:501)])


## SHOWS % OF COUNTS BETWEEN 0 AND 5
((sum(h$counts[c(2:501)]))/(sum(h$counts)))*100



hgram$mids[(which.max(hgram$counts[c(2:501)])+1)]

plot(hgram, xlim = c(0.9,1.1))
abline(v = hgram$mids[(which.max(hgram$counts[c(2:501)])+1)], col = 'red')




pqNorm <- function(X, pcent = T){browser()
  Xm <- apply(X, 2, median)
  dilf <- sapply(1:nrow(X), function(x){
    hgram <- hist((X[x,]/Xm))
    vic <- c(min(hgram$breaks), seq(from = -0.005, to = 5.005, length.out = 502), max(hgram$breaks))
    h <- hist((X[x,]/Xm), breaks = vic)
    d <- h$mids[(which.max(h$counts[c(2:501)])+1)]
    return(d)
  })
  Xn <- sapply(1:nrow(X), function(y){
    if(dilf[y] == 1){
      X[y,]
    } else {
      X[y,]/dilf[y]
    }
  })
  if (pcent){
    p_cent <- sapply(1:nrow(X), function(x){
      hgram <- hist((X[x,]/Xm))
      vic <- c(min(hgram$breaks), seq(from = -0.005, to = 5.005, length.out = 502), max(hgram$breaks))
      h <- hist((X[x,]/Xm), breaks = vic)
      p_cent <- ((sum(h$counts[c(2:501)]))/(sum(h$counts)))*100
      return(p_cent)
    })
    return(list(Xn, dilf, p_cent))
  }else{
    return(list(Xn, dilf))
  }
}


Xpqn = pqNorm(X, pcent = T)

pqn_dilf <- Xpqn[[2]]
pqn_pcent <- Xpqn[[3]]
Xpqn <- t(Xpqn[[1]])
Xm <- apply(X, 2, median)



x1 <- Xpqn[1,idx]/Xm[idx]

idx <- get.idx(c(1,4), ppm)





pqNorm <- function(X, pcent = F){
  Xm <- apply(X, 2, median)
  Xm[Xm <= 0]=1
  dilf <- sapply(1:nrow(X), function(x){
    Xt <- X[x,]
    Xt[Xt<=0]=1
    quo <- (Xt/Xm)
    hgram <- hist(quo)
    vic <- c(0, seq(from = 0.025, to = 5.025, length.out = 501), max(hgram$breaks))
    h <- hist((Xt/Xm), breaks = vic, xlim = c(0,5))
    d <- h$mids[(which.max(h$counts[(c(2:501))]))]
    return(d)
  })

  Xn <- sapply(1:nrow(X), function(y){
      X[y,]/dilf[y]
  })
  return(list(Xn, dilf))
}



Xpqn5 <- pqNorm(X1)
















d <-sapply(1:nrow(X), function(x){
  density((X[x,]/Xm), from = 0, to = 5)
})

sapply(1:ncol(d), function(x){
  plot(d)
  abline(v = d[1,x[which.max(d$y)])
})
d$x[which.max(d$y)]
plot(d)
abline(v = d$x[which.max(d$y)])













