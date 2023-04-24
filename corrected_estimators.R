#testing bias correction
# number of ind in pool
n = 1000
# read coverage
M = 1000
# minimal read coverage
b = 2

# correction for theta

M= 100

b = 2

# theta_py inner sum
tp_inner <- function(m,M,n) {
  # takes running index of m, coverage M and pool size n
  r=seq(1,n-1,by=1)
  PXm_Yr=dbinom(m,M,r/n)
  res = sum(m*(M-m)*PXm_Yr/r)
  return(res)
}


tw_corr <- function(b,M,n,sites=1){
  # takes running minimal count of b, coverage M and pool size n and optional number of sites, returns the bias corrected wattersons theta
  sites=1
  cm = sum(1/seq(1,M-1,by=1))
  r=seq(1,n-1,by=1)
  return(sites / sum( (pbinom(M-b,M,r/n)-pbinom(b-1,M,r/n))/r ) )
}

tw_corr2 <- function(b,M,n, sites=1){
  # takes running minimal count of b, coverage M and pool size n and optional number of sites
  # returns the bias corrected wattersons theta for high n compared to M
  cm = sum(1/seq(1,M-1,by=1))
  cb = sum(1/seq(b,M-b,by=1))
  r=seq(1,n-1,by=1)
  return( sites/cb)
}



n=50
M=150
b=5

theta_corr = choose(M,2) / sum(vapply(seq(b,M-b,by=1), tp_inner, 1, M = M, n = n))


(M-1)/(M-2*b+1)


tw_corr(b,M,n)

tw_corr2(b,M,n)


