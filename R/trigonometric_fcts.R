# Discrete Cosine Transform matrix
#
# Computes the Discrete Cosine Transform I (DCT-I) matrix.
# param N matrix dimension
# return the \code{NxN}-dimensional DCT-I matrix
DCT.matrix<-function(N){
  seq=seq(0,N-1)
  arg=(pi*seq/(N-1))
  DCT=sapply(seq, function(j){cos(arg*j)})
  fac=matrix(1/sqrt(2),N,N)
  fac[2:(N-1), 2:(N-1)]=1
  fac[1,1]=fac[1,N]=fac[N,1]=fac[N,N]=1/2
  DCT=DCT*fac*sqrt(2)/sqrt(N-1)
  return(DCT)
}

# Sinc Function
#
# Input: x vector
# Output: y sinc(x)
sinc <- function(x){
  y=sin(x)/x
  y[is.na(y)]=1
  return(y)
}

# Auto-covariance to Spectral Density Function
#
# Input: x vector with the auto-covariance function
# Output: y  vector with the spectral density function
acf2sdf <- function(x){
  N=length(x)
  y=2*dtt::dtt(x, type="dct", variant = 1) #*sqrt(2/(N-1))
  #y=Re(fft(c(x,x[(N-1):2]))) #the same
  return (y)
}

# Spectral Density to Auto-covariance Function
#
# Input:  x vector with the spectral density function
# Output:  y  vector with the auto-covariance function
sdf2acf <- function(x){
  N=length(x)
  y=dtt::dtt(x, type="dct", variant = 1)/(N-1)
  return (y)
}
