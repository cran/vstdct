#library(nlme) check if needed
#library(MASS) check if needed

#' Data Examples
#'
#' [vstdct::example1], [vstdct::example2] and [vstdct::example3] generate i.i.d. vectors from a given distribution with different Toeplitz covariance matrices.
#' The covariance function \eqn{\sigma} of the Toeplitz covariance matrix of
#'         \itemize{
#'         \item{`example1`: }{has a polynomial decay, \eqn{\sigma(\tau)= sd^2(1+|\tau|)^{-gamma}},}
#'         \item{`example2`: }{follows an \eqn{ARMA(2,2)} model with coefficients \eqn{(0.7,-0.4,-0.2,0.2)} and innovations variance \eqn{sd^2},}
#'         \item{`example3`: }{yields a Lipschitz continuous spectral density \eqn{f} that is not differentiable, i.e. \eqn{f(x)= sd^2({|\sin(x+0.5\pi)|^{gamma}+0.45})}}
#'         }

#' @name Data Examples
#' @param p  vector length
#' @param n  sample size
#' @param sd standard deviation
#' @param gamma polynomial decay of covariance function for `example1` resp. exponent for `example3`
#' @param family distribution of the simulated data. Available distributions are "`Gaussian`", "`Gamma`", "`Uniform`". The default is "`Gaussian`".
#'
#' @return A list containing the following elements:
#'         \itemize{
#'         \item{`Y`:  }{`pxn` dimensional data matrix}
#'         \item{`sdf`:   }{true spectral density function}
#'         \item{`acf`:   }{true covariance function}
#'         }
#' @export
#' @examples example1(p=10, n=1, sd=1, gamma=1.2, family="Gaussian")
example1<-function(p, n, sd, gamma,family="Gaussian"){
  x=seq(0,1,length=p)
  acf=c(sd^2,sd^2/(1+seq(1,p-1))^gamma)
  Sigma=stats::toeplitz(acf)
  sdf=acf2sdf(acf)
  if (family=="Gaussian")
    Y=matrix(MASS::mvrnorm(n, mu=numeric(p), Sigma=Sigma),n,p)
  if (family=="Gamma"){
    Y=matrix(NA,n,p)
    Sigma.sq2=t(chol(Sigma))
    for (i in 1:n){
      Y[i,]=Sigma.sq2%*%(stats::rgamma(n=p,shape=2,scale=1/sqrt(2))-sqrt(2))
    }
  }
  if (family=="Uniform"){
    Y=matrix(NA,n,p)
    Sigma.sq2=t(chol(Sigma))
    for (i in 1:n){
      Y[i,]=Sigma.sq2%*%(stats::runif(n=p,min=-sqrt(3),max=sqrt(3)))
    }
  }
  return(list(Y=Y,sdf=sdf,acf=acf))
}

#'
#' example3(p, n, sd, gamma,family="Gaussian")
#' [vstdct::example1], [vstdct::example2] and [vstdct::example3] generate i.i.d. vectors from a given distribution with different Toeplitz covariance matrices.
#' The covariance function \eqn{\sigma} of the Toeplitz covariance matrix of
#'         \itemize{
#'         \item{`example1`: }{has a polynomial decay, \eqn{\sigma(\tau)= sd^2(1+|\tau|)^{-gamma}},}
#'         \item{`example2`: }{follows an \eqn{ARMA(2,2)} model with coefficients \eqn{(0.7,-0.4,-0.2,0.2)} and innovations variance \eqn{sd^2},}
#'         \item{`example3`: }{yields a Lipschitz continuous spectral density \eqn{f} that is not differentiable, i.e. \eqn{f(x)= sd^2{|\sin(x+0.5\pi)|^{gamma}+0.45}}}
#'         }

#' @name Data Examples
#' @param p  vector length
#' @param n  sample size
#' @param sd standard deviation
#' @param gamma polynomial decay of covariance function for `example1` resp. exponent for `example3`
#' @param family distribution of the simulated data. Available distributions are "`Gaussian`", "`Gamma`", "`Uniform`". The default is "`Gaussian`".
#'
#' @export
#' @examples example2(p=10,n=1,sd=1,family="Gaussian")
example2<-function(p,n,sd,family="Gaussian"){
  x=seq(0,1,length=p)
  pp=2
  qq=2
  coef=c(0.7, -0.4,-0.2, 0.2)
  R.mat=nlme::corMatrix(nlme::Initialize(nlme::corARMA(coef,p=pp,q=qq),data=diag(1:p)))
  Sigma=sd^2*R.mat
  acf=Sigma[1,]
  sdf=acf2sdf(acf)
  if (family=="Gaussian")
    Y=matrix(MASS::mvrnorm(n, mu=numeric(p), Sigma=Sigma),n,p)
  if (family=="Gamma"){
    Y=matrix(NA,n,p)
    for (i in 1:n){
      Y[i,]=sqrt(sd)*stats::arima.sim(n = p, list(ar = c(0.7, -0.4), ma = c(-0.2, 0.2)), innov=stats::rgamma(n=p,shape=2,scale=1/sqrt(2))-sqrt(2))
    }
  }
  if (family=="Uniform"){
    Y=matrix(NA,n,p)
    for (i in 1:n){
      Y[i,]=sqrt(sd)*stats::arima.sim(n = p, list(ar = c(0.7, -0.4), ma = c(-0.2, 0.2)), innov=stats::runif(p,min=-sqrt(3), max=sqrt(3)))
    }
  }
  return(list(Y=Y,sdf=sdf,acf=acf))
}



#' [vstdct::example1], [vstdct::example2] and [vstdct::example3] generate i.i.d. vectors from a given distribution with different Toeplitz covariance matrices.
#' The covariance function \eqn{\sigma} of the Toeplitz covariance matrix of
#'         \itemize{
#'         \item{`example1`: }{has a polynomial decay, \eqn{\sigma(\tau)= sd^2(1+|\tau|)^{-gamma}},}
#'         \item{`example2`: }{follows an \eqn{ARMA(2,2)} model with coefficients \eqn{(0.7,-0.4,-0.2,0.2)} and innovations variance \eqn{sd^2},}
#'         \item{`example3`: }{yields a Lipschitz continuous spectral density \eqn{f} that is not differentiable, i.e. \eqn{f(x)= sd^2{|\sin(x+0.5\pi)|^{gamma}+0.45}}}
#'         }

#' @name Data Examples
#' @param p  vector length
#' @param n  sample size
#' @param sd standard deviation
#' @param gamma polynomial decay of covariance function for `example1` resp. exponent for `example3`
#' @param family distribution of the simulated data. Available distributions are "`Gaussian`", "`Gamma`", "`Uniform`". The default is "`Gaussian`".
#'

#' @export
#' @examples example3(p=10, n=1, sd=1, gamma=2,family="Gaussian")
example3<-function(p, n, sd, gamma,family="Gaussian"){
  x=seq(0,1,length=p)
  sdf=sd^2*(abs(sin(pi*x +0.5*pi))^gamma+0.45)
  acf=sdf2acf(sdf)
  Sigma=stats::toeplitz(acf)
  if (family=="Gaussian")
    Y=matrix(MASS::mvrnorm(n, mu=numeric(p), Sigma=Sigma),n,p)
  if (family=="Gamma"){
    Y=matrix(NA,n,p)
    Sigma.sq2=t(chol(Sigma))
    for (i in 1:n){
      Y[i,]=Sigma.sq2%*%(stats::rgamma(n=p,shape=2,scale=1/sqrt(2))-sqrt(2))
    }
  }
  if (family=="Uniform"){
    Y=matrix(NA,n,p)
    Sigma.sq2=t(chol(Sigma))
    for (i in 1:n){
      Y[i,]=Sigma.sq2%*%(stats::runif(n=p,min=-sqrt(3),max=sqrt(3)))
    }
  }
  return(list(Y=Y,sdf=sdf,acf=acf))
}
