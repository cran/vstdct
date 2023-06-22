#library(matrixcalc) #check if needed


##############################################################
#1.DATA TRANSFORMATION FUNCTIONS
##############################################################
# Data Binning
#
# Each \code{p}-dimensional vector is divided into \code{Te} consecutive bins. The data points in each bin and across the \code{n} samples are summed.
# param w \code{nxp} dimensional data matrix
# param Te number of bins for data binning. \code{Te} should be smaller than the vector length \code{p}.
# param p vector length
# param n sample size
# return A list containing the following elements:
#         \itemize{
#         \item{\code{m}:  }{number of data points per bin, that is \code{m=n*round(p/Te)}}
#          \item{\code{w.star}:   }{\code{Te}-dimensional vector with \code{m} summed data points for each bin.  The bin number \code{Te} may be modified to guarantee at least two data points per bin. If \code{p/Te} is not an integer, the vector dimension can be different than \code{Te} and in the last bin more than \code{m} data points may be summed.}
#          }

Binning<-function(w,Te,p,n){
  if (Te>p){
    stop("Te must be smaller or equal than p.")
  }
  if(sum(dim(w) == c(n,p))!=2){
    stop("dimension of w does not match c(n,p)")
  }
  if(p/Te == 1 && n>1){
    #"n>1 case"
    m=1
    w.star=colSums(w)
  }else{
    #"time series case"
    m=max(2,round(p/Te)) #m=max(2,floor(p/Te))
    start=seq(1,p,by=m)
    start=start[start<=(p-m+1)] #last bin may contains #obs>=m
    end=c(start[-1]-1,p)
    w.star=mapply(function(start,end){sum(w[,start:end])},start,end)
  }
  return(list(m=m*n,w.star=w.star))
}

# Variance Stabilizing Transform
#
# Applies the variance stabilizing transform for gamma distribution to a vector. For details see (Brown et al., \emph{"Nonparametric regression in exponential families."}, The Annals of Statistics 38.4 (2010): 2005-2046).
# param w.star vector with binned data points
# param m number of summed data points per bin
# return vector with binned and variance stabilized data, i.e. \code{log(./m)/sqrt(2)} is applied entrywise to \code{w.star}.

VST<-function(w.star,m){
  y.star=log(w.star/m)/sqrt(2)
  return(y.star)
}


#' Data Transformation
#'
#' Applies the Discrete Cosine I transform, data binning and the variance stabilizing transform function to the data.
#' @param y \code{nxp} dimensional data matrix
#' @param Te number of bins for data binning.  \code{Te} should be smaller than the vector length \code{p}.
#' @param dct.out logical. If \code{TRUE}, the \code{p}-dim. DCT-I matrix is returned. The default is \code{FALSE}.
#' @return A list containing the following elements:
#'         \itemize{
#'         \item{\code{m}:  }{number of data points per bin, that is \code{m=n*round(p/Te)}. If \code{p/Te} is not an integer, the first/last bin may contain more than \code{m} data points.}
#'         \item{\code{y.star}:   }{\code{2Te-2} dimensional vector with binned, variance stabilized and mirrowed data. The bin number \code{Te} may be modified to guarantee at least two data points per bin. If \code{p/Te} is not an integer, the vector dimension is \code{2*floor(p/round(p/Te))-2}.}
#'         \item{\code{dct.matrix}: }{\code{p}-dim. DCT-I matrix (if \code{dct.out}=TRUE)}
#'         }
#' @export
Data.trafo<-function(y,Te,dct.out=FALSE){
  n=dim(y)[1]
  p=dim(y)[2]
  if (Te>p){
    stop("Te must be smaller or equal than p")
  }
  if (is.null(n)){
    p=length(y)
    n=1
  }
  dct.matrix=DCT.matrix(p)
  w=(y%*%dct.matrix)^2
  Bins=Binning(w,Te,p,n)
  y.star=VST(Bins$w.star,Bins$m)
  TTe=length(y.star)
  y.star=c(y.star[2:(TTe-1)],y.star[TTe:1])
  if(dct.out){
    return(list(m=Bins$m,y.star=y.star,dct.matrix=dct.matrix))
  }else{
    return(list(m=Bins$m,y.star=y.star))
  }
}



##############################################################
#2.SPECTRAL DENSITY ESTIMATION
##############################################################
# Periodic Smoothing Spline
#
# computes the smoothing spline estimator of order q using the periodic Demmler-Reinsch basis.
# param x vector with values in \[0,1\] at which the periodic smoothing spline is evaluated
# param y vector with data points on equi-spaced design points
# param lambda smoothing parameter
# param q penalization order
# return f periodic smoothing spline estimator evaluated at x
# keyword Internal
SSper.estimator<-function(x,y,lambda,q){
  if (!(q %in% c(1,2,3,4))){
    stop("q can only attain the vaules 1,2,3,4")
  }
  if (lambda<=0){
    stop("lambda must be positive")
  }
  N=length(y)
  x.equi=seq(1,N)/N
  evals=(2*pi*1:N)^(2*q)*sinc(pi*x.equi)^(2*q)/Q_pminus1(2*q-2,x.equi)
  f=stats::fft(y[N:1])[N:1]/(1+lambda*evals)
  if(identical(x,seq(1,N)/N)){
    f=Re(stats::fft(f[N:1])[N:1])/N
  }else{
    f=c(Re(DR.basis(x,N,q)%*%f/N)) #imaginary part zero by construction
  }
  return(list(f=f,evals=evals))
}

# Spectral Density Estimator
#
# Estimates the spectral density from the binned and transformed data points. A periodic smoothing spline is fitted for the log-spectral density \eqn{H(f)}, where \eqn{ H(y)= \left \{ \phi(m/2) + \log \left ( 2y/m\right ) \right\}/\sqrt{2}} and \eqn{\phi} is the digamma function.
# The spectral density is then obtained by applying the inverse variance stabilizing transform  \eqn{H^{-1}(y)= m\exp \left \{\sqrt{2}y-\phi\left (m/2\right) \right\}/2}.
# param x vector with values in \[0,1\] at which the spectral density estimator is evaluated
# param y.star vector with binned, variance stabilized and mirrowed data
# param lambda smoothing parameter
# param q penalization order
# param m number of data points per bin
# return A list containing the following elements:
#         \itemize{
#          \item{\code{f}:   }{vector with the spectral density estimator evaluated at \code{x}}
#          \item{\code{Hf}:  }{vector with the log spectral density \eqn{H(f)} evaluated at \code{x}}
#          }
sdf.estimator<-function(x,y.star,lambda,q,m){
  if (!(q %in% c(1,2,3,4))){
    stop("q can only attain the vaules 1,2,3,4")
  }
  if (lambda<=0){
    stop("lambda must be positive")
  }
    Hf=SSper.estimator(x,y.star,lambda,q)$f
    f=m/2*exp(sqrt(2)*Hf-digamma(m/2))
    return(list(f=f,Hf=Hf))
}

##############################################################
#3.TOEPLITZ MATRIX ESTIMATION
##############################################################
#' Toeplitz Covariance and Precision Matrix Estimator
#'
#' Estimates the Toeplitz covariance matrix, the inverse matrix and the spectral density from a sample of \code{n} i.i.d. \code{p}-dimensional vectors with mean zero.
#' @param y \code{nxp}  dimensional data matrix
#' @param Te number of bins for data binning.
#' @param q penalization order, \code{q=1,2,3,4} are available
#' @param method to select the smoothing parameter of the smoothing spline. Available methods are restricted maxmimum likelihodd "\code{ML}", generalized cross-validation "\code{GCV}" and the oracle versions "\code{ML-oracle}", "\code{GCV-oracle}".
#' @param f.true Te-dimensional vector with the true spectral density function evaluated at equi-sapced points in \[0,\code{pi}\]. Only required, if an oracle method ("\code{ML-oracle}", "\code{GCV-oracle}") is chosen for \code{method}.
#' @return A list containing the following elements:
#'         \itemize{
#'          \item{\code{toep}:   }{\code{p}-dim. Toeplitz covariance matrix}
#'          \item{\code{toep.inv}:  }{\code{p}-dim. precision matrix}
#'          \item{\code{acf}:  }{\code{p}-dim. vector with the covariance function}
#'          \item{\code{sdf}:  }{\code{p}-dim. vector with the spectral density in the interval \[0,1\]}
#'          }
#' @export
#' @examples
#' #EXAMPLE 1: Simulate Gaussian ARMA(2,2)
#' library(nlme)
#' library(MASS)
#' p=100
#' n=1
#' Sigma=1.44*corMatrix(Initialize(corARMA(c(0.7, -0.4,-0.2, 0.2),p=2,q=2),data=diag(1:p)))
#' Y=matrix(mvrnorm(n, mu=numeric(p), Sigma=Sigma),n,p)
#' fit.toep=Toep.estimator(y=Y,Te=10,q=2,method="GCV")$toep
#'
#'
#' #EXAMPLE 2: AQUAPORIN DATA
#' data(aquaporin)
#' n=length(aquaporin$Y)
#' y.train=aquaporin$Y[1:(0.01*n)]
#' y.train=y.train-mean(y.train)
#' fit.toep=Toep.estimator(y=y.train,Te=10,q=1,method="ML")$toep

Toep.estimator<-function(y,Te,q,method,f.true=NULL){
  n=dim(y)[1]
  p=dim(y)[2]
  if (is.null(n)){
    p=length(y)
    n=1
    y=matrix(y,n,p)
  }
  if (!(q %in% c(1,2,3,4))){
    stop("q can only attain the vaules 1,2,3,4")
  }
  if (Te>p){
    stop("Te must be smaller or equal than p")
  }
  if (! method %in% c("GCV","GCV-oracle","ML","ML-oracle")){
    stop("method is invalid, choose one of GCV, GCV-oracle, ML, ML-oracle")
  }
  Duplicate<-function(x){
    r=length(x)
    return(c(x[2:(r-1)],x[r:1]))
  }
  Data=Data.trafo(y,Te)
  if(method=="GCV")
    lambda=GCV.opt(Data$y.star,q)
  if(method=="GCV-oracle"){
    if (is.null(f.true)){
      stop("true spectral density is missing")
    }
    f.true=Duplicate(f.true)
    log.sdf.true=VST(f.true,m=1)
    lambda=GCV.opt(y=NULL,q,oracle=TRUE,log.sdf.true,sigma=1/Data$m)
  }
  if(method=="ML")
    lambda=ML.opt(Data$y.star,q)
  if(method=="ML-oracle"){
    if (is.null(f.true)){
      stop("true spectral density is missing")
    }
    f.true=Duplicate(f.true)
    log.sdf.true=VST(f.true,m=1)
    lambda=ML.opt(y=NULL,q,oracle=TRUE,log.sdf.true,sigma=1/Data$m)
  }
  xx=seq(0,0.5,length=p) #due to symmetry: evaluation on [0,0.5] sufficient
    sdf.fit=sdf.estimator(xx,Data$y.star,lambda,q,Data$m)$f
    acf.fit=sdf2acf(sdf.fit) #DFT approximation to integral
    #acf=integrate(f.estimator,0,1,Data$y.star,lambda,q,Data$m)[[1]] #R approximation to integral
    acf.fit.inv=sdf2acf(1/sdf.fit) #DFT approximation to integral

  return(list(toep=stats::toeplitz(acf.fit),toep.inv=stats::toeplitz(acf.fit.inv),acf=acf.fit,sdf=sdf.fit))
}

#################################################################
#4. SMOOTHING PARAMETER SELECTION
#################################################################
# Generalized Cross-Validation
#
# param y  N-dim. vector with regression data points
# param q  penalization order
# param oracle logical. If TRUE, then the GCV-oracle smoothing parameter is calculated.
# param f.true N-dim. vector with the true regression function evaluated at equi-spaced points. Only required if oracle=TRUE.
# param sigma  variance of the regression data
# return lambda estimated smoothing parameter
# keyword Internal
GCV.opt<-function(y,q,oracle=FALSE,f.true=NULL,sigma=NULL){

  GCV<-function(lambda,y,q){
    N=length(y)
    x=seq(1,N)/N
    evals=(2*pi*1:N)^(2*q)*sinc(pi*x)^(2*q)/Q_pminus1(2*q-2,x)
    term1=stats::fft(y[N:1])[N:1]*(lambda*evals/(1+lambda*evals))^2
    term1=c(y%*%Re(stats::fft(term1[N:1])[N:1])/N) #inverse=TRUE not needed for symmetric vector
    trace=sum(1/(1+lambda*evals))
    score=term1/(N*(1-trace/N)^2)
    return(score)
  }

  GCV.oracle<-function(lambda,f.true,q,sigma){
    N=length(f.true)
    x=seq(1,N)/N
    evals=(2*pi*1:N)^(2*q)*sinc(pi*x)^(2*q)/Q_pminus1(2*q-2,x)
    norm.term=stats::fft(f.true[N:1])[N:1]*lambda*evals/(1+lambda*evals)
    norm.term=Re(stats::optimize(norm.term[N:1])[N:1])/N
    trace=sum(1/(1+lambda*evals)^2)
    score=(sum(norm.term^2)+trace*sigma)/N
    return(score)
  }

  if(oracle)
    lambda=stats::optimize(GCV.oracle,interval=c(.Machine$double.eps,1),f.true=f.true,q=q,sigma=sigma,tol=.Machine$double.eps)$minimum #1
  else
    lambda=stats::optimize(GCV,interval=c(.Machine$double.eps,1),y=y,q=q,tol=.Machine$double.eps)$minimum #1
    #lambda=uniroot(GCV2,interval=c(.Machine$double.eps,1),y=y,q=q,tol=.Machine$double.eps)$root
  return(lambda)
}

# Maximum Likelihood method
#
# param y  N-dim. vector with regression data points
#param q  penalization order
# param oracle logical. If TRUE, then the ML-oracle smoothing parameter is calculated.
# param f.true N-dim. vector with the true regression function evaluated at equi-spaced points. Only required if oracle=TRUE.
# param sigma  variance of the regression data
# return estimated smoothing parameter
# keyword Internal
ML.opt<-function(y,q,oracle=FALSE,f.true=NULL,sigma=NULL){

  ML<-function(lambda,y,q){
    N=length(y)
    x=seq(1,N)/N
    evals=(2*pi*1:N)^(2*q)*sinc(pi*x)^(2*q)/Q_pminus1(2*q-2,x)
    term1=stats::fft(y[N:1])[N:1]*(1/(1+lambda*evals)-1/(1+lambda*evals)^2)
    term1=c(y%*%Re(stats::fft(term1[N:1])[N:1])/N) #inverse=TRUE not needed for symmetric vector
    fac=(sum(1/(1+lambda*evals))-q)/(N-q)
    term2=stats::fft(y[N:1])[N:1]*(lambda*evals/(1+lambda*evals))
    term2=c(y%*%Re(stats::fft(term2[N:1])[N:1])/N) #inverse=TRUE not needed for symmetric vector
    score=(term1-term2*fac)/N
    return(score)
  }

  ML.oracle<-function(lambda,q,f.true,sigma){
    N=length(f.true)
    x=seq(1,N)/N
    evals=(2*pi*1:N)^(2*q)*sinc(pi*x)^(2*q)/Q_pminus1(2*q-2,x)
    term1=stats::fft(f.true[N:1])[N:1]*(1/(1+lambda*evals)-1/(1+lambda*evals)^2)
    term1=c(f.true%*%Re(stats::fft(term1[N:1])[N:1])/N) #inverse=TRUE not needed for symmetric vector
    term2=sum(1/(1+lambda*evals)^2)
    score=(term1-sigma*(term2-q))/N
    return(score)
  }
  if(oracle)
    lambda=stats::uniroot(ML.oracle,interval=c(.Machine$double.eps,1),f.true=f.true,q=q,sigma=sigma,tol=.Machine$double.eps, extendInt = "yes")$root #1
  else
    lambda=stats::uniroot(ML,interval=c(.Machine$double.eps,1),y=y,q=q,tol=.Machine$double.eps, extendInt = "yes")$root #1
  return(lambda)
}


