##############################################################
#periodic DEMMLER-REINSCH BASIS FUNCTIONS
##############################################################

# Euler-Frobenius Polynomial
#
# param p degree of polynomial
# param t r-dim. vector of evaluation points
# return res r-dim. vector with Euler-Frobenius poly. of degree p at t, for details see (Schwarz, Krivobokova; 2016).
# keyword Internal
EF.poly <-function(p,t){
  if (p<0){
    stop ("p must be non-negative")
  }
  res=switch(p+1, #+1 to get cases 0,1,... instead of 1,2,...
             t^0, #p=0
             t^0, #p=1
             t+1, #p=2
             t^2+4*t+1, #p=3
             t^3+11*t^2+11*t+1, #p=4
             t^4+26*t^3+66*t^2+26*t+1) #p=5
  if (p>5){
    res=c()
    for (te in t){
      sum=0
      for (i in 0:p){
        j=seq(0,p-i)
        sum=sum+(choose(p+1,i)*(-1)^i*sum(j^p*te^(p-i-j)))
      }
      res=c(res,sum)
    }
  }
  return(res)
}


# Q_p-Polynomials
#
# param pminus1 degree of polynomial
# param x r-dim. vector of evaluation points
# return res r-dim. vector with Q polynomial of degree pminus1 at x, for details see (Schwarz, Krivobokova; 2016).
# keyword Internal
Q_pminus1 <- function(pminus1,x){
  if (pminus1<0){
    stop ("pminus1 must be non-negative")
  }
  res=switch(pminus1+1, #+1 to get cases 0,1,... instead of 1,2,...
             x^0, #Q_(p-1)=Q_0
             1/2+cos(pi*x)^2/2,        #Q_(p-1)=Q1
             1/3+2*cos(pi*x)^2/3,      #Q_(p-1)=Q2
             5/24+3*cos(pi*x)^2/4+cos(pi*x)^4/24, #Q_(p-1)=Q3
             2/15+11*cos(pi*x)^2/15+2*cos(pi*x)^4/15, #Q_(p-1)=Q4
             (61+479*cos(pi*x)^2+179*cos(pi*x)^4+cos(pi*x)^6)/720, #Q_(p-1)=Q5
             17/315+4*cos(pi*x)^2/7+38*cos(pi*x)^4/105+4*cos(pi*x)^6/315) #Q_(p-1)=Q6
  return (res)
}

#INPUT
##x:
##N:
##q:  of smoothing spline
#OUTPUT
##basis:
#' Periodic Demmler-Reinsch Basis
#'
#' Calculates the periodic Demmler-Reinsch basisfor a given smoothness and a given vector of grid points. For details see (Schwarz, Krivobokova; 2016).
#' @param x \code{m}-dim. vector with grid values in \[0,1\]
#' @param n dimension of the basis
#' @param q penalization order, \code{q=1,2,3,4} are available
#' @return \code{mxn} dimensional matrix with the \code{n} DR basis functions evaluated at grid points \code{x}
#' @export
#' @examples DR.basis(seq(1,10)/10,5,2)
DR.basis<-function(x,n,q){
  if (!(q %in% c(1,2,3,4))){
    stop("q can only attain the vaules 1,2,3,4")
  }
  p=2*q-1
  t=n*x+q
  Q.M=Q_pminus1(p-1,(1:n)/n)

  #build exponential spline
  j=seq(0,p)
  z=exp(-2*pi*1i*(1:n)/n)
  basis=c()
  for (te in t){
    res=c()
    t.int=floor(te)
    idx=0
    for(ze in z){
      idx=idx+1
      if(ze==1){
        res=c(res,1/Q.M[idx])
      }else{
        res=c(res,ze^t.int*(1-1/ze)^p*sum(choose(p,j)*(te-t.int)^(p-j)*sapply(j,EF.poly,ze)/(factorial(p)*(ze-1)^j))/Q.M[idx])
      }
    }
    basis=matrix(rbind(basis,res),ncol=n)
  }
  # if(all(abs(Im(basis))<10^(-8))){
  #   basis=Re(basis)}
  return(basis)
}


