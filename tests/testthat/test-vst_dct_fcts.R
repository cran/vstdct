test_that("valid arguments are passed - Binning", {
  n=100
  p=200
  Te=201
  w=matrix(1,n,p)
  expect_error(Binning(w=w,Te=Te,p=p,n=n))
  Te=20
  expect_error(Binning(w=t(w),Te=Te,p=p,n=n))
})


test_that("output type and length are correct - Binning", {
  n=100
  p=200
  Te=4
  w=matrix(1,n,p)
  expect_type(Binning(w=w,Te=Te,p=p,n=n), "list")
  expect_length(Binning(w=w,Te=Te,p=p,n=n), 2)
})

test_that("number of points per bin is correct - Binning", {
  #Case: n>1 and p/Te integer > 1
  n=100
  p=200
  Te=10
  w=matrix(1,n,p)
  expect_equal(Binning(w=w,Te=Te,p=p,n=n)$m, max(2,round(p/Te))*n)
  #Case: n>1 and p/Te not an integer
  n=100
  p=200
  Te=11
  w=matrix(1,n,p)
  expect_equal(Binning(w=w,Te=Te,p=p,n=n)$m,  max(2,round(p/Te))*n)
  #Case: n>1 and p/Te = 1
  n=100
  p=200
  Te=p
  w=matrix(1,n,p)
  expect_equal(Binning(w=w,Te=Te,p=p,n=n)$m,  n)
  #Case: n=1 and p/Te not an integer
  n=1
  p=200
  Te=6
  w=matrix(1,n,p)
  expect_equal(Binning(w=w,Te=Te,p=p,n=n)$m,  max(2,round(p/Te))*n)
})

test_that("length of binned vector is correct - Binning", {
  #Case: n>1 and p/Te integer > 1
  n=100
  p=200
  Te=10
  w=matrix(1,n,p)
  expect_length(Binning(w=w,Te=Te,p=p,n=n)$w.star,Te)
  #Case: n>1 and p/Te not an integer
  n=100
  p=200
  Te=11
  w=matrix(1,n,p)
  mm=max(2,round(p/Te))
  expect_length(Binning(w=w,Te=Te,p=p,n=n)$w.star,floor(p/mm))
  #Case: n>1 and p/Te = 1
  n=100
  p=200
  Te=p
  w=matrix(1,n,p)
  expect_length(Binning(w=w,Te=Te,p=p,n=n)$w.star,Te)
  #Case: n=1 and p/Te not an integer
  n=1
  p=200
  Te=6
  w=matrix(1,n,p)
  mm=max(2,round(p/Te))
  expect_length(Binning(w=w,Te=Te,p=p,n=n)$w.star,floor(p/mm))
})

test_that("output is correct - VST", {
  w.star=seq(15)
  m=10
  expect_equal(VST(w.star=w.star,m=m),log(w.star/m)/sqrt(2))
})

test_that("valid arguments are passed - Data.trafo", {
  n=100
  p=200
  Te=p+1
  y=matrix(1,n,p)
  expect_error(Data.trafo(y,Te))
})

test_that("output type and length are correct - Data.trafo", {
  n=100
  p=200
  Te=10
  y=matrix(1,n,p)
  expect_type(Data.trafo(y=y,Te=Te), "list")
  expect_length(Data.trafo(y=y,Te=Te), 2)
})

test_that("output is correct - Data.trafo",{
  #Case: n>1 and p/Te integer > 1
  n=100
  p=200
  Te=10
  y=matrix(1,n,p)
  expect_length(Data.trafo(y=y,Te=Te)$y.star, 2*Te-2)
  expect_equal(Data.trafo(y=y,Te=Te)$m,max(2,round(p/Te))*n)
  #Case: n>1 and p/Te not an integer
  n=100
  p=200
  Te=11
  y=matrix(1,n,p)
  mm=max(2,round(p/Te))
  Tee=floor(p/mm)
  expect_length(Data.trafo(y=y,Te=Te)$y.star, 2*Tee-2)
  expect_equal(Data.trafo(y=y,Te=Te)$m,max(2,round(p/Te))*n)
  #Case: n>1 and p/Te = 1
  n=100
  p=200
  Te=p
  y=matrix(1,n,p)
  expect_length(Data.trafo(y=y,Te=Te)$y.star,2*Te-2)
  expect_equal(Data.trafo(y=y,Te=Te)$m,n)
  #Case: n=1 and p/Te not an integer
  n=1
  p=200
  Te=6
  y=matrix(1,n,p)
  mm=max(2,round(p/Te))
  Tee=floor(p/mm)
  expect_length(Data.trafo(y=y,Te=Te)$y.star, 2*Tee-2)
  expect_equal(Data.trafo(y=y,Te=Te)$m,max(2,round(p/Te))*n)
})

test_that("valid arguments are passed - SSper.estimator",{
  n=100
  x=seq(1,n)/n
  y.star=seq(1,100)
  expect_error(SSper.estimator(x,y.star,lambda=1,q=5))
  expect_error(SSper.estimator(x,y.star,lambda=-1,q=2))
})

test_that("output type and length are correct - SSper.estimator",{
  Te=100
  x=seq(1,Te)/Te
  y.star=seq(1,100)
  q=3
  m=10
  expect_type(SSper.estimator(x,y.star,lambda=1,q=q), "list")
  expect_length(SSper.estimator(x,y.star,lambda=1,q=q), 2)
  expect_length(SSper.estimator(x,y.star,lambda=1,q=q)$f, Te)
  expect_length(SSper.estimator(x,y.star,lambda=1,q=q)$evals, Te)
})

test_that("valid arguments are passed - sdf.estimator",{
  n=100
  x=seq(1,n)/n
  y.star=seq(1,100)
  m=10
  expect_error(sdf.estimator(x,y.star,lambda=1,q=5,m))
  expect_error(sdf.estimator(x,y.star,lambda=-1,q=2,m))
})

test_that("output type and length are correct - sdf.estimator",{
  Te=100
  x=seq(1,Te)/Te
  y.star=seq(1,100)
  q=3
  m=10
  expect_type(sdf.estimator(x,y.star,lambda=1,q=q,m), "list")
  expect_length(sdf.estimator(x,y.star,lambda=1,q=q,m), 2)
  expect_length(sdf.estimator(x,y.star,lambda=1,q=q,m)$f, Te)
  expect_length(sdf.estimator(x,y.star,lambda=1,q=q,m)$Hf, Te)
})

test_that("valid arguments are passed - Toep.estimator",{
  n=100
  p=200
  Te=10
  set.seed(11)
  acf=c(1.44,1.44/(1+seq(1,p-1))^2.1)
  Sigma=toeplitz(acf)
  y=matrix(MASS::mvrnorm(n, mu=numeric(p), Sigma=Sigma),n,p)
  expect_error(Toep.estimator(y=y,Te=Te,q=5,method="GCV",f.true=NULL))
  expect_error(Toep.estimator(y=y,Te=p+1,q=2,method="ML",f.true=NULL))
  expect_error(Toep.estimator(y=y,Te=Te,q=2,method="rML",f.true=NULL))
  expect_error(Toep.estimator(y=y,Te=Te,q=2,method="ML-oracle",f.true=NULL))
})

test_that("output type and length are correct - Toep.estimator",{
  n=100
  p=200
  Te=10
  q=3
  method="GCV"
  set.seed(11)
  acf=c(1.44,1.44/(1+seq(1,p-1))^2.1)
  Sigma=toeplitz(acf)
  y=matrix(MASS::mvrnorm(n, mu=numeric(p), Sigma=Sigma),n,p)
  expect_type(Toep.estimator(y=y,Te=Te,q=q,method=method,f.true=NULL), "list")
  expect_length(Toep.estimator(y=y,Te=Te,q=q,method=method,f.true=NULL), 4)
})

