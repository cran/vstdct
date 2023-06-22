
test_that("valid arguments are passed - EF.poly", {
  n=30
  x=seq(1,n)/n
  expect_error(EF.poly(-1,x))
})

test_that("valid arguments are passed - Q_pminus1 ", {
  n=30
  x=seq(1,n)/n
  expect_error(Q_pminus1 (-1,x))
})

test_that("valid arguments are passed - DR.basis", {
  n=30
  x=seq(1,n)/n
  expect_error(DR.basis(x,n,q=1.1))
})

test_that("EF-polynomial of order 0 and 1 is constant 1", {
  n=100
  x=seq(1,n)/n
  expect_equal(EF.poly(0,x),rep(1,n))
  expect_equal(EF.poly(1,x),rep(1,n))
})

test_that("Qpminus-polynomial of order 0 is constant 1", {
  n=100
  x=seq(1,n)/n
  expect_equal(Q_pminus1(0,x),rep(1,n))
})

test_that("valid arguments are passed - DR.basis", {
  n=100
  x=seq(1,n)/n
  expect_error(DR.basis(x,n,q=1.1))
})

test_that("DR basis matrix is orthogonal", {
  n=100
  x=seq(1,n)/n
  expect_true(
    abs(sum(DR.basis(x,n,q=1)%*%t(Conj(DR.basis(x,n,q=1))) - n*diag(n)))< .Machine$double.eps*n^2
  )
})

