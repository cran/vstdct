test_that("output type is correct", {
  expect_true(
    is.matrix(DCT.matrix(10))
  )
})

test_that("matrix is orthogonal", {
  expect_true(
    sum(DCT.matrix(10)%*%DCT.matrix(10) - diag(10))< .Machine$double.eps
  )
})

test_that("acf2sdf is inverse operation of sdf2acf", {
  x=rnorm(20)
  expect_equal(
    sdf2acf(acf2sdf(x))
    ,x)
  expect_equal(
    acf2sdf(sdf2acf(x))
    ,x)
})
