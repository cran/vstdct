test_that("valid arguments are passed - example 1", {
  expect_error(
    example1(
      p=10,
      n=20,
      sd=1
    )
  )
})

test_that("valid arguments are passed - example 2", {
  expect_error(
    example2(
      p=10,
      n=20,
      family="Gamma"
    )
  )
})

test_that("valid arguments are passed - example 2", {
  expect_error(
    example3(
      p=10,
      n=20,
      sd=1.2,
      gamma=1.3,
      family="uniform"
    )
  )
})

test_that("output type and length are correct - example 1", {
  expect_type(
    example1(
      p=10,
      n=20,
      sd=1,
      gamma=2,
      family="Gaussian"
    ), "list"
  )
  expect_length(
    example1(
      p=10,
      n=20,
      sd=1,
      gamma=2,
      family="Gaussian"
    ), 3
  )

})

test_that("output type and length are correct - example 2", {
  expect_type(
    example2(
      p=10,
      n=20,
      sd=1,
      family="Gamma"
    ), "list"
  )
  expect_length(
    example2(
      p=10,
      n=20,
      sd=1,
      family="Gamma"
    ),3
   )
})

test_that("output type and length are correct - example 3", {
  expect_type(
  example3(
    p=10,
    n=20,
    sd=1,
    gamma=3.1,
    family="Uniform"
  ), "list"
)
  expect_length(
    example3(
      p=10,
      n=20,
      sd=1,
      gamma=3.1,
      family="Uniform"
    ), 3
  )
})
