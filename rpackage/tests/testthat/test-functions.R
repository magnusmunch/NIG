context("miscellaneous functions")

test_that("ratio of modified Bessel functions", {

  # setting parameters
  x <- rchisq(1, 1)
  nu <- rchisq(1, 1)

  # implemented function
  test.ratio <- ratio_besselK_cpp(x, nu)

  # actual calculation (see supplement)
  ratio <- besselK(x, nu - 1)/besselK(x, nu)

  # check
  expect_equal(test.ratio, ratio)

})
