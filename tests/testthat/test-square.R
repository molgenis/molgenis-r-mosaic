context("A sample test")

test_that("Some test", {
  res = square(3)

  expect_equal(res, 9)
})

