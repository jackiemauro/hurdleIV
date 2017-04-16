context("test simulation fn")
test_that("makes dataframe", {
  expect_that(hurdle.IV.sim(),is_a("data.frame"))
  expect_that(hurdle.IV.sim(type="cragg"),is_a("data.frame"))
  expect_equal(dim(hurdle.IV.sim(n=10))[1],10)
})
