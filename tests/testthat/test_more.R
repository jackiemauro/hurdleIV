context("testing the reconstitute covariance matrix fn")
test_that("our function works", {
  expect_equal(reconstitute.cov(c(2,4,5,6,7),num = 1,chol = F)
               , matrix(c(1,2,5,2,4,6,5,6,7),ncol =3,byrow =F))

})

context("test make dataframe function gives dataframe")
test_that("data frame function gives dataframes", {
  expect_that(make.df(c(1,2,3),c(1,1,1),100),is_a("data.frame"))
})

context("test make dataframe function fails on wrong dim")
test_that("data frame function fails if mean and sd diff lengths", {
  expect_that(make.df(c(1,2,3),c(1,1),100),throws_error())
})

