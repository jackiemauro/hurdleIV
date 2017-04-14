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

context("test cragg error functions")
test_that("makes vector", {
  cov = matrix(c(1,.2,.1,.2,4,.3,.1,.3,16),ncol = 3, byrow = T)
  n = 100
  x1 = as.data.frame(rnorm(n, 3, 1))
  z = as.data.frame(rnorm(n, -2,1))
  expect_that(cragg_errs(cov=cov,pi=c(1,2,3)
                         ,x1 = x1, gamma = c(.1,.2,.3)
                         ,beta = c(.2,.3,.5)
                         ,n=n,z = z)
              ,is_a("list"))
  expect_that(cragg_errs2(cov=cov,pi=c(1,2,3)
                         ,x1 = x1, gamma = c(.1,.2,.3)
                         ,beta = c(.2,.3,.5)
                         ,n=n,z = z)
              ,is_a("list"))
})
