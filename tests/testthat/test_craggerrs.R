context("test cragg error functions")
test_that("makes list", {
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

test_that("right mean/sd", {
  cov = matrix(c(1,.2,.1,.2,4,.3,.1,.3,16),ncol = 3, byrow = T)
  n = 100
  x1 = as.data.frame(rnorm(n, 3, 1))
  z = as.data.frame(rnorm(n, -2,1))
  expect_equal(sd(cragg_errs(cov=cov,pi=c(1,2,3)
                             ,x1 = x1, gamma = c(.1,.2,.3)
                             ,beta = c(.2,.3,.5)
                             ,n=n,z = z)$endog)
               ,4
               ,tolerance = 2)
  expect_equal(min(cragg_errs2(cov=cov,pi=c(1,2,3)
                               ,x1 = x1, gamma = c(.1,.2,.3)
                               ,beta = c(.2,.3,.5)
                               ,n=n,z = z)$yStar),
               0,
               tolerance = 3)
})
