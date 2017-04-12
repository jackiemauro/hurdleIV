context("testing the reconstitute covariance matrix fn")
test_that("our function works", {
  expect_equal(reconstitute.cov(c(2,4,5,6,7),num = 1,chol = F)
               , matrix(c(1,2,5,2,4,6,5,6,7),ncol =3,byrow =F))
  
})

context("test start function")
test_that("start function works", {
  expect_equal(reconstitute.cov(c(2,4,5,6,7),num = 1,chol = F)
               , matrix(c(1,2,5,2,4,6,5,6,7),ncol =3,byrow =F))
  
})