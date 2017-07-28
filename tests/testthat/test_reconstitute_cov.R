context("testing the reconstitute covariance matrix fn")
#len_cov = 2; num_betas = 1; num_pis=1; num_endog = 1
test_that("our function works", {
  expect_equal(reconstitute.cov(c(2,4,5,6,7),num = 1,chol = F)
               , matrix(c(1,2,5,2,4,6,5,6,7),ncol =3,byrow =F))
  
})