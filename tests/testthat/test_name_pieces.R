context("testing the name pieces function")
#len_cov = 2; num_betas = 1; num_pis=1; num_endog = 1
test_that("our function works", {
  expect_equal(name.pieces(c(1,2,3,4,5),
                           a1 = 2, a2 = 1, a3 = 1, a4=1)
               , list(cov_start=c(1,2),beta = 3,gamma = 4,pi = list(5)))
  expect_equal(name.pieces(c(1,2,3,4,5,6,7,8),
                           a1 = 4, a2 = 1, a3 = 1, a4=2)
               , list(cov_start=c(1,2,3,4)
                      ,beta = 5
                      ,gamma = 6
                      ,pi = list(c(7,8))))
  expect_equal(name.pieces(c(1,2,3,4,5,6,7,8),
                           a1 = 4, a2 = 1, a3 = 2, a4=list(1,1))
               , list(cov_start=c(1,2,3,4)
                      ,beta = 5
                      ,gamma = 6
                      ,pi = list(7,8)))

})

