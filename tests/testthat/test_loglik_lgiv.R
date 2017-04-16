context("test lognormal likelihood fn")
test_that("makes value", {
  df = hurdle.IV.sim()
  outcome = df$y
  ER_mat = list(as.matrix(data.frame(endog = df$endog, exog1 = df$exog1, inst1 = df$inst1)))
  y_mat = as.matrix(data.frame(outcome = outcome, exog1 = df$exog1, endog = df$endog))
  endog_mat = as.matrix(df$endog)
  num_endog = 1; myChol = T; len_cov = 5; num_betas = 3; numpis = 3
  g = c(.2,3,.1,.4,9,1,1,1,2,2,2,3,3,3)

  expect_that(loglik_lgiv(g),is_a("numeric"))

  expect_error(loglik_lgiv(g[-1]))

})
g[ a]
