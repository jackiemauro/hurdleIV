context("test lognormal likelihood fn")
test_that("makes value", {
  df = hurdle.IV.sim()
  outcome = df$y
  ER_mat = list(as.matrix(data.frame(endog = df$endog, exog1 = df$exog1, inst1 = df$inst1)))
  y_mat = as.matrix(data.frame(outcome = outcome, exog1 = df$exog1, endog = df$endog))
  endog_mat = as.matrix(df$endog)
  num_endog = 1; myChol = T; len_cov = 5; num_betas = 3; numpis = 3
  g = c(.2,3,.1,.4,9,1,1,1,2,2,2,3,3,3)

  att = list(df=df,outcome=outcome,ER_mat=ER_mat
             ,y_mat=y_mat,endog_mat=endog_mat,num_endog=num_endog
             ,myChol=myChol,len_cov=len_cov,num_betas=num_betas,numpis=numpis)
  attach(att)

  expect_that(loglik_lgiv(g),is_a("numeric"))

  expect_error(loglik_lgiv(g[-1]))

  detach(att)

})

