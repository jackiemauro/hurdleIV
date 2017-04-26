context("test cragg likelihood fn")
test_that("makes value", {
  df = hurdle.IV.sim(type = "cragg")
  outcome = df$y
  ER_mat = list(as.matrix(data.frame(endog = df$endog, exog1 = df$exog1, inst1 = df$inst1)))
  y_mat = as.matrix(data.frame(outcome = outcome, exog1 = df$exog1, endog = df$endog))
  endog_mat = as.matrix(df$endog)
  num_endog = 1; myChol = T; len_cov = 5; num_betas = 3; numpis = 3

  att = list(df=df,outcome=outcome,ER_mat=ER_mat
             ,y_mat=y_mat,endog_mat=endog_mat,num_endog=num_endog
             ,myChol=myChol,len_cov=len_cov,num_betas=num_betas,numpis=numpis
             ,endog_names = "endog1")
  attach(att)

  g = c(.2,3,.1,.4,9,1,1,1,2,2,2,3,3,3)

  expect_that(loglik_craggiv(g),is_a("numeric"))

  expect_error(loglik_craggiv(g[-1]))

  # test that the pmvnorm function does roughly what you want
  muMat = matrix(rep(0,8),ncol=2)
  l = c(-Inf,-Inf); u = c(Inf,1.96)
  f <- function(x){pmvnorm(lower = l, upper = u, mean = x, sigma = diag(2))}
  expect_equal(apply(muMat,1,f),rep(pnorm(1.96),4))

  detach(att)
})
