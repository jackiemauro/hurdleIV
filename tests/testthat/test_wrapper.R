context("test wrapper fn")
test_that("wrapper runs", {
#   df = hurdle.IV.sim(type = "cragg")
#
#   expect_that(hurdle.IV(y~exog1+endog
#                         ,inst = inst1
#                         ,endog = endog
#                         ,exog = exog1
#                         ,data = df
#                         ,type = "cragg"),is_a("list"))

  df = hurdle.IV.sim(type = "lognormal")

  expect_that(hurdle.IV(y~exog1+endog
                        ,inst = inst1
                        ,endog = endog
                        ,exog = exog1
                        ,data = df
                        ,type = "lognormal"),is_a("list"))

})
