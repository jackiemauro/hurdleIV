### figure with iv tobit coefficients ##
# correct generalization needs to be figured out for new dgp
require(censReg)


parameters = list(pi = c(1,-1,3), #intercept, exog, inst
                  gamma = c(-.2,.8,.07), #intercept, exog, endog
                  beta = c(.05,.06,.02), #intercept, exog, endog
                  endog_reg = list(),
                  exog_mean = 1,
                  exog_sd = 1,
                  z_mean = 3,
                  z_sd = 1,
                  endog_sd = 1,
                  y_sd = 2,
                  rho = 0.2,
                  tau0 = .3,
                  tau1 = .1)

#tobit assumption values
Sig = make.covTrans(list(rho = parameters$rho
                         , y_sd = parameters$y_sd
                         , endog_sd = parameters$endog_sd
                         , tau0 = parameters$tau0
                         , tau1 = parameters$tau1)
                    , num_endog = 1
                    , gamma = parameters$gamma
                    , beta = parameters$beta
                    , option = "parameters", noname = T)
newpar = parameters
newpar$tau1 = 0
newpar$rho = 0
newpar$y_sd <- 1
sig_y0_x2 <- Sig[1,1] - Sig[1,3]%*%solve(Sig[3,3])%*%Sig[3,1]
assmp_b2 <- (parameters$gamma[3] - parameters$tau0/parameters$endog_sd^2)/sig_y0_x2
assmp_b1 <- parameters$gamma[1:2]/sig_y0_x2
assmp_g1 <- parameters$beta[1:2]*sig_y0_x2
assmp_g2 <- parameters$beta[3]*sig_y0_x2 + parameters$tau0/parameters$endog_sd^2
newpar$gamma = c(assmp_g1, assmp_g2)

gRange2 = seq(from = assmp_g2*.5, to = assmp_g2*2, length = 10)

tobit_beta2 = rep(NA,length(gRange2))
sd = rep(NA,length(gRange2))
n = 1000; m = 10
for(ii in 1:length(gRange2)){
  
  print(paste("round ",ii))
  gamma = parameters$gamma
  gamma[3] = gRange2[ii]
  vals = c(rep(NA,m))
  for(jj in 1:m){
    # simulate data
    dat = hurdle.IV.sim(formula = F,
                        n = n,
                        pi = newpar$pi,
                        gamma = gamma,
                        beta = newpar$beta,
                        exog_mean = newpar$exog_mean,
                        z_mean = newpar$z_mean,
                        exog_sd = newpar$exog_sd,
                        z_sd = newpar$z_sd,
                        endog_sd = newpar$endog_sd,
                        y_sd = newpar$y_sd,
                        rho = newpar$rho,
                        tau0 = newpar$tau0,
                        tau1 = newpar$tau1,
                        #options = list(silent = T, cragg_errors = 4),
                        type = "cragg")
    
    # iv tobit
    require(censReg)
    reduced.form <- lm(endog~exog1 + inst1, dat)
    consistent.tobit <- censReg(y~fitted(reduced.form)+residuals(reduced.form) + exog1,data = dat)
    vals[jj] = coef(consistent.tobit)['fitted(reduced.form)']
  }
  tobit_beta2[ii] = mean(vals, na.rm = T)
  sd[ii] = sd(vals, na.rm = T)
}

dat_tob = data.frame(Gamma = gRange2, Beta = tobit_beta2, err = sd)
g <- ggplot(dat_tob,aes(x = Gamma, y = Beta)) +
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=Beta-2*err, ymax=Beta+2*err), width=.1) +
  geom_hline(yintercept=parameters$beta[3], lty = 2)+
  geom_vline(xintercept=assmp_g2) +
  theme_bw()
g