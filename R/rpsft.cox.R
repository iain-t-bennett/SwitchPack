


RPSFT.cox <- function(rpsft.input, rpsft.output, Grho = 0, use.latent.only = TRUE){



  # combine the observed and latent times
  df <- cbind(rpsft.input, transmute(rpsft.output$latent, latent.event.time = event.time, latent.censor.ind = censor.ind, psi = rpsft.output$psi.chosen))

  # derive counterfactual times from latent event times and as T vs U
  if (use.latent.only){
    df <- mutate(df,
                 lat.on  = pmin(t.on * exp(psi), latent.event.time),
                 lat.off = ifelse((latent.event.time - lat.on) > 0, latent.event.time - lat.on, 0),
                 cfact.time       = ifelse(trt.ind == 1, lat.on * exp(-psi) + lat.off,  latent.event.time ),
                 cfact.censor.ind = latent.censor.ind
    )
  } else {
    df <- mutate(df,
                 cfact.time       = ifelse(trt.ind == 1, event.time,  latent.event.time ),
                 cfact.censor.ind = ifelse(trt.ind == 1, censor.ind,  latent.censor.ind )
    )
  }

  # derive corrected hazard ratios
  cox.rpsft <- coxph(Surv(cfact.time, cfact.censor.ind) ~ trt.ind, data = df)

  # correct standard error using Z value
  cox.rpsft$var.orig <- cox.rpsft$var
  z0 <- survdiff(Surv(event.time, censor.ind) ~ trt.ind, data = df, rho = Grho)$chisq ^ 0.5
  cox.rpsft$var <- matrix((abs(coef(cox.rpsft)) / z0)^2)

  rc <- list(cox.rpsft = cox.rpsft, cfact.time = df$cfact.time, cfact.censor.ind = df$cfact.censor.ind, psi.unique = rpsft.output$psi.unique)

  return(rc)

}
