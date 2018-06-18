
#' cox.RPSFT
#'
#' Annalyzes a RPSFT object using a cox model. Corrects the standard errors based on rank test statistic used in G estimation.
#' @param x An RPSFT object created with RPSFT()
#' @param use.latent.only If TRUE experimental arm is imputed from latent time. If FALSE then observed times used for experimental arm.
#' @keywords RPSFT, Survival
#' @import dplyr
#' @import survival
#' @export
#' @examples
#'
#' sim.df <- simStudy()
#' x <- SurvExt(Surv(os.t, os.e) ~ I(x.trt==1),
#'    Exposure = ifelse(x.trt == 1, os.t, ifelse(x.switch == 1, os.t - t.switch, 0)),
#'    AdminCensTime = t.censor,
#'    data = sim.df)
#' y <- RPSFT(x)
#' cox.RPSFT(y)

cox.RPSFT <- function(x, use.latent.only = FALSE){

  # combine the observed and latent times
  df <- transmute(x$latent, latent.event.time = event.time, latent.censor.ind = censor.ind, psi = x$psi.chosen) %>%
    cbind(as.data.frame(x$input)) %>%
    mutate(event.time = t.on + t.off)

  # derive counterfactual times from latent event times and as T vs U
  df <- mutate(df,
              # need to convert CFExposure to latent event time scale for censoring purposes
              lat.on  = pmin(CFExposure * exp(psi), latent.event.time),
              lat.off = ifelse((latent.event.time - lat.on) > 0, latent.event.time - lat.on, 0),
              # convert back to real time
              cfact.time       = lat.on * exp(-psi) + lat.off,
              cfact.censor.ind = latent.censor.ind
              )

  # select

  # derive corrected hazard ratios
  cox.rpsft <- coxph(Surv(cfact.time, cfact.censor.ind) ~ trt.ind, data = df)

  cox.rpsft$call <- "RPSFT Analysis. Std Err adjusted per White method."

  # correct standard error using Z value
  cox.rpsft$var.orig <- cox.rpsft$var
  z0 <- survdiff(Surv(event.time, censor.ind) ~ trt.ind, data = df, rho = x$Grho)$chisq ^ 0.5
  cox.rpsft$var <- matrix((abs(coef(cox.rpsft)) / z0)^2)

  rc <- list(cox.rpsft = cox.rpsft, cfact.time = df$cfact.time, cfact.censor.ind = df$cfact.censor.ind, psi.unique = x$psi.unique)

  return(cox.rpsft)

}
