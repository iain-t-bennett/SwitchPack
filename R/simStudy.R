#' A function to simulate Survival data
#'
#' This function simulates survival data with correlated time to progression and overall survival times.
#' @param simid Identifer added to simulated dataframe. Defaults to 1.
#' @param seed Seed used for random number generator. Defaults to 1234.
#' @param rho correlation coefficient between TTP and OS. Defaults to 0.
#' @param p.switch proportion of patients who switch. Defaults to 0.25.
#' @param prop.pd.int proportion of patients with PFS before switch allowed. Defaults to 0.
#' @param beta_1a treatment effect (as log(Hazard Ratio)). Defaults to log(0.7).
#' @param beta_1b treatment effect (as log(Hazard Ratio)). Defaults to log(0.7).
#' @param beta_2a treatment effect (as log(Hazard Ratio)). Defaults to log(0.7).
#' @param beta_2b treatment effect (as log(Hazard Ratio)). Defaults to log(0.7).
#' @param beta_pd treatment effect on prog (as log(HR)). Defaults to log(0.4).
#' @param arm.n patients per arm. Defaults to 250.
#' @param enroll.start start of enrollment. Defaults to 0.
#' @param enroll.end end of enrollment. Defaults to 2.
#' @param admin.censor end of trial. Defaults to 3.
#' @param os.gamma weibull shape - for OS. Defaults to 1.2.
#' @param os.lambda weibull scale - for OS. Defaults to 0.3.
#' @param ttp.gamma weibull shape - for TTP. Defaults to 1.5.
#' @param ttp.lambda weibull scale - for TTP. Defaults to 2.
#' @param include.metadata Should returned dataframe include simulation parameters. Defaults to False.
#' @keywords RPSFT, Survival, Simulation
#' @import Hmisc
#' @import dplyr
#' @export
#' @examples
#' require(survival)
#' require(dplyr)
#'
#' tst.df <- simStudy()
#' survfit(Surv(pfs.t, event = pfs.e) ~ x.trt, data = tst.df) %>%
#'   plot()
#'
#' survfit(Surv(os.t, event = os.e) ~ x.trt, data = tst.df) %>%
#'   plot()

simStudy <- function(
  # params modified by sim
  simid         = 1,         # identifier
  seed          = 1234,      # seed
  rho           = 0,         # correlation coefficient between switch time and os
  p.switch      = 0.25,      # proportion switch
  prop.pd.int   = 0,         # proportion of patients with PFS before switch allowed
  beta_1a       = log(0.7),  # treatment effect (as log(Hazard Ratio))
  beta_1b       = log(0.7),  # treatment effect (as log(Hazard Ratio))
  beta_2a       = log(0.7),  # treatment effect (as log(Hazard Ratio))
  beta_2b       = log(0.7),  # treatment effect (as log(Hazard Ratio))
  # params not changed across all sims
  beta_pd       = log(0.4),  # treatment effect on prog (as log(HR))
  arm.n         = 250,       # patients per arm
  enroll.start  = 0,         # start of enrollment
  enroll.end    = 2,         # end of enrollment
  admin.censor  = 3,         # end of trial
  os.gamma      = 1.2,       # weibull shape - for OS
  os.lambda     = 0.3,       # weibull scale - for OS
  ttp.gamma     = 1.5,       # weibull shape - for TTP
  ttp.lambda    = 2,         # weibull scale - for TTP,
  include.metadata = FALSE
){
  set.seed(seed)
  # combine so can select based on x
  n_switch <- arm.n * p.switch
  # u1 is just used to setup then not used
  # u2 and u3 are correlated uniform random variables used for OS and switch time
  # u4 used for censoring time (entry time)
  # u5 used for selecting switchers
  u <- matrix(nrow = arm.n*2, ncol = 5,data = runif(arm.n*2*5,0,1))
  # setup correlations using normal variables then convert back to uniform
  u[,3] <- pnorm(rho* qnorm(u[,2],0,1) + ((1-rho^2)^0.5)*qnorm(u[,1],0,1),0,1)
  # generate covariates
  x.trt    <- rep(c(0,1), each = arm.n)
  # generate os without treatment
  os.basis.t  <- (-(log(u[,2]) / (os.lambda)))^(1/os.gamma)
  # generate ttp
  ttp.untreated.t <- (-(log(u[,3]) / (ttp.lambda)))^(1/ttp.gamma)
  # apply treatment effect
  ttp.treated.t <-(-(log(u[,3]) / (ttp.lambda*exp(beta_pd))))^(1/ttp.gamma)
  # deive basis ttp
  ttp.basis.t <- ifelse(x.trt==1, ttp.treated.t, ttp.untreated.t)
  # generate pfs
  pfs.basis.t <- pmin(os.basis.t, ttp.basis.t)
  # generate entry time
  t.entry <- enroll.start + runif(u[,4])*(enroll.end-enroll.start)
  # generate follow up time
  t.censor <- admin.censor - t.entry
  # define a time in trial after which switch happens
  t.start.switch <- as.numeric(quantile(t.entry + pfs.basis.t, prop.pd.int))
  # who can switch - have ttp before os and censor
  can.switch <- as.numeric(
    ttp.basis.t < os.basis.t &
      ttp.basis.t < t.censor &
        t.entry+ttp.basis.t >= t.start.switch
  )[x.trt==0]

  rank <- rep(arm.n,arm.n)
  # choose randomly from those who can switch to get the reqd number of switchers
  rank[can.switch==1] <- order(u[x.trt==0 & can.switch==1,5])
  x.switch <- c(as.numeric(rank<=n_switch), rep(0,arm.n))
  # calculate how many actually switch
  actual.p.switch <- sum(x.switch) / arm.n
  # generate switch time
  t.switch <- rep(NA,arm.n*2)
  t.switch[x.switch==1] <- ttp.basis.t[x.switch==1]

  # get control non switch os
  os.orig.t <- rep(NA,arm.n*2)
  os.orig.t[x.trt == 0] <- (-(log(u[x.trt==0,2])) / (os.lambda ))^(1/os.gamma)

  # get experimental arm os
  t0  <- ttp.basis.t
  t0v <- t0 ^ os.gamma
  c1 <- os.lambda * exp(beta_1a)
  c2 <- os.lambda * exp(beta_1b)
  v1 <- 1/os.gamma
  os.orig.t[x.trt == 1] <- ifelse(-log(u[x.trt==1,2])  < c1 * t0v[x.trt==1],
                                 (-log(u[x.trt==1,2]) / c1) ^ v1,
                                 ((-log(u[x.trt==1,2]) - c1*t0v[x.trt==1] + c2*t0v[x.trt==1]) / c2) ^ v1
                                  )

  # get switch treatment (start from orig)
  os.cross.t <- os.orig.t
  t1 <- ttp.basis.t
  t2 <- ttp.basis.t + ttp.treated.t
  t1v <- t1^os.gamma
  t2v <- t2^os.gamma
  c3 <- os.lambda*exp(beta_2a)
  c4 <- os.lambda*exp(beta_2b)
  d1 <- x.switch==1 & -log(u[,2]) < os.lambda*t1v
  d3 <- x.switch==1 & -log(u[,2]) >= os.lambda*t1v + c3*(t2-t1)
  d2 <- x.switch==1 & -log(u[,2]) >= os.lambda*t1v & !d3
  os.cross.t[d1] <- (-log(u[d1,2]) / os.lambda )^v1
  os.cross.t[d2] <- ((-log(u[d2,2])-os.lambda*t1v[d2] + c3*t1v[d2]) / c3 )^v1
  os.cross.t[d3] <- ((-log(u[d3,2])-os.lambda*t1v[d3] - c3*(t2v[d3] - t1v[d3]) +c4*t2v[d3]) / c4 )^v1

  # apply censoring
  os.b.t <- pmin(os.basis.t, t.censor)
  os.b.c <- as.numeric(os.b.t==t.censor)
  os.o.t <- pmin(os.orig.t, t.censor)
  os.o.c <- as.numeric(os.o.t==t.censor)
  os.t <- pmin(os.cross.t, t.censor)
  os.c <- as.numeric(os.t==t.censor)

  # store ttp
  ttp.t <- pmin(ttp.basis.t, t.censor)
  ttp.c <- as.numeric(ttp.t==t.censor)

  # derive PFS (min of os and ttp)
  pfs.t <- pmin(os.t, ttp.t)
  pfs.c <- as.numeric(pfs.t==t.censor)

  # create a data frame
  patid <- 1:(arm.n*2)
  df <- data.frame(patid    = patid,
                   x.trt    = x.trt,
                   x.switch = x.switch,
                   t.switch = t.switch,
                   t.censor = t.censor,
                   pfs.t    = pfs.t,
                   pfs.e    = 1 - pfs.c,
                   os.t     = os.t,
                   os.e     = 1 - os.c,
                   ttp.t    = ttp.t,
                   ttp.e    = 1 - ttp.c,
                   os.b.t   = os.b.t,
                   os.b.e   = 1 - os.b.c,
                   os.o.t   = os.o.t,
                   os.o.e   = 1 - os.o.c,
                   seed     = seed,
                   rho      = rho,
                   prop.pd.int = prop.pd.int,
                   beta_1a = beta_1a,
                   beta_1b = beta_1b,
                   beta_2a = beta_2a,
                   beta_2b = beta_2b,
                   beta_pd = beta_pd,
                   ttp.lambda = ttp.lambda,
                   ttp.gamma = ttp.gamma,
                   os.lambda = os.lambda,
                   os.gamma = os.gamma,
                   simid   = simid,
                   p.switch = p.switch
  )

  # apply labels

  var.labels <- c(
    patid    = "Patient ID",
    x.trt    = "Randomized Treatment (0 = control, 1 = experimental)",
    x.switch = "Switch Indicator (0 = switch to experimental)",
    t.switch = "Time of switch",
    t.censor = "Time of administrative censoring",
    pfs.t    = "Progression Free Survival (PFS) Time",
    pfs.e    = "PFS Censoring indicator (1 = event)",
    os.t     = "Overall Survival (OS) Time",
    os.e     = "OS Censoring indicator (1 = event)",
    ttp.t    = "Time To Progression (TTP)",
    ttp.e    = "TTP censoring indicator (1 = event)",
    os.b.t   = "Unobserved Basis Overall Survivial (OS B)",
    os.b.e   = "Unobserved OS B censoring indicator (1 = event)",
    os.o.t   = "Unobserved Overall Survivial without switch (OS TRUE)",
    os.o.e   = "Unobserved OS TRUE censoring indicator (1 = event)",
    seed     = "Seed for simulation",
    rho      = "Simulation Parameter: rho - Correlation between TTP and OS",
    prop.pd.int = "Simulation Parameter: prop.pd.int",
    beta_1a = "Simulation Parameter: beta_1a",
    beta_1b = "Simulation Parameter: beta_1b",
    beta_2a = "Simulation Parameter: beta_2a",
    beta_2b = "Simulation Parameter: beta_2b",
    beta_pd = "Simulation Parameter: beta_pd",
    os.lambda = "Simulation Parameter: weibull scale for OS",
    os.gamma = "Simulation Parameter: weibull shape for OS",
    ttp.lambda  = "Simulation Parameter: weibull scale for TTP",
    ttp.gamma = "Simulation Parameter: weibull shape for TTP",
    simid   = "Simulation Parameter: simid",
    p.switch ="Simulation Parameter: p.switch"
  )

  df <- upData(df, labels = var.labels)

  # select variables to return

  if (include.metadata) {
    rc <- df
  } else {
    rc <- transmute(df,  patid, x.trt, x.switch, t.switch, t.censor, pfs.t, pfs.e, os.t, os.e)
  }

  return(rc)
}



