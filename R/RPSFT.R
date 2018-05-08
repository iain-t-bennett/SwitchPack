#' A function to do RPSFT
#'
#' This function applies RPSFT method
#' @param rpsft.input An object created using rpsft.input()
#' @param psimat A matrix of values to try in grid search
#' @param Grho Defines the test used in the grid search (0 = logrank, 1 = wilcoxon)
#' @param twopass Defines if two searches are performed (wide then fine)
#' @keywords RPSFT, Survival
#' @import dplyr
#' @import survival
#' @import GGally
#' @export
#' @examples
#'
#' sim.df <- simStudy()
#' x <- rpsft.input(Surv(os.t, os.e) ~ I(x.trt==1),
#'    Exposure = ifelse(x.trt == 1, os.t, ifelse(x.switch == 1, os.t - t.switch, 0)),
#'    AdminCensTime = t.censor,
#'    data = sim.df)
#' plot(x)
#' RPSFT(x)


RPSFT <- function(rpsft.input,
                  psimat = matrix(ncol = 1, data = seq(from=-3, to=3, by=0.1)),
                  Grho   = 0,
                  twopass = TRUE){
  if (twopass){
    rc <- RPSFT.2pass(rpsft.input = rpsft.input, pass1.psimat = psimat, Grho = Grho)
  } else{
    rc <- RPSFT.1pass(rpsft.input = rpsft.input, psimat = psimat, Grho = Grho)
  }

  rc$input <- rpsft.input
  rc$Grho <- Grho

  attr(rc, "class") <- "RPSFT"

  return(rc)

}


#' print.RPSFT
#'
#' Prints a summary of an RPSFT object to screen.
#' @param x An RPSFT object created with RPSFT()
#' @keywords RPSFT, Survival
#' @import dplyr
#' @import survival
#' @export
#' @examples
#'
#' sim.df <- simStudy()
#' x <- rpsft.input(Surv(os.t, os.e) ~ I(x.trt==1),
#'    Exposure = ifelse(x.trt == 1, os.t, ifelse(x.switch == 1, os.t - t.switch, 0)),
#'    AdminCensTime = t.censor,
#'    data = sim.df)
#' y <- RPSFT(x)
#' y

print.RPSFT <- function(x){

  cat("RPSFT analysis performed using rank test with rho =",x$Grho,"\n")
  cat(length(x$psi.tried), " values for psi tested. ")
  if(x$psi.unique){
    cat("A unique solution of psi =", x$psi.chosen, " was found.")
  } else{
    cat("A unique solution for psi was not found. Possible solutions include: ", x$psi.found, ". Suggested value for further analysis is: ", x$psi.chosen)
  }

}

# function for the two pass

RPSFT.2pass <- function(rpsft.input, pass1.psimat = matrix(ncol = 1, data = c(seq(from=-3, to=3, by=0.1))), Grho = 0) {

  # apply RPSFT with wide grid
  pass1 <- RPSFT.1pass(rpsft.input, Grho=Grho, psimat = pass1.psimat)

  # get interesting places to look (within 95 percent CI for z)
  if (any(abs(pass1$z) <= 1.96)){
    pass1.limits <- pass1$psi.tried[abs(pass1$z) <= 1.96]
    pass2.psimat <- matrix(ncol = 1, data = c(seq(from = min(pass1.limits), to = max(pass1.limits), by = 0.01)))
  } else{
    pass2.psimat <- matrix(ncol = 1, data = c(seq(from = -5, to = 5, by = 0.05)))
  }

  # apply RPSFT again
  pass2 <- RPSFT.1pass(rpsft.input, Grho=Grho, psimat = pass2.psimat)

  # combine the grids searched for return
  pass.df <- rbind(data.frame(pass = 1, psi=pass1$psi.tried, z = pass1$z),
                   data.frame(pass = 2, psi=pass2$psi.tried, z = pass2$z))

  pass.df <- summarise(group_by(pass.df, psi , z), pass = min(pass))

  # return the second pass but add in the values tried in the first pass
  rc <- pass2
  rc$psi.tried <- pass.df$psi
  rc$z         <- pass.df$z
  rc$pass      <- pass.df$pass

  return(rc)
}

# does actual RPSFT grid search

RPSFT.1pass<-function(rpsft.input,
                psimat      = matrix(ncol = 1,
                                     data = seq(from=-2, to=2, by=0.01)),
                Grho        = 0

){

  # derive z values for a range of psi
  z <- apply(psimat,c(1),RPSFT.trypsi, Grho = Grho, rpsft.input = as.data.frame(rpsft.input))

  # chosen value for psi is that which has rank test Z value of 0
  # select based on when sign changes
  x <- sign(z)
  x1 <- x[-1]
  x2 <- x[-length(x)]

  # selected psi index (either is 0 or is between a sign change)
  indx         <- which(z==0)
  indx_chg_low <- which(x1!= x2 & x1!=0 & x2!=0)
  indx_chg_hgh <- indx_chg_low+1
  psi.found    <- c(psimat[indx,1], 0.5*(psimat[indx_chg_low,1]+psimat[indx_chg_hgh,1]))

  # is it a unique solution?
  psi.unique <- as.numeric(length(psi.found)==1)

  # apply proposal of White to take weighted average if is multiple
  psi.chosen <- sum(c(1,rep(c(-1,1),0.5*(length(psi.found)-1)))*psi.found)

  # get latent survival time for chosen value
  latent <- RPSFT.latent(psi.chosen, rpsft.input)

  # return the tested values, associated Z values, chosen psi and latent survival
  rc <- list(psi.tried  = as.numeric(psimat),
             z          = z,
             psi.found  = psi.found,
             psi.unique = psi.unique,
             psi.chosen = psi.chosen,
             latent     = latent,
             pass       = 1
  )
  return(rc)
}

# returns latent survival times for given psi and standard model
# inputs are:
#   psi - values of psi to try
#   rpsft.input - dataframe/list containing:
#     t.on        = time on experimental treatment
#     t.off       = time off experimental treatment
#     censor.ind  = censor indicator (1 = event)
#     trt.ind     = treatmemt indoicator (1 = randomized to experimental)
#     cutofftime  = time of admin censor


RPSFT.latent <- function(psi, rpsft.input){
  # recensor time
  latent<-data.frame(cstar = pmin(rpsft.input$cutofftime*exp(psi), rpsft.input$cutofftime))
  # counterfactual time
  latent$event.time <- pmin(rpsft.input$t.off + rpsft.input$t.on*exp(psi), latent$cstar )
  # apply recensoring
  latent$censor.ind <- as.numeric((latent$event.time < latent$cstar) & rpsft.input$censor.ind)
  return(latent)
}

# inputs are:
#   psi   - value of psi to try
#   Grho  -  defines the test used in the grid search (0 = logrank, 1 = wilcoxon)
#   rpsft.input - dataframe/list containing:
#     t.on        = time on experimental treatment
#     t.off       = time off experimental treatment
#     censor.ind  = censor indicator (1 = event)
#     trt.ind     = treatmemt indoicator (1 = randomized to experimental)
#     cutofftime  = time of admin censor

# tests for equality of latent survival times
RPSFT.trypsi <- function(psi, Grho, rpsft.input){

  # define latent survival time given psi
  latent <- RPSFT.latent(psi, rpsft.input)

  # derive Z value for chosen test
  sv <- Surv(latent$event.time, latent$censor.ind)
  sd <- survdiff(sv ~ rpsft.input$trt.ind, rho = Grho)
  rc <- sign(sd$obs[2]-sd$exp[2])*sd$chisq^0.5
  return(rc)
}
