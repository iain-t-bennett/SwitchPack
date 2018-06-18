
#' SurvCF - Counter Factual Survival Object
#'
#' Creates an SurvCF object that combines survival and counter factual exposure information. Derived from an RPSFT model.
#' @param x An RPSFT object created with RPSFT().
#' @param CFExposure Defines the counter factual exposure to study drug per patient.
#' @param UseObserved TRUE/FALSE indicator per patient if Observed values should be used for a patient.
#' @keywords RPSFT, Survival
#' @import dplyr
#' @import survival
#' @export
#' @examples
#'
#' sim.df <- simStudy()
#' x <- SurvExt(Surv(os.t, os.e) ~ I(x.trt==1),
#'    CFExposure = ifelse(x.trt == 1, os.t, ifelse(x.switch == 1, os.t - t.switch, 0)),
#'    AdminCensTime = t.censor,
#'    data = sim.df)
#' y <- RPSFT(x)
#' z <- SurvCF(y, CFExposure = ifelse(sim.df$x.trt == 0, 0, sim.df$os.t))
#' z <- SurvCF(y, CFExposure = 0, UseObserved = sim.df$x.trt == 1)

SurvCF <- function(x, CFExposure, UseObserved = FALSE){

  if (class(x)!= "RPSFT"){
    stop("x must be an RPSFT object")
  }

  # process the variables
  #indx2 <- match(c("CFExposure"), names(Call), nomatch = 0)

  #if (indx2[1] == 0) {
  #  stop("a CFExposure argument is required")
  #}

  # combine the latent and observed survival times
  df <- transmute(x$latent,
                  latent.event.time = event.time,
                  latent.censor.ind = censor.ind,
                  psi = x$psi.chosen) %>%
    cbind(as.data.frame(x$input)) %>%
    mutate(event.time = t.on + t.off)

  # derive counter factual times based on CFExposure
  df <- mutate(df,
               CFExpo           = CFExposure,
               type             = ifelse(UseObserved==TRUE, "Observed", "Counterfactual"),
               # derive latent time on and off
               lat.on           = pmin(t.on * exp(psi), latent.event.time),
               lat.off          = ifelse((latent.event.time - lat.on) > 0, latent.event.time - lat.on, 0),
               # derive counter factual time on and off (checked against the latent time scale)
               cfact.on         = pmin(CFExposure * exp(psi), latent.event.time) * exp(-psi),
               cfact.off        = pmax(latent.event.time * exp(-psi) - cfact.on,0),
               # derive the counterfactuals/observed
               cfact.time       = ifelse(type == "Counterfactual", cfact.on + cfact.off,  event.time),
               cfact.censor.ind = ifelse(type == "Counterfactual", latent.censor.ind, censor.ind)
               )

  #create object

  rc <- df

  attr(rc, "class") <- "SurvCF"

  return(rc)
}

#' as.data.frame.SurvCF
#'
#' Converts an SurvCF object to a data frame.
#' @param x An SurvCF object created with SurvCF()
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
#' as.data.frame(x)

as.data.frame.SurvCF <- function(x){
  data_frame(trt.ind = x$trt.ind,
             type = x$type,
             t.on = x$t.on,
             cfact.on = x$cfact.on,
             cfact.event.time = x$cfact.time,
             cfact.censor.ind = x$cfact.censor.ind,
             event.time = x$event.time,
             censor.ind = x$censor.ind,
             latent.event.time = x$event.time,
             latent.censor.ind = x$censor.ind
  )
}

#' plot.SurvCF
#'
#' Plots a SurvCF object using GGally and some default settings.
#' @param x An SurvCF object created with SurvCF()
#' @keywords RPSFT, Survival
#' @import dplyr
#' @import survival
#' @import GGally
#' @export
#' @examples
#'
#' sim.df <- simStudy()
#' x <- SurvExt(Surv(os.t, os.e) ~ I(x.trt==1),
#'    Exposure = ifelse(x.trt == 1, os.t, ifelse(x.switch == 1, os.t - t.switch, 0)),
#'    AdminCensTime = t.censor,
#'    data = sim.df)
#' plot(x)

plot.SurvCF <- function(x){
  SurvCF.KMPlot(x)
}


#' print.SurvCF
#'
#' Prints a summary of n rpsft.input object to screen.
#' @param x An rpsft.input object created with rpsft.input()
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
#' x

print.SurvCF <- function(x){

  rc <- mutate(as.data.frame(x), TreatmentInd = trt.ind) %>%
    group_by(TreatmentInd) %>%
    summarise(NPatients = n(),
              NEvents = sum(censor.ind),
              MeanExposure = mean(t.on),
              MeanCFExposure = mean(cfact.on)
              )

  cat("An analysis ready RPSFT object\n")
  print(rc)

}
