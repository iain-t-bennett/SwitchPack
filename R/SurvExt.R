
#' SurvExt - Extended Survival Object
#'
#' Creates an SurvExt object that combines survival and exposure information. Can be used for RPSFT analysis.
#' @param formula A formula for survival time on randomized treatment. Can only include single regression variable.
#' @param data Required. A data frame in which to interpret the variables named in the formula and other arguments.
#' @param Exposure Defines the exposure to study drug per patient.
#' @param AdminCensTime Defines the administrative censoring time.
#' @param FirstSwitchTime Defines the first time of switch
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
#' plot(x)
#' RPSFT(x)
SurvExt <- function(formula, data , Exposure, AdminCensTime, FirstSwitchTime){

  Call <- match.call()

  indx <- match(c("formula", "data"), names(Call), nomatch = 0)
  if (indx[1] == 0) {
    stop("a formula argument is required")
  }

  temp <- Call[c(1,indx)]
  temp[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(temp)

  Terms <- terms(formula)
  ll <- attr(Terms, "term.labels")
  if (length(ll) != 1)
    stop("Only randomization indicator can be included in formula")

  # get response and Surv object
  Y <- model.extract(m, "response")

  if (!is.Surv(Y))
    stop("Response must be a survival object")

  sv <- with(data, Y)

  # process the additional variables
  indx2 <- match(c("Exposure", "AdminCensTime", "data", "FirstSwitchTime"), names(Call), nomatch = 0)

  if (indx2[1] == 0) {
    stop("an Exposure argument is required")
  }
  if (indx2[2] == 0) {
    stop("an AdminCensTime argument is required")
  }
  temp2 <- Call[indx2]

  t.on = transmute_(data, temp2[[1]])[,1]
  t.cens = transmute_(data, temp2[[2]])[,1]
  t.off = sv[,1]-t.on
  trt.ind <- as.numeric(m[,2])

  # process the optional variables
  indx3 <- match(c("FirstSwitchTime"), names(Call), nomatch = 0)

  # First switch time is provided
  if (indx3[1] != 0) {
    temp3 <- Call[indx3]
    t.start = transmute_(data, temp3[[1]])[,1] %>%
      as.numeric
  } else{
    t.start <- NA
  }

  #create object

  rc <- list(trt.ind = as.numeric(trt.ind),
             t.on = as.numeric(t.on),
             t.off = as.numeric(t.off),
             censor.ind = as.numeric(sv[,2]),
             cutofftime = as.numeric(t.cens),
             t.start = t.start,
             orig.sv = sv)

  attr(rc, "class") <- "SurvExt"

  return(rc)
}

#' as.data.frame.SurvExt
#'
#' Converts an SurvExt object to a data frame.
#' @param x An SurvExt object created with SurvExt()
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

as.data.frame.SurvExt <- function(x){
  data_frame(trt.ind = x$trt.ind,
             t.on = x$t.on,
             t.off = x$t.off,
             censor.ind = x$censor.ind,
             cutofftime = x$cutofftime,
             t.start = x$t.start
  )
}

#' plot.SurvExt
#'
#' Plots a SurvExt object using GGally and some default settings.
#' @param x An SurvExt object created with SurvExt()
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

plot.SurvExt <- function(x){
  SurvExt.KMPlot(x)
}


#' print.SurvExt
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

print.SurvExt <- function(x){

  rc <- mutate(as.data.frame(x), TreatmentInd = trt.ind) %>%
    group_by(TreatmentInd) %>%
    summarise(NPatients = n(),
              NEvents = sum(censor.ind),
              MeanExposure = mean(t.on))
  cat("An analysis ready RPSFT object\n")
  print(rc)

}
