
#' rpsft.input
#'
#' Creates an rpsft.input object that combines surival and exposure information.
#' @param formula A formula for survival time on randomized treatment. Can only include single regression variable.
#' @param data Required. A data frame in which to interpret the variables named in the formula and other arguments.
#' @param Exposure Defines the exposure to study drug per patient.
#' @param AdminCensTime Defined the administrative censoring time.
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
#' plot(x)
#' RPSFT(x)
rpsft.input <- function(formula, data , Exposure, AdminCensTime){
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
  indx2 <- match(c("Exposure", "AdminCensTime", "data"), names(Call), nomatch = 0)

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

  rc <- list(trt.ind = as.numeric(trt.ind),
             t.on = as.numeric(t.on),
             t.off = as.numeric(t.off),
             censor.ind = as.numeric(sv[,2]),
             cutofftime = as.numeric(t.cens),
             orig.sv = sv)

  attr(rc, "class") <- "rpsft.input"

  return(rc)
}

#' as.data.frame.rpsft.input
#'
#' Converts an rpsft.input object to a data frame.
#' @param x An rpsft.input object created with rpsft.input()
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
#' as.data.frame(x)

as.data.frame.rpsft.input <- function(x){
  data_frame(trt.ind = x$trt.ind,
             t.on = x$t.on,
             t.off = x$t.off,
             censor.ind = x$censor.ind,
             cutofftime = x$cutofftime
  )
}

#' plot.rpsft.input
#'
#' Plots an rpsft.input object using GGally and some default settings.
#' @param x An rpsft.input object created with rpsft.input()
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

plot.rpsft.input <- function(x){
  Treatment <- factor(x$trt.ind, levels = c(1,0), labels = c("Experimental", "Control Unadjusted"), ordered = TRUE)

  ggsurv(survfit(x$orig.sv ~ Treatment)) +
    theme_bw() +
    coord_cartesian(ylim = c(0,1))
}


#' print.rpsft.input
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
#' x <- rpsft.input(Surv(os.t, os.e) ~ I(x.trt==1),
#'    Exposure = ifelse(x.trt == 1, os.t, ifelse(x.switch == 1, os.t - t.switch, 0)),
#'    AdminCensTime = t.censor,
#'    data = sim.df)
#' x

print.rpsft.input <- function(x){

  rc <- mutate(as.data.frame(x), TreatmentInd = trt.ind) %>%
    group_by(TreatmentInd) %>%
    summarise(NPatients = n(),
              NEvents = sum(censor.ind),
              MeanExposure = mean(t.on))
  cat("An analysis ready RPSFT object\n")
  print(rc)

}
