#' SurvExt.KMPlot
#'
#' KM plot of a SurvExt object using GGally and some default settings.
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
#' SurvExt.KMPlot(x)

SurvExt.KMPlot <- function(x){
  Treatment <- factor(x$trt.ind, levels = c(1,0), labels = c("Experimental Observed", "Control Observed"), ordered = TRUE)

  ggsurv(survfit(x$orig.sv ~ Treatment)) +
    theme_bw() +
    coord_cartesian(ylim = c(0,1))
}

#' SurvCF.KMPlot
#'
#' KM plot of a SurvCF object using GGally and some default settings.
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
#' SurvExt.KMPlot(x)

SurvCF.KMPlot <- function(x){
  Treatment <- factor(x$trt.ind, levels = c(1,0), labels = c("Experimental", "Control Unadjusted"), ordered = TRUE)

  ggsurv(survfit(x$orig.sv ~ Treatment)) +
    theme_bw() +
    coord_cartesian(ylim = c(0,1))
}
