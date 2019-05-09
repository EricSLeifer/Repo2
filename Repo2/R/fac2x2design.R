#' Power for the 2/3-1/3, 1/3-1/3-1/3, 1/2-1/2 procedures
#'
#' @param n  total sample size
#' @param rateC group C one year event rate
#' @param hrA   group A to group C hazard ratio
#' @param hrB  group B to group C hazard ratio
#' @param hrAB  group AB to group C hazard ratio
#' @param mincens  minimum censoring time
#' @param maxcens  maximum censoring time
#' @param dig number of decimal places to \code{\link{roundDown}} the critical value to
#' @param alpha two-sided significance level
#' @param corAa	correlation between the overall A and simple A log hazard ratio estimates
#' @param corAab correlation between the overall A and simple AB log hazard ratio estimates
#' @param coraab	correalation between the simple A and simple AB log hazard ratio estimates
#' @return \item{powerA}{power to detect the overall A effect at the two-sided \code{alpha} level}
#' @return \item{powerB}{power to detect the overall B effect at the two-sided \code{alpha} level}
#' @return \item{power23.13 }{power to detect the overall A or simple AB effects using the 2/3-1/3 procedure}
#' @return \item{power13.13.13 }{power to detect the overall A, simple A,  or simple AB effects using
#' the 1/3-1/3-1/3 procedure}
#' @return \item{power12.12 }{power to detect the simple A or simple AB effects using
#' the 1/2-1/2 procedure}
#' @return \item{events}{expected number of events}
#' @return \item{evtprob}{event probabilities for the C, A, B, and AB groups, respectively}
#' @seealso  \code{\link{eventProb}}, \code{\link{crit2x2}}, \code{\link{strLgrkPower}}
#' \code{\link{power23_13}}, \code{\link{power13_13_13}}, \code{\link{power12_12}}
#' @references Leifer, E.S., Troendle, J.F., Kolecki, A., Follmann, D.
#' Joint testing of overall and simple effect for the two-by-two factorial design. (2019). Submitted.
#' @references Slud, E.V. Analysis of factorial survival experiments. Biometrics. 1994; 50: 25-38.
#' @export fac2x2design
#' @examples
#' # Corresponds to scenario 5 in Table 2 from Leifer, Troendle, et al. (2019)
#' n <- 4600
#' rateC <- 0.0445
#' hrA <- 0.80
#' hrB <- 0.80
#' hrAB <- 0.72
#' mincens <- 4.0
#' maxcens <- 8.4
#'
#' fac2x2design(n, rateC, hrA, hrB, hrAB, mincens, maxcens, dig = 2, alpha = 0.05)
#' # $powerA
#' # [1] 0.7182932
#'
#' # $power23.13
#' # [1] 0.9290271
#'
#' # $power13.13.13
#' # [1] 0.9302084
#'
#' # $power12.12
#' # [1] 0.9411688
#'
fac2x2design <- function(n, rateC, hrA, hrB, hrAB, mincens, maxcens, dig = 2, alpha = 0.05,
                         corAa = 1/sqrt(2), corAab = 1/sqrt(2), coraab = 1/2){
  evtprob <- eventProb(rateC, hrA, hrB, hrAB, mincens, maxcens)
  avgprob <- evtprob$avgprob
  probA_C <- evtprob$probA_C
  probAB_C <- evtprob$probAB_C
  critvals <- crit2x2(corAa, corAab, coraab, dig, alpha)
  crit23A <- critvals$crit23A
  crit23ab <- critvals$crit23ab
  crit13 <- critvals$crit13
  crit12 <- critvals$crit12

  # compute power for the overall A test
  powerA <- strLgrkPower(n, hrA, hrB, hrAB, avgprob, dig, alpha)$power

  # compute the power for the overall B test be reversing the roles of
  # of hrA and hrB in the previous line
  powerB <- strLgrkPower(n, hrB, hrA, hrAB, avgprob, dig, alpha)$power

  power23.13 <- power23_13(n, hrA, hrB, hrAB, avgprob, probAB_C,
    crit23A, crit23ab, dig, cormat =
    matrix(c(1, sqrt(0.5), sqrt(0.5), 1), byrow = T,
    nrow = 2), niter = 100, rseed = 047477, abseps = 1e-05)$power23.13

  power13.13.13 <- power13_13_13(n, hrA, hrB, hrAB, dig, avgprob, probA_C,
    probAB_C, crit13, cormat12 = matrix(c(1, sqrt(0.5),
    sqrt(0.5), 1), byrow = T, nrow = 2), cormat23 = matrix(c(1, 0.5, 0.5, 1),
    byrow = T, nrow = 2),
    cormat123 = matrix(c(1, sqrt(0.5), sqrt(0.5), sqrt(0.5), 1, 0.5,
    sqrt(0.5), 0.5, 1), byrow=T, nrow = 3), niter = 100, rseed = 047477,
    abseps = 1e-05)$power13.13.13

  power12.12 <- power12_12(n, hrA, hrAB, probA_C, probAB_C,
    crit12, cormat = matrix(c(1,0.5,0.5,1), byrow =T, nrow =2),
    niter = 100, rseed = 047477, abseps = 1e-05)$power12.12

  list(powerA = powerA, powerB = powerB, power23.13 = power23.13,
       power13.13.13 = power13.13.13,
       power12.12 = power12.12, events = n * avgprob,
       evtprob = c(unlist(evtprob))[2:5])
}
