#' Significance testing for the 2/3-1/3, 1/3-1/3-1/3, 1/2-1/2 procedures
#'
#' Performs significance testing for the 2/3-1/3, 1/3-1/3-1/3, 1/2-1/2 procedures.
#' Also reports the hazard ratios, 95\% confidence intervals, p-values,
#' nominal significance levels, and correlations for the overall A,
#' simple A, and simple AB test statistics.
#'
#' @param time follow-up times
#' @param event event indicators (0/1)
#' @param indA treatment A indicators (0/1)
#' @param indB treatment B indicators (0/1)
#' @param covmat matrix of covariates; one row per subject
#' @param niter number of interations passed to \code{crit2x2} function call
#'
#' @return \item{hrA }{overall A hazard ratio}
#' @return \item{ciA }{95\% confidence interval for the overall A hazard ratio}
#' @return \item{pvalA }{two-sided p-value for the overall A hazard ratio}
#' @return \item{hra }{simple A hazard ratio}
#' @return \item{cia }{95\% confidence interval for the simple A hazard ratio}
#' @return \item{pvala }{two-sided p-value for the simple A hazard ratio}
#' @return \item{hrab }{simple AB hazard ratio}
#' @return \item{ciab}{95\% confidence interval for the simple AB hazard ratio}
#' @return \item{pvalab }{two-sided p-value for the simple AB hazard ratio}
#' @return \item{sig23A }{2/3-1/3 procedure's p-value rejection criterion for the overall A null hypothesis}
#' @return \item{sig23ab }{2/3-1/3 procedure's p-value rejection criterion for the simple A null hypothesis}
#' @return \item{result23 }{2/3-1/3 procedure's accept/reject for the overall A and simple A null hypotheses results}
#' @return \item{sig13 }{1/3-1/3-1/3 procedure's p-value rejection criterion for the overall A, simple A, and simple AB null hypotheses}
#' @return \item{result13 }{1/3-1/3-1/3 procedure's accept/reject criterion for the overall A, simple A, and simple AB null hypotheses results}
#' @return \item{sig12 }{1/2-1/2 procedure's p-value rejection criterion for the simple A and simple AB null hypotheses}
#' @return \item{result12 }{1/2-1/2 procedure's accept/reject criterion for the simple A and simple AB null hypotheses results}
#' @return \item{corAa }{correlation between the overall A and simple A Wald statistics}
#' @return \item{corAab }{correlation between the overall A and simple AB Wald statistics}
#' @return \item{coraab }{correlation between the simple A and simple AB Wald statistics}
#' @author Eric Leifer, James Troendle
#' @references Leifer, E.S., Troendle, J.F., Kolecki, A., Follmann, D.
#' Joint testing of overall and simple effect for the two-by-two factorial design. (2019). Submitted.
#' @export fac2x2analyze
#' @examples
#'  # First load the simulated data variables. The "simdata" file is
#'  # a 4600-by-9 matrix which is loaded with the Factorial2x2 package.
#'  time <- simdata[, "time"]
#'  event <- simdata[, "event"]
#'  indA <- simdata[, "indA"]
#'  indB <- simdata[, "indB"]
#'  covmat <- simdata[, 6:10]
#' fac2x2analyze(time, event, indA, indB, covmat, niter = 5)
#' # $hrA
#' # [1] 0.8895135
#'
#' # $ciA
#' # [1] 0.786823 1.005607
#'
#' # $pvalA
#' # [1] 0.06139083
#'
#' # $hra
#' # [1] 0.8096082
#'
#' # $cia
#' # [1] 0.6832791 0.9592939
#'
#' # $pvala
#' # [1] 0.01468184
#'
#' # $hrab
#' # [1] 0.7583061
#'
#' # $ciab
#' # [1] 0.6389355 0.8999785
#'
#' # $pvalab
#' # [1] 0.001545967
#'
#' # $sig23A
#' # [1] 0.03333294
#'
#' # $sig23ab
#' # [1] 0.02560439
#'
#' # $result23
#' # [1] "accept overall A" "reject simple AB"
#'
#' # $sig13
#' # [1] 0.02091308
#'
#' # $result13
#' # [1] "accept overall A" "reject simple A"  "reject simple AB"
#'
#' # $sig12
#' # [1] 0.02665044
#'
#' # $result12
#' # [1] "reject simple A"  "reject simple AB"
#'
#' # $corAa
#' # [1] 0.7274961
#'
#' # $corAab
#' # [1] 0.7164075
#'
#' # $coraab
#' # [1] 0.4572905
#'
fac2x2analyze <- function(time, event, indA, indB, covmat, niter = 5){
  # Performs the Wald significance tests for the 2/3-1/3,
  # 1/3-1/3-1/3, and 1/2-1/2 procedures.  NEED to consider the case
  # when covmat = NULL.  It calls crit2x2 which calls
  # cor2x2.
  # time =  follow-up time
  # event = event indicator: 0=censoring, 1=event
  # indA = treatment A indicator (0 = no, 1 = yes)
  # indB = treatment B indicator (0 = no, 1 = yes)
  # covmat = covariate matrix.  NOTE!! Factor variables have to
  #		use 0/1 dummy variables
  # niter = number of iterations in the crit2x2 call
  # It computes the overall A, simple A, and simple AB p-values for
  # the corresponding Cox model log hazard ratio parameter estimates
  # (where the overall A Cox model is stratified
  # on B) to determine the statistical significance of each effect for
  # each of the three procedures.

  coraux <- cor2x2(time, event, indA, indB, covmat)
  corAa <- coraux$corAa
  corAab <- coraux$corAab
  coraab <- coraux$coraab
  hrA <- coraux$hrA
  ciA <- coraux$ciA
  pvalA <- coraux$pvalA
  hra <- coraux$hra
  cia <- coraux$cia
  pvala <- coraux$pvala
  hrab <- coraux$hrab
  ciab <- coraux$ciab
  pvalab <- coraux$pvalab

  critaux <- crit2x2(corAa, corAab, coraab, niter)
  sig23A <- critaux$sig23A
  sig23ab <- critaux$sig23ab
  sig13 <- critaux$sig13
  sig12 <- critaux$sig12

  # Make accept/reject statements.
  accA <- "accept overall A"
  rejA <- "reject overall A"
  acca <- "accept simple A"
  reja <- "reject simple A"
  accab <- "accept simple AB"
  rejab <- "reject simple AB"

  # Statement for 2/3-1/3 procedure
  if (pvalA <= sig23A) {
    res23A <- rejA
  } else {
    res23A <- accA
  }
  if (pvalab <= sig23ab) {
    res23ab <- rejab
  } else {
    res23ab <- accab
  }
  result23 <- c(res23A, res23ab)

  # Statement for 1/3-1/3-1/3 procedure
  if (pvalA <= sig13) {
    res13A <- rejA
  } else {
    res13A <- accA
  }
  if (pvala <= sig13) {
    res13a <- reja
  } else {
    res13a <- acca
  }
  if (pvalab <= sig13) {
    res13ab <- rejab
  } else {
    res13ab <- accab
  }
  result13 <- c(res13A, res13a, res13ab)

  # Statement for 1/2-1/2 procedure
  if (pvala <= sig12) {
    res12a <- reja
  } else {
    res12a <- acca
  }
  if (pvalab <= sig12) {
    res12ab <- rejab
  } else {
    res12ab <- accab
  }
  result12 <- c(res12a, res12ab)

  list(hrA = hrA, ciA = ciA, pvalA = pvalA,
       hra = hra, cia = cia, pvala = pvala,
       hrab = hrab, ciab = ciab, pvalab = pvalab,
       sig23A = sig23A, sig23ab = sig23ab, result23 = result23,
       sig13 = sig13, result13 = result13,
       sig12 = sig12, result12 = result12,
       corAa = corAa, corAab = corAab, coraab = coraab
       )
}

