context("Univariate ACE model")
# Simulate a univariate twin model to test functions with
# Based on Nivard's lavaan tutorial 1/3 in each ACE

# Create basic model
mod <- univTwinModel()
mod <- regress(mod, P~1)
expected <- "P_1 ~ 1\nP_2 ~ 1\nA_1 =~ c(a, a)*P_1\nA_2 =~ c(a, a)*P_2\nA_1 ~~ c(1, 0.5)*A_2\nC_1 =~ c(c, c)*P_1\nC_2 =~ c(c, c)*P_2\nC_1 ~~ c(1, 1)*C_2\nP_1 ~~ c(e, e)*P_1\nP_2 ~~ c(e, e)*P_2\nV_P := (a^2 + c^2 + e)\nA_P_share := a^2/V_P\nC_P_share := c^2/V_P\nE_P_share := e/V_P\n"
test_that("Expected regression", {
  expect_equal(as.character(mod), expected)
})

getSimUnivTwinData <- function (nMZ=5000, nDZ=7500, long=F, seed=12345) {
  set.seed(seed)
  require(MASS)
  A <- matrix(1, 2, 2)
  C <- matrix(1, 2, 2)
  E <- diag(2)
  corrDZ <- matrix(c(1, .5, .5, 1), 2, 2)
  # MZ pairs
  dtaMZ <- mvrnorm(nMZ, mu=c(0, 0), Sigma = A + C + E)
  dtaMZ <- cbind.data.frame('MZ', dtaMZ)
  colnames(dtaMZ) <- c('zyg', 'P_1', 'P_2')
  # DZ
  dtaDZ <- mvrnorm(nDZ, mu=c(0, 0), Sigma=corrDZ + C + E)
  dtaDZ <- cbind.data.frame('DZ', dtaDZ)
  colnames(dtaDZ) <- c('zyg', 'P_1', 'P_2')
  if (long) {
    dtaMZ$pid <- paste0('MZ', 1:nMZ)
    dtaDZ$pid <- paste0('DZ', 1:nDZ)
    dta <- data.frame(P=dtaMZ$P_1, zyg=dtaMZ$zyg, pid=dtaMZ$pid, stringsAsFactors = F)
    dta <- rbind.data.frame(dta, data.frame(P=dtaMZ$P_2, zyg=dtaMZ$zyg, pid=dtaMZ$pid, stringsAsFactors = F), stringsAsFactors=F)
    dta <- rbind.data.frame(dta, data.frame(P=dtaDZ$P_1, zyg=dtaDZ$zyg, pid=dtaDZ$pid, stringsAsFactors = F), 
                            data.frame(P=dtaDZ$P_2, zyg=dtaDZ$zyg, pid=dtaDZ$pid, stringsAsFactors = F), stringsAsFactors=F)
    colnames(dta) <- c('P','zyg','pid')
  } else {
    dta <- rbind(dtaMZ, dtaDZ)
  }
  return(dta)
}
dtaLong <- getSimUnivTwinData(nMZ=1000, nDZ=1500, long=T)

#summary(twinlm(P~1, data=subset(dtaLong,zyg != 'zyg'),id='pid',zyg='zyg',DZ='DZ'))
library(lavaan)
library(tidySEM)
dta <- getSimUnivTwinData(nMZ=3000, nDZ=4500)
#lavMod <- univTwinMod('P1', 'P2')
#cat(lavMod)
ft1 <- lavaan(as.character(mod), data=dta, group='zyg',group.label=c('MZ','DZ'), std.lv=T, std.ov=F)
A <- subset(parameterestimates(ft1), lhs=='A_P_share')$est
C <- subset(parameterestimates(ft1), lhs=='C_P_share')$est
E <- subset(parameterestimates(ft1), lhs=='E_P_share')$est
tol <- 0.02

#ft2 <- lavaan(as.character(mod), data=dta, group='zyg',group.label=c('MZ','DZ'), std.lv=T, std.ov=F)
test_that('Expected variance for A=1/3', {
          expect_equal(1/3, A + tol, A - tol)
          expect_equal(1/3, C + tol, C - tol)
          expect_equal(1/3, E + tol, E - tol)
          })
