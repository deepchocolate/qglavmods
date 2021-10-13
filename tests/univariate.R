# Simulate a univariate twin model to test functions with
# Based on Nivard's lavaan tutorial 1/3 in each ACE
library(MASS)
library(lavaan)
library(tidySEM)
A <- matrix(1, 2, 2)
C <- matrix(1, 2, 2)
E <- diag(2)
corrDZ <- matrix(c(1, .5, .5, 1), 2, 2)
# MZ pairs
dtaMZ <- mvrnorm(1000, mu=c(0, 0), Sigma = A + C + E)
dtaMZ <- cbind.data.frame('MZ', dtaMZ)
colnames(dtaMZ) <- c('zyg', 'P1', 'P2')
# DZ
dtaDZ <- mvrnorm(1500, mu=c(0, 0), Sigma=corrDZ + C + E)
dtaDZ <- cbind.data.frame('DZ', dtaDZ)
colnames(dtaDZ) <- c('zyg', 'P1', 'P2')
dta <- rbind(dtaMZ, dtaDZ)

lavMod <- '
A1 =~ c(a, a)*P1
A2 =~ c(a, a)*P2
A1 ~~ c(1, 0.5)*A2
C1 =~ c(c, c)*P1
C2 =~ c(c, c)*P2
C1 ~~ 1*C2
P1 ~~ c(e2, e2)*P1
P2 ~~ c(e2, e2)*P2

add := a^2/(a^2+c^2+e2)
sha := c^2/(a^2+c^2+e2)
uni := e2/(a^2+c^2+e2)
'
ft <- lavaan(lavMod, data=dta, group='zyg', std.lv=T, std.ov=F)
summary(ft)

binDta <- dta
binDta[binDta[,2] < 1, 2] <- 0
binDta[binDta[,2] > 1, 2] <- 1
binDta[binDta[,3] < 1, 3] <- 0
binDta[binDta[,3] > 1, 3] <- 1

lavMod <- '
A1 =~ c(a, a)*P1
A2 =~ c(a, a)*P2
A1 ~~ c(1, 0.5)*A2
C1 =~ c(c, c)*P1
C2 =~ c(c, c)*P2
C1 ~~ 1*C2
P1 ~~ c(e2, e2)*P1
P2 ~~ c(e2, e2)*P2

#add := a^2/(a^2+c^2+e2)
#sha := c^2/(a^2+c^2+e2)
#uni := e2/(a^2+c^2+e2)
P1 | 1*t1
P2 | 1*t1
'
ft <- lavaan(lavMod, data=binDta, group='zyg', std.lv=T, ordered=c('P1','P2'), parameterization='theta',auto.th=F)
summary(ft)
lay <- get_layout("", "", "", "",   "", "", "",
                  "A1", "", "C1", "",   "A2", "", "C2",
                  "", "P1", "", "",   "", "P2", "", rows=3)
graph_sem(model=ft, layout=lay, variance_diameter=.3, angle=180, rect_height=.35, ellipses_height=.35, spacing_x=1, spacing_y=0.7)
