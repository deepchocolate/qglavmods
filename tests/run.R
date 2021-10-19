# Build the twinModel
devtools::document()
mod <- univTwinModel()
mod <- setSuffix(mod, '1','2')
mod <- regress(mod, P~1)
mod

ft <- lavaan(as.character(mod), data=dta, group='zyg',group.label=c('MZ','DZ'), std.lv=T, std.ov=F)
m <- 'P1 ~ c(int,int)*1
P2 ~ c(int,int)*1
A1 =~ c(a, a)*P1
A2 =~ c(a, a)*P2
A1 ~~ c(1, 0.5)*A2
C1 =~ c(c, c)*P1
C2 =~ c(c, c)*P2
C1 ~~ c(1, 1)*C2
P1 ~~ c(e, e)*P2
'
devtools::document()
mod1 <- univTwinModel()
mod1 <- setLatentFactors(mod1, A='aa', C='bb', E='ee')
mod1 <- setSuffix(mod1, '1','2')
mod1 <- regress(mod1, P~1)
cat(as.character(mod1))
mod2 <- regress(mod1, PA~1)
mod2 <- bivTwinModel(mod1, mod2)
cat(as.character(mod2))
