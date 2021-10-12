# Get one of the ACE matrixes
getMat <- function (ft, measures, factor='A') {
  fitm <- parameterTable(ft)
  mat <- paste0(factor,'_')
  labs <- c(paste0(mat,measures[1]), paste0(mat,measures[1],'_',measures[2]), paste0(mat,measures[1],'_',measures[3]))
  params <- subset(fitm, label %in% labs)$est
  labs <- c(paste0(mat,measures[2]), paste0(mat,measures[2],'_',measures[3]))
  params <- c(params, 0, subset(fitm, label %in% labs)$est)
  labs <- c(paste0(mat,measures[3]))
  params <- c(params, 0, 0, subset(fitm, label %in% labs)$est)
  mtrx <- matrix(params, ncol=3)
  rownames(mtrx) <- measures
  colnames(mtrx) <- c(paste0(mat,measures[1]),paste0(mat,measures[2]),paste0(mat,measures[3]))
  return(mtrx)
}

# Get the variance explained
getVarExp <- function (ft, measures, factor='A') {
  mat <- getMat(ft, measures, factor)
  mat <-  mat %*% t(mat)
  A <- getMat(ft, measures, 'A')
  C <- getMat(ft, measures, 'C')
  E <- getMat(ft, measures, 'E')
  V <- A %*% t(A) + C %*% t(C) + E %*% t(E)
  mat <- mat/V
  rownames(mat) <- measures
  colnames(mat) <- measures
  return(mat)
}

### TEST ORDER
# Pass measurements in an indexed list, variable names as indexes and regressions as values, eg:
# list(var1="c(int,int)*1 + ...", var2="c(int,int)*1 + ...", var3="c(int,int)*1 + ...")
trivCholLavMod <- function (measures) {
  chol <- '
  # Regressions
  {REGRESSIONS}
  
  # Additive genetics
  # {M1}
  A1_{M1} =~ c(a_{M1}, a_{M1})*{M1}_1
  A2_{M1} =~ c(a_{M1}, a_{M1})*{M1}_2          # A in twin 2
  # Note that the order of parameters matters depending on how lavaan has ordered them
  A1_{M1} ~~ c(0.5, 1)*A2_{M1}         # The correlation between twins is 1 in MZ and .5 in DZ
  
  # {M2}
  A1_{M2} =~ c(a_{M2}, a_{M2})*{M2}_1 + c(a_{M1}_{M2}, a_{M1}_{M2})*{M1}_1
  A2_{M2} =~ c(a_{M2}, a_{M2})*{M2}_2 + c(a_{M1}_{M2}, a_{M1}_{M2})*{M1}_2          # A in twin 2
  # Note that the order of parameters matters depending on how lavaan has ordered them
  A1_{M2} ~~ c(0.5, 1)*A2_{M2}         # The correlation between twins is 1 in MZ and .5 in DZ
  
  # {M3}
  A1_{M3} =~ c(a_{M3}, a_{M3})*{M3}_1 +  c(a_{M2}_{M3}, a_{M2}_{M3})*{M2}_1 + c(a_{M1}_{M3}, a_{M1}_{M3})*{M1}_1
  A2_{M3} =~ c(a_{M3}, a_{M3})*{M3}_2 + c(a_{M2}_{M3}, a_{M2}_{M3})*{M2}_2 + c(a_{M1}_{M3}, a_{M1}_{M3})*{M1}_2          # A in twin 2
  # Note that the order of parameters matters depending on how lavaan has ordered them
  A1_{M3} ~~ c(0.5, 1)*A2_{M3}         # The correlation between twins is 1 in MZ and .5 in DZ
  
  # Shared environment
  # {M1}
  C1_{M1} =~ c(c_{M1}, c_{M1})*{M1}_1
  C2_{M1} =~ c(c_{M1}, c_{M1})*{M1}_2         # E/Residual variation in twin 2
  C1_{M1} ~~ c(1,1)*C2_{M1}                # Assumption of no correlation across twins
  # {M2}
  C1_{M2} =~ c(c_{M2}, c_{M2})*{M2}_1 + c(c_{M1}_{M2}, c_{M1}_{M2})*{M1}_1
  C2_{M2} =~ c(c_{M2}, c_{M2})*{M2}_2 + c(c_{M1}_{M2}, c_{M1}_{M2})*{M1}_2          # E/Residual variation in twin 2
  C1_{M2} ~~ c(1,1)*C2_{M2}                # Assumption of no correlation across twins
  # {M3}
  C1_{M3} =~ c(c_{M3}, c_{M3})*{M3}_1 + c(c_{M2}_{M3}, c_{M2}_{M3})*{M2}_1 + c(c_{M1}_{M3}, c_{M1}_{M3})*{M1}_1
  C2_{M3} =~ c(c_{M3}, c_{M3})*{M3}_2 + c(c_{M2}_{M3}, c_{M2}_{M3})*{M2}_2 + c(c_{M1}_{M3}, c_{M1}_{M3})*{M1}_2
  C1_{M3} ~~ c(1,1)*C2_{M3}
  
  
  # E: Unique environment (DC,MZ)
  # {M1}
  E1_{M1} =~ c(e_{M1}, e_{M1})*{M1}_1
  E2_{M1} =~ c(e_{M1}, e_{M1})*{M1}_2
  
  # {M2}
  E1_{M2} =~ c(e_{M2}, e_{M2})*{M2}_1 + c(e_{M1}_{M2}, e_{M1}_{M2})*{M1}_1
  E2_{M2} =~ c(e_{M2}, e_{M2})*{M2}_2 + c(e_{M1}_{M2}, e_{M1}_{M2})*{M1}_2          # E/Residual variation in twin 2
  
  # {M3}
  E1_{M3} =~ c(e_{M3}, e_{M3})*{M3}_1 + c(e_{M2}_{M3}, e_{M2}_{M3})*{M2}_1 + c(e_{M1}_{M3}, e_{M1}_{M3})*{M1}_1
  E2_{M3} =~ c(e_{M3}, e_{M3})*{M3}_2 + c(e_{M2}_{M3}, e_{M2}_{M3})*{M2}_2 + c(e_{M1}_{M3}, e_{M1}_{M3})*{M1}_2
  #E1_mUmp ~~ 0*E2_mUmp
  
  
  # Definitions
  # These definitions are to reduce the parameters to 1 for each path (labels above generate two per groups)
  A_{M1} := a_{M1}
  A_{M1}_{M2} := a_{M1}_{M2}
  A_{M1}_{M3} := a_{M1}_{M3}
  A_{M2} := a_{M2}
  A_{M2}_{M3} := a_{M2}_{M3}
  A_{M3} := a_{M3}
  C_{M1} := c_{M1}
  C_{M1}_{M2} := c_{M1}_{M2}
  C_{M1}_{M3} := c_{M1}_{M3}
  C_{M2} := c_{M2}
  C_{M2}_{M3} := c_{M2}_{M3}
  C_{M3} := c_{M3}
  E_{M1} := e_{M1}
  E_{M1}_{M2} := e_{M1}_{M2}
  E_{M1}_{M3} := e_{M1}_{M3}
  E_{M2} := e_{M2}
  E_{M2}_{M3} := e_{M2}_{M3}
  E_{M3} := e_{M3}
  
  var_{M1} := a_{M1}^2 + c_{M1}^2 + e_{M1}^2
  var_A_{M1} := a_{M1}^2
  var_C_{M1} := c_{M1}^2
  var_E_{M1} := e_{M1}^2
  sh_A_{M1} := (a_{M1}^2)/var_{M1}
  sh_C_{M1} := (c_{M1}^2)/var_{M1}
  sh_E_{M1} := (e_{M1}^2)/var_{M1}
  var_{M2} := a_{M2}^2 + c_{M2}^2 + e_{M2}^2 + a_{M1}_{M2}^2 + c_{M1}_{M2}^2+ e_{M1}_{M2}^2
  var_A_{M2} := a_{M1}_{M2}^2 + a_{M2}^2
  var_C_{M2} := c_{M1}_{M2}^2 + c_{M2}^2
  var_E_{M2} := e_{M1}_{M2}^2 + e_{M2}^2
  sh_A_{M2} := (a_{M2}^2)/var_{M2}
  sh_C_{M2} := (c_{M2}^2)/var_{M2} 
  sh_E_{M2} := (e_{M2}^2)/var_{M2}
  sh_A_{M1}_{M2} := (a_{M1}_{M2}^2)/var_{M2}
  sh_C_{M1}_{M2} := (c_{M1}_{M2}^2)/var_{M2}
  sh_E_{M1}_{M2} := (e_{M1}_{M2}^2)/var_{M2}
  var_{M3} := (a_{M1}_{M3})^2 + c_{M1}_{M3}^2 + (e_{M1}_{M3})^2 + (a_{M2}_{M3})^2 + c_{M2}_{M3}^2  + (e_{M2}_{M3})^2 + (a_{M3})^2 + (c_{M3})^2 + (e_{M3})^2
  var_A_{M3} := a_{M1}_{M3}^2 + a_{M2}_{M3}^2 + a_{M3}^2
  var_C_{M3} := c_{M1}_{M3}^2 + c_{M2}_{M3}^2 + c_{M3}^2
  var_E_{M3} := (e_{M1}_{M3})^2  + (e_{M2}_{M3})^2 + (e_{M3})^2
  sh_A_{M3} := var_A_{M3}/var_{M3}
  sh_C_{M3} := var_C_{M3}/var_{M3}
  sh_E_{M3} := var_E_{M3}/var_{M3}
  sh_A_{M1}_{M3} := (a_{M1}_{M3})^2/var_{M3}
  sh_C_{M1}_{M3} := (c_{M1}_{M3})^2/var_{M3}
  sh_E_{M1}_{M3} := (e_{M1}_{M3})^2/var_{M3}
  sh_A_{M2}_{M3} := (a_{M2}_{M3})^2/var_{M3}
  sh_C_{M2}_{M3} := (c_{M2}_{M3})^2/var_{M3}
  sh_E_{M2}_{M3} := (e_{M2}_{M3})^2/var_{M3}
  # Genetic and environmental correlations
  rg_{M1}_{M2} := (a_{M1}*a_{M1}_{M2})/sqrt(var_A_{M1} * var_A_{M2})
  rc_{M1}_{M2} := c_{M1}*c_{M1}_{M2}/sqrt(var_C_{M2}*var_C_{M1})
  re_{M1}_{M2} := e_{M1}*e_{M1}_{M2}/sqrt(var_E_{M2}*var_E_{M1})
  rg_{M1}_{M3} := a_{M1}*a_{M1}_{M3}/sqrt(var_A_{M1}*var_A_{M3})
  rc_{M1}_{M3} := c_{M1}*c_{M1}_{M3}/sqrt(var_C_{M1}*var_C_{M3})
  re_{M1}_{M3} := e_{M1}*e_{M1}_{M3}/sqrt(var_E_{M1}*var_E_{M3})
  rg_{M2}_{M3} := (a_{M1}_{M3}*a_{M1}_{M2} + a_{M2}*a_{M2}_{M3})/sqrt(var_A_{M2}*var_A_{M3})
  rc_{M2}_{M3} := (c_{M1}_{M3}*c_{M1}_{M2} + c_{M2}*c_{M2}_{M3})/sqrt(var_C_{M2}*var_C_{M3})
  re_{M2}_{M3} := (e_{M1}_{M3}*e_{M1}_{M2} + e_{M2}*e_{M2}_{M3})/sqrt(var_E_{M2}*var_E_{M3})
  # Phenotypic correlations
  ph_{M1}_{M2} := rg_{M1}_{M2}*sqrt(sh_A_{M1})*sqrt(sh_A_{M2}) + rc_{M1}_{M2}*sqrt(sh_C_{M1})*sqrt(sh_C_{M2}) + re_{M1}_{M2}*sqrt(sh_E_{M1})*sqrt(sh_E_{M2})
  ph_{M1}_{M3} := rg_{M1}_{M3}*sqrt(sh_A_{M1})*sqrt(sh_A_{M3}) + rc_{M1}_{M3}*sqrt(sh_C_{M1})*sqrt(sh_C_{M3}) + re_{M1}_{M3}*sqrt(sh_E_{M1})*sqrt(sh_E_{M3})
  ph_{M2}_{M3} := rg_{M2}_{M3}*sqrt(sh_A_{M2})*sqrt(sh_A_{M3}) + rc_{M2}_{M3}*sqrt(sh_C_{M2})*sqrt(sh_C_{M3}) + re_{M2}_{M3}*sqrt(sh_E_{M2})*sqrt(sh_E_{M3})
  '
  varNames <- names(measures)
  m1 <- varNames[1]
  m2 <- varNames[2]
  m3 <- varNames[3]
  regressions <- ''
  if (measures[[m1]] != '') regressions <- paste0(regressions, m1,'_1 + ',m1,'_2 ~ ',measures[[m1]],'\n');
  if (measures[[m2]] != '') regressions <- paste0(regressions, m2,'_1 + ',m2,'_2 ~ ',measures[[m2]],'\n');
  if (measures[[m3]] != '') regressions <- paste0(regressions, m3,'_1 + ',m3,'_2 ~ ',measures[[m3]]);
  chol <- gsub(x=chol,pattern = '\\{REGRESSIONS\\}', replacement=regressions)
  chol <- gsub(x=chol, pattern = '\\{M1\\}',replacement = m1)
  chol <- gsub(x=chol, pattern = '\\{M2\\}',replacement = m2)
  chol <- gsub(x=chol, pattern = '\\{M3\\}',replacement = m3)
  return(chol)
}


meas <- list(stdGPA='c(int2,int2)*1 + c(fem2,fem2)*female_1 + c(fem2,fem2)*female_2',
             stdAdh='c(int1,int1)*1 + c(fem1,fem1)*female_1 + c(fem1,fem1)*female_2',
             stdmUmp='c(int3,int3)*1 + c(fem3,fem3)*female_1 + c(fem3,fem3)*female_2')
meas <- list(adh='c(int1,int1)*1 + c(fem1,fem1)*female_1 + c(fem1,fem1)*female_2',
             GPA='c(int2,int2)*1 + c(fem2,fem2)*female_1 + c(fem2,fem2)*female_2',
             mUmp='c(int3,int3)*1 + c(fem3,fem3)*female_1 + c(fem3,fem3)*female_2')
ft3 <- lavaan(trivCholLavMod(meas),data=subset(studentsWide,!is.na(zyg) & Y_1<=2010 & Y_2<=2010),group='zyg',group.label=c('DZ','MZ'),
              missing='listwise',information='observed',std.lv=T,std.ov=F)
summary(ft3,standardized=T)
getMat(ft3, c('stdAdh','stdGPA','stdmUmp'), factor='A')
getVarExp(ft3, c('stdAdh','stdGPA','stdmUmp'),factor='A')

meas <- list(stdGPA='c(fem2,fem2)*female_1 + c(fem2,fem2)*female_2',
             stdAdh='c(fem1,fem1)*female_1 + c(fem1,fem1)*female_2',
             stdmUmp='c(fem3,fem3)*female_1 + c(fem3,fem3)*female_2 + c(y1,y1)*year_1 + c(y1,y1)*year_2')
ft3 <- lavaan(trivCholLavMod(meas),data=subset(studentsWide,!is.na(zyg)),group='zyg',group.label=c('DZ','MZ'),
              missing='listwise',information='observed',std.lv=T,std.ov=F)
getMat(ft3, c('stdGPA','stdAdh','stdmUmp'))
getVarExp(ft3, c('stdGPA','stdAdh','stdmUmp'))
