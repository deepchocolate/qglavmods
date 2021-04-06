#' Generate a Cholesky decomposition for a multivariate twin model
#'
#' @param measures A list of measurments and regressions, eg list(M1=list(1, female))
#' @param suffixGroup A vector of measurment suffix for MZ and DZ group.
#' @export
choleskyMod <- function (measures, suffixGroup, mod=list(A='a',C='c',E='e')) {
  m <- ''
  for (meas in names(measures)) {
    # No direct correlation between twin measures
    measCorr <- paste0(meas,suffixGroup, collapse=' ~~ 0*')
    m <- paste0(m, measCorr,'\n')
    # Regressions
    regr <- paste0(meas, suffixGroup, collapse=' + ')
    regr <- paste0(regr, ' ~ ')
    for (var in measures[[meas]]) {
      lab <- paste0('param_', var)
      lab <- paste0('c(',lab,',',lab,')*',var)
      lab <- paste0(lab, suffixGroup,collapse=' + ' )
      regr <- paste0(regr,' + ', lab)
    }
    m <- paste0(m, regr,'\n')
  }
  # Store parameter names by measure
  parameters <- list()
  # Latent variables
  for (latVar in names(mod)) {
    #parameters[[latVar]] <- list()
    nMeasures <- length(measures)
    for (suffix in suffixGroup) {
      paramBaseLab <- paste0('param_',mod[[latVar]])
      maxMeas <- 1
      for (measLat in names(measures)) {
        latVarName <- paste0(latVar,suffix, '_', measLat)
        latReg <- paste0(latVarName,' =~ ')
        i <- 1
        loadings <- c()
        while (i <= maxMeas) {
          latCovs <- ''
          meas <- names(measures)[i]
          paramLab <- paste0(paramBaseLab, '_',meas)
          if (i != maxMeas) paramLab <- paste0(paramLab, '_', names(measures)[i+1])
          print(paramLab)
          parameters[[measLat]] <- c(parameters[[measLat]], paramLab)
          loadings <- c(loadings, paste0('c(',paramLab,',', paramLab, ')*',meas))
          i <- i + 1
        }
        latReg <- paste0(latReg, paste(loadings, collapse=' + '))
        m <- paste0(m, latReg, '\n')
        maxMeas <- maxMeas + 1
      }
    }
    # Add up covariance restrictions
    for (measLat in names(measures)) {
      latVars <- c()
      for (suffix in suffixGroup) {
        latVarName <- paste0(latVar,suffix, '_', measLat)
        latVars <- c(latVars, latVarName)
      }
      if (latVar %in% c('A','C')) {
        if (latVar == 'A') latCov <- paste(latVars,collapse = '~~ c(1, 0.5)*')
        else latCov <- paste(latVars,collapse = '~~ c(1, 1)*')
        m <- paste0(m, latCov, "\n")
      }
    }
  }
  # Parameter definitions
  # Variances
  for (measure in names(measures)) {
    params <- paste0(parameters[[measure]], collapse=')^2 + (')
    m <- paste0(m, 'var_', measure, ' := (',params, ")^2\n")
  }
  return(m)
}
lmod <- choleskyMod(list(M1=list(1, 'female'), M2=list(1, 'female')), c('_1','_2'))
cat(lmod)
myChol3 <- '
var_adh := (a_Adh)^2 + (c_Adh)^2 + (e_Adh)^2
sh_A_Adh := (a_Adh)^2/var_adh
sh_C_Adh := (c_Adh)^2/var_adh
sh_E_Adh := (e_Adh)^2/var_adh
var_gpa := (a_GPA)^2 + (c_GPA)^2 + (e_GPA)^2 + (a_Adh_GPA)^2 + (c_Adh_GPA)^2 + (e_Adh_GPA)^2
sh_A_gpa := (a_GPA)^2/var_gpa
sh_C_gpa := (c_GPA)^2/var_gpa
sh_E_gpa := (e_GPA)^2/var_gpa
sh_A_Adh_gpa := (a_Adh_GPA)^2/var_gpa
sh_C_Adh_gpa := (c_Adh_GPA)^2/var_gpa
sh_E_Adh_gpa := (e_Adh_GPA)^2/var_gpa
var_ump := (a_Adh_mUmp)^2 + (c_Adh_mUmp)^2 + (e_Adh_mUmp)^2 +(a_GPA_mUmp)^2 + (c_GPA_mUmp)^2 + (e_GPA_mUmp)^2 + (a_mUmp)^2 + (c_mUmp)^2 + (e_mUmp)^2
sh_A_ump := (a_mUmp)^2/var_ump
sh_C_ump := (c_mUmp)^2/var_ump
sh_E_ump := (e_mUmp)^2/var_ump
sh_A_adh_ump := (a_Adh_mUmp)^2/var_ump
sh_C_adh_ump := (c_Adh_mUmp)^2/var_ump
sh_E_adh_ump := (e_Adh_mUmp)^2/var_ump
sh_A_gpa_ump := (a_GPA_mUmp)^2/var_ump
sh_C_gpa_ump := (c_GPA_mUmp)^2/var_ump
sh_E_gpa_ump := (e_GPA_mUmp)^2/var_ump

rg_adh_gpa := a_Adh_GPA/(a_GPA*a_Adh)^.5
rc_adh_gpa := c_Adh_GPA/(c_GPA*c_Adh)^.5
re_adh_gpa := e_Adh_GPA/(e_GPA*a_Adh)^.5
rg_adh_ump := a_Adh_mUmp/(a_mUmp*a_Adh)^.5
rc_adh_ump := c_Adh_mUmp/(c_mUmp*c_Adh)^.5
re_adh_ump := e_Adh_mUmp/(e_mUmp*e_Adh)^.5
rg_gpa_ump := a_GPA_mUmp/(a_GPA*a_mUmp)^.5
rc_gpa_ump := c_GPA_mUmp/(c_GPA*c_mUmp)^.5
re_gpa_ump := e_GPA_mUmp/(e_GPA*e_mUmp)^.5

# Additive genetics
# ADHD symptoms
A1_Adh =~ c(a_Adh, a_Adh)*adh_1 + c(a_Adh_GPA, a_Adh_GPA)*GPA_1 + c(a_Adh_mUmp, a_Adh_mUmp)*mUmp_1
A2_Adh =~ c(a_Adh, a_Adh)*adh_2 + c(a_Adh_GPA, a_Adh_GPA)*GPA_2 + c(a_Adh_mUmp, a_Adh_mUmp)*mUmp_1          # A in twin 2
# Note that the order of parameters matters depending on how lavaan has ordered them
A1_Adh ~~ c(0.5, 1)*A2_Adh         # The correlation between twins is 1 in MZ and .5 in DZ

# GPA
A1_GPA =~ c(a_GPA, a_GPA)*GPA_1 + c(a_GPA_mUmp, a_GPA_mUmp)*mUmp_1
A2_GPA =~ c(a_GPA, a_GPA)*GPA_2 + c(a_GPA_mUmp, a_GPA_mUmp)*mUmp_2          # A in twin 2
# Note that the order of parameters matters depending on how lavaan has ordered them
A1_GPA ~~ c(0.5, 1)*A2_GPA         # The correlation between twins is 1 in MZ and .5 in DZ

# Unemployment
A1_mUmp =~ c(a_mUmp, a_mUmp)*mUmp_1
A2_mUmp =~ c(a_mUmp, a_mUmp)*mUmp_2          # A in twin 2
# Note that the order of parameters matters depending on how lavaan has ordered them
A1_mUmp ~~ c(0.5, 1)*A2_mUmp         # The correlation between twins is 1 in MZ and .5 in DZ

# Shared environment
# ADHD
C1_Adh =~ c(c_Adh, c_Adh)*adh_1 + c(c_Adh_GPA, c_Adh_GPA)*GPA_1 + c(c_Adh_mUmp, c_Adh_mUmp)*mUmp_1
C2_Adh =~ c(c_Adh, c_Adh)*adh_2 + c(c_Adh_GPA, c_Adh_GPA)*GPA_2 + c(c_Adh_mUmp, c_Adh_mUmp)*mUmp_2          # E/Residual variation in twin 2
C1_Adh ~~ c(1,1)*C2_Adh                # Assumption of no correlation across twins
# GPA
C1_GPA =~ c(c_GPA, c_GPA)*GPA_1 + c(c_GPA_mUmp, c_GPA_mUmp)*mUmp_1
C2_GPA =~ c(c_GPA, c_GPA)*GPA_2 + c(c_GPA_mUmp, c_GPA_mUmp)*mUmp_1          # E/Residual variation in twin 2
C1_GPA ~~ c(1,1)*C2_GPA                # Assumption of no correlation across twins
# Unemployment
C1_mUmp =~ c(c_mUmp, c_mUmp)*mUmp_1
C2_mUmp =~ c(c_mUmp, c_mUmp)*mUmp_2
C1_mUmp ~~ c(1,1)*C2_mUmp


# E: Unique environment (DC,MZ)
# ADHD
E1_Adh =~ c(e_Adh, e_Adh)*adh_1 + c(e_Adh_GPA, e_Adh_GPA)*GPA_1 + c(e_Adh_mUmp, e_Adh_mUmp)*mUmp_1
E2_Adh =~ c(e_Adh, e_Adh)*adh_2 + c(e_Adh_GPA, e_Adh_GPA)*GPA_2 + c(e_Adh_mUmp, e_Adh_mUmp)*mUmp_2          # E/Residual variation in twin 2
#E1_Adh ~~ 0*E2_Adh                # Assumption of no correlation across twins
# GPA
E1_GPA =~ c(e_GPA, e_GPA)*GPA_1 + c(e_GPA_mUmp, e_GPA_mUmp)*mUmp_1
E2_GPA =~ c(e_GPA, e_GPA)*GPA_2 + c(e_GPA_mUmp, e_GPA_mUmp)*mUmp_2          # E/Residual variation in twin 2
#E1_GPA ~~ 0*E2_GPA                # Assumption of no correlation across twins
# Unemployment
E1_mUmp =~ c(e_mUmp, e_mUmp)*mUmp_1
E2_mUmp =~ c(e_mUmp, e_mUmp)*mUmp_2
#E1_mUmp ~~ 0*E2_mUmp

#A1_Adh ~~ 0*C1_Adh + 0*C2_Adh + 0*A1_GPA + 0*A2_GPA + 0*A1_mUmp + 0*A2_mUmp + 0*E1_Adh + 0*E2_Adh + 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp +0*E2_mUmp + 0*C1_GPA + 0*C2_GPA + 0*C1_mUmp + 0*C2_mUmp
#A2_Adh ~~ 0*C1_Adh + 0*C2_Adh + 0*A1_GPA + 0*A2_GPA + 0*A1_mUmp + 0*A2_mUmp + 0*E1_Adh + 0*E2_Adh + 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp +0*E2_mUmp + 0*C1_GPA + 0*C2_GPA + 0*C1_mUmp + 0*C2_mUmp
#A1_GPA ~~ 0*C1_Adh + 0*C2_Adh + 0*A1_mUmp + 0*A2_mUmp + 0*E1_Adh + 0*E2_Adh + 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp +0*E2_mUmp + 0*C1_GPA + 0*C2_GPA + 0*C1_mUmp + 0*C2_mUmp
#A2_GPA ~~ 0*C1_Adh + 0*C2_Adh + 0*A1_mUmp + 0*A2_mUmp + 0*E1_Adh + 0*E2_Adh + 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp +0*E2_mUmp + 0*C1_GPA + 0*C2_GPA + 0*C1_mUmp + 0*C2_mUmp
#A1_mUmp ~~ 0*C1_Adh + 0*C2_Adh + 0*E1_Adh + 0*E2_Adh + 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp +0*E2_mUmp + 0*C1_GPA + 0*C2_GPA + 0*C1_mUmp + 0*C2_mUmp
#A2_mUmp ~~ 0*C1_Adh + 0*C2_Adh + 0*E1_Adh + 0*E2_Adh + 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp +0*E2_mUmp + 0*C1_GPA + 0*C2_GPA + 0*C1_mUmp + 0*C2_mUmp
#E1_Adh ~~ 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp + 0*E2_mUmp
#E2_Adh ~~ 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp + 0*E2_mUmp
#E1_GPA ~~ 0*E1_mUmp + 0*E2_mUmp
#E2_GPA ~~ 0*E1_mUmp + 0*E2_mUmp

#C1_Adh ~~ 0*C1_GPA + 0*C2_GPA + 0*C1_mUmp + 0*C2_mUmp + 0*E1_Adh + 0*E2_Adh + 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp + 0*E2_mUmp
#C2_Adh ~~ 0*C1_GPA + 0*C2_GPA + 0*C1_mUmp + 0*C2_mUmp + 0*E1_Adh + 0*E2_Adh + 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp + 0*E2_mUmp
#C1_GPA ~~ 0*C1_mUmp + 0*C2_mUmp + 0*E1_Adh + 0*E2_Adh + 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp + 0*E2_mUmp
#C2_GPA ~~ 0*C1_mUmp + 0*C2_mUmp + 0*E1_Adh + 0*E2_Adh + 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp + 0*E2_mUmp
#C1_mUmp ~~ 0*E1_Adh + 0*E2_Adh + 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp + 0*E2_mUmp
#C2_mUmp ~~ 0*E1_Adh + 0*E2_Adh + 0*E1_GPA + 0*E2_GPA + 0*E1_mUmp + 0*E2_mUmp

#
#adh_1 + adh_2 ~ c(A1_adh,A1_adh)*A_Adh + c(E1_adh,E1_adh)*E_Adh #+ c(int1,int1)*1 + c(fem1,fem1)*female_1 + c(fem1,fem1)*female_2
#GPA_1 + GPA_2 ~ c(A2_adh,A2_adh)*A_Adh + c(E2_adh,E2_adh)*E_Adh + c(A1_gpa,A1_gpa)*A_GPA + c(E1_gpa,E1_gpa)*E_GPA #+ c(int2,in2)*1 + c(fem2,fem2)*female_1 + c(fem2,fem2)*female_2  + c(t11,t11)*year_1 + c(t11,t11)*year_2
#mUmp_1 + mUmp_2 ~ c(A3_adh,A3_adh)*A_Adh + c(E3_adh,E3_adh)*E_Adh + c(A2_gpa,A2_gpa)*A_GPA + c(E2_gpa,E2_gpa)*E_GPA + c(A_mump,A_mump)*A_mUmp + c(E_mump,E_mump)*E_mUmp #+ c(int3,int3)*1 + c(fem3,fem3)*female_1 + c(fem3,fem3)*female_2  + c(t22,t22)*year_1 + c(t22,t22)*year_2
#adh_1 + adh_2 ~ c(A1,A1)*A1_Adh + c(A1,A1)*A2_Adh + c(E1,E1)*E1_Adh + c(E1,E1)*E2_Adh
#GPA_1 + GPA_2 ~ c(A2,A2)*A1_GPA + c(A2,A2)*A2_GPA + c(E2,E2)*E1_GPA + c(E2,E2)*E2_GPA + c(A12,A12)*A1_Adh + c(A12,A12)*A2_Adh + c(E12,E12)*E1_Adh + c(E12,E12)*E2_Adh
#mUmp_1 + mUmp_2 ~ c(A23,A23)*A1_GPA + c(A23,A23)*A2_GPA + c(E23,E23)*E1_GPA + c(E23,E23)*E2_GPA + c(A13,A13)*A1_Adh + c(A13,A13)*A2_Adh + c(E13,E13)*E1_Adh + c(E13,E13)*E2_Adh + c(A3,A3)*A1_mUmp + c(A3,A3)*A2_mUmp + c(E3,E3)*E1_mUmp + c(E3,E3)*E2_mUmp

# Definitions
var_adh := (a_Adh)^2 + (c_Adh)^2 + (e_Adh)^2
sh_A_Adh := (a_Adh)^2/var_adh
sh_C_Adh := (c_Adh)^2/var_adh
sh_E_Adh := (e_Adh)^2/var_adh
var_gpa := (a_GPA)^2 + (c_GPA)^2 + (e_GPA)^2 + (a_Adh_GPA)^2 + (c_Adh_GPA)^2 + (e_Adh_GPA)^2
sh_A_gpa := (a_GPA)^2/var_gpa
sh_C_gpa := (c_GPA)^2/var_gpa
sh_E_gpa := (e_GPA)^2/var_gpa
sh_A_Adh_gpa := (a_Adh_GPA)^2/var_gpa
sh_C_Adh_gpa := (c_Adh_GPA)^2/var_gpa
sh_E_Adh_gpa := (e_Adh_GPA)^2/var_gpa
var_ump := (a_Adh_mUmp)^2 + (c_Adh_mUmp)^2 + (e_Adh_mUmp)^2 +(a_GPA_mUmp)^2 + (c_GPA_mUmp)^2 + (e_GPA_mUmp)^2 + (a_mUmp)^2 + (c_mUmp)^2 + (e_mUmp)^2
sh_A_ump := (a_mUmp)^2/var_ump
sh_C_ump := (c_mUmp)^2/var_ump
sh_E_ump := (e_mUmp)^2/var_ump
sh_A_adh_ump := (a_Adh_mUmp)^2/var_ump
sh_C_adh_ump := (c_Adh_mUmp)^2/var_ump
sh_E_adh_ump := (e_Adh_mUmp)^2/var_ump
sh_A_gpa_ump := (a_GPA_mUmp)^2/var_ump
sh_C_gpa_ump := (c_GPA_mUmp)^2/var_ump
sh_E_gpa_ump := (e_GPA_mUmp)^2/var_ump

rg_adh_gpa := a_Adh_GPA/(a_GPA*a_Adh)^.5
rc_adh_gpa := c_Adh_GPA/(c_GPA*c_Adh)^.5
re_adh_gpa := e_Adh_GPA/(e_GPA*a_Adh)^.5
rg_adh_ump := a_Adh_mUmp/(a_mUmp*a_Adh)^.5
rc_adh_ump := c_Adh_mUmp/(c_mUmp*c_Adh)^.5
re_adh_ump := e_Adh_mUmp/(e_mUmp*e_Adh)^.5
rg_gpa_ump := a_GPA_mUmp/(a_GPA*a_mUmp)^.5
rc_gpa_ump := c_GPA_mUmp/(c_GPA*c_mUmp)^.5
re_gpa_ump := e_GPA_mUmp/(e_GPA*e_mUmp)^.5
'
