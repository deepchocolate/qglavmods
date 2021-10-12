# measures List of variablename=covariates
# labels List of varname=label
# model List of model for each trait
bivCorrFacLavMod <- function (measures, labels, model=list(A=1,C=1,E=1)) {
  mod <- '# Traits in twins are assumed to be uncorrelated
  # Regressions
  {REGRESSIONS}
  
  # Additive genetics
  A_{L1}_1 =~ c(a_{L1}, a_{L1})*{M1}_1
  A_{L2}_1 =~ c(a_{L2}, a_{L2})*{M2}_1
  A_{L1}_1 ~~ c(1, 0.5)*A_{L1}_2
  # Twin 2 measures
  A_{L1}_2 =~ c(a_{L1}, a_{L1})*{M1}_2
  A_{L2}_2 =~ c(a_{L2}, a_{L2})*{M2}_2
  A_{L2}_1 ~~ c(1, 0.5)*A_{L2}_2
  # Cross-twin cross-trait correlations
  A_{L1}_1 ~~ c(1*rA_{L1}_{L2}, 0.5*rA_{L1}_{L2})*A_{L2}_2
  A_{L2}_1 ~~ c(1*rA_{L1}_{L2}, 0.5*rA_{L1}_{L2})*A_{L1}_2
  # Cross-trait correlations
  A_{L1}_1 ~~ c(rA_{L1}_{L2}, rA_{L1}_{L2})*A_{L2}_1
  A_{L1}_2 ~~ c(rA_{L1}_{L2}, rA_{L1}_{L2})*A_{L2}_2
  # Shared environment
  # Twin 1 measures
  C_{L1}_1 =~ c(c_{L1},c_{L1})*{M1}_1
  C_{L2}_1 =~ c(c_{L2},c_{L2})*{M2}_1
  C_{L1}_1 ~~ c(1, 1)*C_{L1}_2
  # Twin 2 measures
  C_{L1}_2 =~ c(c_{L1}, c_{L1})*{M1}_2
  C_{L2}_2 =~ c(c_{L2}, c_{L2})*{M2}_2
  C_{L2}_1 ~~ c(1, 1)*C_{L2}_2
  # Cross-twin cross-trait correlations
  C_{L1}_1 ~~ c(rC_{L1}_{L2}, rC_{L1}_{L2})*C_{L2}_2
  C_{L2}_1 ~~ c(rC_{L1}_{L2}, rC_{L1}_{L2})*C_{L1}_2
  # Cross-trait correlations
  C_{L1}_1 ~~ c(rC_{L1}_{L2}, rC_{L1}_{L2})*C_{L2}_1
  C_{L1}_2 ~~ c(rC_{L1}_{L2}, rC_{L1}_{L2})*C_{L2}_2
  # Unique environment
  # Twin 1 measures
  E_{L1}_1 =~ c(e_{L1}, e_{L1})*{M1}_1
  E_{L2}_1 =~ c(e_{L2}, e_{L2})*{M2}_1
  # Twin 2 measures
  E_{L1}_2 =~ c(e_{L1}, e_{L1})*{M1}_2
  E_{L2}_2 =~ c(e_{L2}, e_{L2})*{M2}_2
  # Cross-trait correlations
  E_{L1}_1 ~~ c(rE_{L1}_{L2}, rE_{L1}_{L2})*E_{L2}_1
  E_{L1}_2 ~~ c(rE_{L1}_{L2}, rE_{L1}_{L2})*E_{L2}_2
  
  # Defined parameters
  path_a_{L1} := a_{L1}
  path_c_{L1} := c_{L1}
  path_e_{L1} := e_{L1}
  path_a_{L2} := a_{L2}
  path_c_{L2} := c_{L2}
  path_e_{L2} := e_{L2}
  # Variance explained in each trait
  v_{L1} := a_{L1}^2 + c_{L1}^2 + e_{L1}^2
  a_{L1}_share := a_{L1}^2/v_{L1}
  c_{L1}_share := c_{L1}^2/v_{L1}
  e_{L1}_share := e_{L1}^2/v_{L1}
  v_{L2} := a_{L2}^2 + c_{L2}^2 + e_{L2}^2
  a_{L2}_share := a_{L2}^2/v_{L2}
  c_{L2}_share := c_{L2}^2/v_{L2}
  e_{L2}_share := e_{L2}^2/v_{L2}
  # Genetic and environmental correlations between traits
  corr_A := rA_{L1}_{L2}
  corr_C := rC_{L1}_{L2}
  corr_E := rE_{L1}_{L2}
  # Genetic and environmental contributions to phenotypic correlations
  contr_A := rA_{L1}_{L2}*sqrt(a_{L1}_share)*sqrt(a_{L2}_share)
  contr_C := rC_{L1}_{L2}*sqrt(c_{L1}_share)*sqrt(c_{L2}_share)
  contr_E := rE_{L1}_{L2}*sqrt(e_{L1}_share)*sqrt(e_{L2}_share)
  pheno := contr_A + contr_E
  contr_A_share := contr_A/pheno
  contr_C_share := contr_C/pheno
  contr_E_share := contr_E/pheno'
  varNames <- names(measures)
  m1 <- varNames[1]
  m2 <- varNames[2]
  regressions <- ''
  if (measures[[m1]] != '') regressions <- paste0(regressions, m1,'_1 + ',m1,'_2 ~ ',measures[[m1]],'\n');
  if (measures[[m2]] != '') regressions <- paste0(regressions, m2,'_1 + ',m2,'_2 ~ ',measures[[m2]],'\n');
  if (model[['C']]==0) {
    mod <- paste0(mod, '\nc_{L1} == 0\nrC_{L1}_{L2} == 0\n')
  }
  # Insert regressions
  mod <- gsub(x=mod,pattern = '\\{REGRESSIONS\\}', replacement=regressions)
  # Insert measurements
  mod <- gsub(x=mod, pattern = '\\{M1\\}',replacement = m1)
  mod <- gsub(x=mod, pattern = '\\{M2\\}',replacement = m2)
  # Insert labels
  mod <- gsub(x=mod, pattern = '\\{L1\\}',replacement = labels[1])
  mod <- gsub(x=mod, pattern = '\\{L2\\}',replacement = labels[2])
  return(mod)
}

triCorrFacLavMod <- function (measures, labels, model=list(A=1,C=1,E=1)) {
  mod <- '
  # Regressions
  {REGRESSIONS}
  
  # Additive genetics
  # Trait {M1}
  A1_{L1} =~ c(a_{L1}, a_{L1})*{M1}_1          # A in twin 1
  A2_{L1} =~ c(a_{L1}, a_{L1})*{M1}_2          # A in twin 2
  # Note that the order of parameters matters depending on how lavaan has ordered them
  A1_{L1} ~~ c(1,0.5)*A2_{L1}         # The correlation between twins is 1 in MZ and .5 in DZ
  
  # Trait {M2}
  A1_{L2} =~ c(a_{L2}, a_{L2})*{M2}_1          # A in twin 1
  A2_{L2} =~ c(a_{L2}, a_{L2})*{M2}_2          # A in twin 2
  # Note that the order of parameters matters depending on how lavaan has ordered them
  A1_{L2} ~~ c(1,0.5)*A2_{L2}         # The correlation between twins is 1 in MZ and .5 in DZ
  
  # School performance
  A1_{L3} =~ c(a_{L3}, a_{L3})*{M3}_1          # A in twin 1
  A2_{L3} =~ c(a_{L3}, a_{L3})*{M3}_2          # A in twin 2
  # Note that the order of parameters matters depending on how lavaan has ordered them
  A1_{L3} ~~ c(1, 0.5)*A2_{L3}         # The correlation between twins is 1 in MZ and .5 in DZ
  
  # Correlation A1_S <-> A2_Att
  # Cross twin cross trait
  A1_{L1} ~~ c(1*rA_{L1}_{L3}, 0.5*rA_{L1}_{L3})*A2_{L3}
  A1_{L1} ~~ c(1*rA_{L1}_{L2}, 0.5*rA_{L1}_{L2})*A2_{L2}
  A1_{L2} ~~ c(1*rA_{L2}_{L3}, 0.5*rA_{L2}_{L3})*A2_{L3}
  A1_{L2} ~~ c(1*rA_{L1}_{L2}, 0.5*rA_{L1}_{L2})*A2_{L1}
  A1_{L3} ~~ c(1*rA_{L2}_{L3}, 0.5*rA_{L2}_{L3})*A2_{L2}
  A1_{L3} ~~ c(1*rA_{L1}_{L3}, 0.5*rA_{L1}_{L3})*A2_{L1}
  ### Cross trait
  # Twin 1
  A1_{L3} ~~ c(rA_{L2}_{L3}, rA_{L2}_{L3})*A1_{L2}
  A1_{L3} ~~ c(rA_{L1}_{L3}, rA_{L1}_{L3})*A1_{L1}
  A1_{L2} ~~ c(rA_{L1}_{L2}, rA_{L1}_{L2})*A1_{L1}
  # Twin 2
  A2_{L3} ~~ c(rA_{L2}_{L3}, rA_{L2}_{L3})*A2_{L2}
  A2_{L3} ~~ c(rA_{L1}_{L3}, rA_{L1}_{L3})*A2_{L1}
  A2_{L2} ~~ c(rA_{L1}_{L2}, rA_{L1}_{L2})*A2_{L1}
  
  # Shared environment
  # Inattention
  C1_{L1} =~ c(c_{L1}, c_{L1})*{M1}_1
  C2_{L1} =~ c(c_{L1}, c_{L1})*{M1}_2
  C1_{L1} ~~ 1*C2_{L1}                # C/Shared environment is perfectly correlated between twins
  # GPA
  C1_{L2} =~ c(c_{L2}, c_{L2})*{M2}_1
  C2_{L2} =~ c(c_{L2}, c_{L2})*{M2}_2
  C1_{L2} ~~ 1*C2_{L2}                # C/Shared environment is perfectly correlated between twins
  # School
  C1_{L3} =~ c(c_{L3}, c_{L3})*{M3}_1
  C2_{L3} =~ c(c_{L3}, c_{L3})*{M3}_2
  C1_{L3} ~~ 1*C2_{L3}                # C/Shared environment is perfectly correlated between twins
  # Correlation C1<->C2
  # Cross-twin cross trait
  # School Performance
  C1_{L3} ~~ c(rC_{L1}_{L3}, rC_{L1}_{L3})*C2_{L1}
  C1_{L3} ~~ c(rC_{L2}_{L3}, rC_{L2}_{L3})*C2_{L2}
  C1_{L2} ~~ c(rC_{L2}_{L3}, rC_{L2}_{L3})*C2_{L3}
  C1_{L2} ~~ c(rC_{L1}_{L2}, rC_{L1}_{L2})*C2_{L1}
  C1_{L1} ~~ c(rC_{L1}_{L3}, rC_{L1}_{L3})*C2_{L3}
  C1_{L1} ~~ c(rC_{L1}_{L2}, rC_{L1}_{L2})*C2_{L2}
  ### Cross-trait
  # Twin 1
  C1_{L1} ~~ c(rC_{L1}_{L2}, rC_{L1}_{L2})*C1_{L2}
  C1_{L3} ~~ c(rC_{L2}_{L3}, rC_{L2}_{L3})*C1_{L2}
  C1_{L3} ~~ c(rC_{L1}_{L3}, rC_{L1}_{L3})*C1_{L1}
  # Twin 2
  C2_{L1} ~~ c(rC_{L1}_{L2}, rC_{L1}_{L2})*C2_{L2}
  C2_{L3} ~~ c(rC_{L2}_{L3}, rC_{L2}_{L3})*C2_{L2}
  C2_{L3} ~~ c(rC_{L1}_{L3}, rC_{L1}_{L3})*C2_{L1}
  
  # E: Unique environment (DC,MZ)
  # ADHD
  E1_{L1} =~ c(e_{L1}, e_{L1})*{M1}_1          # E/Residual variation in twin 1
  E2_{L1} =~ c(e_{L1}, e_{L1})*{M1}_2          # E/Residual variation in twin 2
  # GPA
  E1_{L2} =~ c(e_{L2}, e_{L2})*{M2}_1          # E/Residual variation in twin 1
  E2_{L2} =~ c(e_{L2}, e_{L2})*{M2}_2          # E/Residual variation in twin 2
  # School performance
  E1_{L3} =~ c(e_{L3}, e_{L3})*{M3}_1          # E/Residual variation in twin 1
  E2_{L3} =~ c(e_{L3}, e_{L3})*{M3}_2          # E/Residual variation in twin 2
  ### Cross traits
  # Twin 1
  E1_{L3} ~~ c(rE_{L1}_{L3}, rE_{L1}_{L3})*E1_{L1}
  E1_{L3} ~~ c(rE_{L2}_{L3}, rE_{L2}_{L3})*E1_{L2}
  E1_{L2} ~~ c(rE_{L1}_{L2}, rE_{L1}_{L2})*E1_{L1}
  # Twin 2
  E2_{L3} ~~ c(rE_{L1}_{L3}, rE_{L1}_{L3})*E2_{L1}
  E2_{L3} ~~ c(rE_{L2}_{L3}, rE_{L2}_{L3})*E2_{L2}
  E2_{L2} ~~ c(rE_{L1}_{L2}, rE_{L1}_{L2})*E2_{L1}
  
  ### Defined paremeters 
  # Var explained by ACE
  v_{L1} := a_{L1}^2 + c_{L1}^2 + e_{L1}^2
  h_{L1} := a_{L1}^2/v_{L1}
  a_{L1}_share := a_{L1}^2/v_{L1}
  c_{L1}_share := c_{L1}^2/v_{L1}
  e_{L1}_share := e_{L1}^2/v_{L1}
  
  v_{L2} := a_{L2}^2 + c_{L2}^2 + e_{L2}^2
  h_{L2} := a_{L2}^2/v_{L2}
  a_{L2}_share := a_{L2}^2/v_{L2}
  c_{L2}_share := c_{L2}^2/v_{L2}
  e_{L2}_share := e_{L2}^2/v_{L2}
  
  v_{L3} := a_{L3}^2 + c_{L3}^2 + e_{L3}^2
  h_{L3} := a_{L3}^2/v_{L3}
  a_{L3}_share := a_{L3}^2/v_{L3}
  c_{L3}_share := c_{L3}^2/v_{L3}
  e_{L3}_share := e_{L3}^2/v_{L3}
  
  # Genetic correlations
  corr_A_{L1}_{L3} := rA_{L1}_{L3}
  corr_A_{L1}_{L2} := rA_{L1}_{L2}
  corr_A_{L2}_{L3} := rA_{L2}_{L3}
  # Shared environmental correlations
  corr_C_{L1}_{L3} := rC_{L1}_{L3}
  corr_C_{L1}_{L2} := rC_{L1}_{L2}
  corr_C_{L2}_{L3} := rC_{L2}_{L3}
  # Unique environmental correlations
  corr_E_{L1}_{L3} := rE_{L1}_{L3}
  corr_E_{L1}_{L2} := rE_{L1}_{L2}
  corr_E_{L2}_{L3} := rE_{L2}_{L3}
  
  # Genetic and environmental contributions to ACE
  contr_{L1}_{L2}_rA := rA_{L1}_{L2}*sqrt(a_{L2}_share) * sqrt(a_{L1}_share)
  contr_{L1}_{L2}_rC := rC_{L1}_{L2}*sqrt(c_{L2}_share) * sqrt(c_{L1}_share)
  contr_{L1}_{L2}_rE := rE_{L1}_{L2}*sqrt(e_{L2}_share) * sqrt(e_{L1}_share)
  contr_{L1}_{L3}_rA := rA_{L1}_{L3}*sqrt(a_{L3}_share) * sqrt(a_{L1}_share)
  contr_{L1}_{L3}_rC := rC_{L1}_{L3}*sqrt(c_{L3}_share) * sqrt(c_{L1}_share)
  contr_{L1}_{L3}_rE := rE_{L1}_{L3}*sqrt(e_{L3}_share) * sqrt(e_{L1}_share)
  contr_{L2}_{L3}_rA := rA_{L2}_{L3}*sqrt(a_{L3}_share) * sqrt(a_{L2}_share)
  contr_{L2}_{L3}_rC := rC_{L2}_{L3}*sqrt(c_{L3}_share) * sqrt(c_{L2}_share)
  contr_{L2}_{L3}_rE := rE_{L2}_{L3}*sqrt(e_{L3}_share) * sqrt(e_{L2}_share)
  # Check that all calculations equal phenotypic correlations
  pheno_{L1}_{L3} := contr_{L1}_{L3}_rA + contr_{L1}_{L3}_rC + contr_{L1}_{L3}_rE 
  pheno_{L1}_{L2} := contr_{L1}_{L2}_rA + contr_{L1}_{L2}_rC + contr_{L1}_{L2}_rE
  pheno_{L2}_{L3} := contr_{L2}_{L3}_rA + contr_{L2}_{L3}_rC + contr_{L2}_{L3}_rE
  # AS a percentage
  rA_perc_{L1}_{L2} := 100*contr_{L1}_{L2}_rA/pheno_{L1}_{L2}
  rC_perc_{L1}_{L2} := 100*contr_{L1}_{L2}_rC/pheno_{L1}_{L2}
  rE_perc_{L1}_{L2} := 100*contr_{L1}_{L2}_rE/pheno_{L1}_{L2}
  rA_perc_{L1}_{L3} := 100*contr_{L1}_{L3}_rA/pheno_{L1}_{L3}
  rC_perc_{L1}_{L3} := 100*contr_{L1}_{L3}_rC/pheno_{L1}_{L3}
  rE_perc_{L1}_{L3} := 100*contr_{L1}_{L3}_rE/pheno_{L1}_{L3}
  rA_perc_{L2}_{L3} := 100*contr_{L2}_{L3}_rA/pheno_{L2}_{L3}
  rC_perc_{L2}_{L3} := 100*contr_{L2}_{L3}_rC/pheno_{L2}_{L3}
  rE_perc_{L2}_{L3} := 100*contr_{L2}_{L3}_rE/pheno_{L2}_{L3}
  '
  varNames <- names(measures)
  m1 <- varNames[1]
  m2 <- varNames[2]
  m3 <- varNames[3]
  regressions <- ''
  if (measures[[m1]] != '') regressions <- paste0(regressions, m1,'_1 + ',m1,'_2 ~ ',measures[[m1]],'\n');
  if (measures[[m2]] != '') regressions <- paste0(regressions, m2,'_1 + ',m2,'_2 ~ ',measures[[m2]],'\n');
  if (measures[[m3]] != '') regressions <- paste0(regressions, m3,'_1 + ',m3,'_2 ~ ',measures[[m3]],'\n');
  if (model[['C']]==0) {
    mod <- paste0(mod, '\nc_{L1} == 0\nrC_{L1}_{L2} == 0\nrC_{L1}_{L3} == 0\nrC_{L2}_{L3} == 0\n')
  }
  # Insert regressions
  mod <- gsub(x=mod,pattern = '\\{REGRESSIONS\\}', replacement=regressions)
  # Insert measurements
  mod <- gsub(x=mod, pattern = '\\{M1\\}',replacement = m1)
  mod <- gsub(x=mod, pattern = '\\{M2\\}',replacement = m2)
  mod <- gsub(x=mod, pattern = '\\{M3\\}',replacement = m3)
  # Insert labels
  mod <- gsub(x=mod, pattern = '\\{L1\\}',replacement = labels[1])
  mod <- gsub(x=mod, pattern = '\\{L2\\}',replacement = labels[2])
  mod <- gsub(x=mod, pattern = '\\{L3\\}',replacement = labels[3])
  return(mod)
}