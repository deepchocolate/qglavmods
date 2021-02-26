#' Generate syntax for a univariate twin model
#'
#' @param measT1 Variable name of measurement in twin one.
#' @param measT2 Variable name of measurement in twin two.
#' @param model A list of variance components (A,C,D,E) and labels, eg list(A='a',C='c',E='e').
#' @param append Any additional syntax.
#' @param varLabels Currently not used.
#' @export
univTwinMod <- function (measT1, measT2, model=list(A='a',C='c',E='e'), append='',varLabels=list(A='add',C='com',D='dom',E='uni'),
                            independentT1=list(int='1'),independentT2=list(int='1')) {
  corrs <- list(A=c(0.5,1),C=c(1,1),D=c(0.25,1))
  o <- '# Measurments are uncorrelated\n'
  o <- paste0(measT1,'~~0*',measT2,'\n')
  if (length(independentT1) > 0 | length(independentT2) > 0) {
    o <- paste0(o,'# Regressions\n')
    r <- c()
    for (lab in names(independentT1)) {
      r <- c(r,paste0('c(',lab,',',lab,')*',independentT1[lab]))
    }
    o <- paste0(o,measT1,' ~ ',paste(r,collapse=' + '),'\n')
    r <- c()
    for (lab in names(independentT2)) {
      r <- c(r,paste0('c(',lab,',',lab,')*',independentT2[lab]))
    }
    o <- paste0(o,measT2,' ~ ',paste(r,collapse=' + '),'\n')
  }
  totalVars <- c()
  for (fac in names(model)) {
    fac1 <- paste0(fac,'1')
    fac2 <- paste0(fac,'2')
    facName <- model[fac]
    totalVars <- c(totalVars,paste0(facName,'2'))
    o <- paste0(o,fac1,' =~ c(',facName,',',facName,')*',measT1,'\n')
    o <- paste0(o,fac2,' =~ c(',facName,',',facName,')*',measT2,'\n')
    if (fac=='A') {
      o <- paste0(o,fac1,' ~~ c(0.5,1)*',fac2,'\n')
    }
    if (fac=='C') {
      o <- paste0(o,fac1,'~~1*',fac2,'\n')
    }
    if (fac == 'D') {
      o <- paste0(o,fac1,' ~~ c(0.25,1)*',fac2,'\n')
    }
    if (fac=='E') {
      # Residual variation is uncorrelated
      o <- paste0(o,fac1,'~~0*',fac2,'\n')
    }
    o <- paste0(o,facName,'2 := ',facName,'*',facName,'\n')
  }
  totVar <- paste0('(',paste0(totalVars,collapse='+'),')')
  geneFacs <- c()
  for (fac in names(model)) {
    if (fac %in% c('A','D')) geneFacs <- c(geneFacs,paste0(model[fac],'2'))
    o <- paste0(o,varLabels[fac],' := ',model[fac],'2/',totVar,'\n')
  }
  o <- paste0(o,'h2 := (',paste0(geneFacs,collapse='+'),')/',totVar,'\n')
  return(paste0(o,append))
}

