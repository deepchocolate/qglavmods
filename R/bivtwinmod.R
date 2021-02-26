#' Generate syntax for a bivariate twin model
#'
#' @param measT1 A list of variable name and a designated label eg label(varName1="label1") for measure 1.
#' @param measT2 A list of variable name and a designated label eg label(varName1="label2") for measure 2.
#' @param regressions A list of lists for the regressions for measT1 and measT1.
#' @param model A list of variance components (A,C,D,E) and labels, eg list(A='a',C='c',E='e').
#' @export
bivTwinMod <- function (measT1, measT2,regressions=list(list(int=1,fem='female'),list(int=1,fem='female')), model=list(A='a',C='c',E='e')) {
  corrs <- list(A=c(0.5,1),C=c(1,1),D=c(0.25,1))
  defs <- c()
  defsVar <- list()
  defsCorr <- list()
  m <- '# Traits in twins are assumed to be uncorrelated\n'
  for (i in 1:2) m <- paste0(m,names(measT1[i]),' ~~ 0*',names(measT2[i]),'\n');
  m <- paste0(m,'# Regressions\n')
  for (i in 1:2) {
    reg <- c()
    for (x in names(regressions[[i]])) {
      reg <- c(reg,paste0('c(',x,',',x,')*',regressions[[i]][[x]],ifelse(!is.numeric(regressions[[i]][[x]]),'_1','')))
    }
    reg <- paste(reg,collapse='+')
    m <- paste0(m,names(measT1[i]),' ~ ',paste(reg,collapse=' + '),'\n')
    reg <- c()
    for (x in names(regressions[[i]])) reg <- c(reg,paste0('c(',x,',',x,')*',regressions[[i]][[x]],ifelse(!is.numeric(regressions[[i]][[x]]),'_2','')))
    reg <- paste(reg,collapse='+')
    m <- paste0(m,names(measT2[i]),' ~ ',paste(reg,collapse=' + '),'\n')
  }
  for (fac in names(model)) {
    defsCorr[[fac]] <- c()
    if (fac == 'A') {m <- paste0(m,'# Genetics\n');}
    else if (fac== 'C') {m <- paste0(m,'# Shared environment\n');}
    else if (fac == 'E') {m <- paste0(m,'# Unique environment\n');}
    m <- paste0(m,'# Twin 1 measures\n')
    for (out in names(measT1)) {
      paraLab <- paste0(model[fac],measT1[out])
      m <- paste0(m,fac,out,' =~ c(',paraLab,',',paraLab,')*',out,'\n')
    }
    corr <- paste(corrs[fac],collapse=',')
    if (fac == 'E') corr <- 0;
    m <- paste0(m, fac,names(measT1[1]),' ~~ ',corr,'*',fac,names(measT2[1]),'\n')
    m <- paste0(m, '# Twin 2 measures\n')
    for (out in names(measT2)) {
      paraLab <- paste0(model[fac],measT2[out])
      m <- paste0(m,fac,out,' =~ c(',paraLab,',',paraLab,')*',out,'\n')
    }
    corr <- paste(corrs[fac],collapse=',')
    if (fac == 'E') corr <- 0;
    m <- paste0(m, fac,names(measT1[2]),' ~~ ',corr,'*',fac,names(measT2[2]),'\n');
    m <- paste0(m,'# Cross-twin cross-trait correlations\n')
    corrParam <- paste0('*r',fac)
    corrParam <- paste0('c(',paste0(corrs[[fac]],c(corrParam,corrParam),collapse=','),')')
    if (fac == 'E') corrParam <- 0;
    if (fac == 'C') corrParam <- paste0('c(r',fac,',r',fac,')');
    #for (i in 1:2) {
    m <- paste0(m,fac,names(measT1[1]),' ~~ ',corrParam,'*',fac,names(measT2[2]),'\n')
    m <- paste0(m,fac,names(measT1[2]),' ~~ ',corrParam,'*',fac,names(measT2[1]),'\n')
    #}
    m <- paste0(m,'# Cross-trait correlations\n')
    corrParam <- paste0('c(','r',fac,',r',fac,')')
    #if (fac == 'E') corrParam <- 0;
    m <- paste0(m,fac,names(measT1[1]),' ~~ ',corrParam,'*',fac,names(measT1[2]),'\n')
    m <- paste0(m,fac,names(measT2[1]),' ~~ ',corrParam,'*',fac,names(measT2[2]),'\n')
    for (i in 1:2) {
      if (is.null(defsVar[[measT1[[i]]]])) defsVar[measT1[[i]]] <- c()
      corrParam <- paste0(model[fac],measT1[i])
      sqParam <- paste0(corrParam,'SQ')
      defs <- c(defs, paste0(corrParam,'SQ := ',corrParam,'*',corrParam))
      defsVar[[measT1[[i]]]] <- c(defsVar[[measT1[[i]]]],sqParam)
      defsCorr[[fac]] <- c(defsCorr[[fac]],corrParam)
    }
  }
  m <- paste0(m,'# Defined parameters\n# Squared correlations\n')
  m <- paste0(m,paste0(defs,collapse='\n'),'\n')
  m <- paste0(m,'# Variance explained in each trait\n')
  #print(defsVar)
  for (x in names(defsVar)) {
    tot <- paste0(defsVar[[x]],collapse='+')
    for (y in defsVar[[x]]) {
      sh <- paste0(y,'_share := ')
      m <- paste0(m,sh,y,'/(', tot,')\n')
    }
    #m <- paste0(m,sh,'\n')
  }
  m <- paste0(m, '# Genetic and environmental correlations between traits\n')
  # Add these definitions so that all parameters are easily accessible through the lhs column in results
  m <- paste0(m,'corr_A := rA\ncorr_C := rC\ncorr_E := rE\n')
  m <- paste0(m,'# Genetic and environmental contributions to phenotypic correlations\n')
  corrLabs <- c()
  for (x in names(defsCorr)) {
    corrLab <- paste0('contr_',x)
    corrLabs <- c(corrLabs,corrLab)
    m <- paste0(m,corrLab,' := r',x,'*',paste0(defsCorr[[x]],collapse='*'),'\n')
  }
  tot <- paste0(corrLabs,collapse='+')
  for (x in corrLabs) {
    m <- paste0(m,x,'_share := ',x,'/(',tot,')','\n')
  }
  return(m)
}

