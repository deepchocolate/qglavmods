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
    print(parameters[[measure]])
    params <- paste0(parameters[[measure]], collapse=')^2 + (')
    m <- paste0(m, 'var_', measure, ' := (',params, ")^2\n")
  }
  return(m)
}