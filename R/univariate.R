# Univariate twin model
#' Create a twin model
#' @export
univTwinModel <- function () {
  return(new('Univariate'))
}

#' The univariate class
#' 
#' @slot factors A, C, D, E factors in a twin model.
#' @slot regressions Regressions to perform on measurements
#' @slot suffix Suffix for the measurement in twin and co twin
#' @slot latentSuffix Suffix for latent variables.
#' @export
setClass('Univariate', contains='twinModel', 
         slots=list(
           measure = 'character',
           factors = 'list',
           regressions = 'formula',
           suffix = 'vector',
           latentSuffix='vector'
           ))
setMethod('initialize', 'twinModel',
          function (.Object) {
            .Object@factors = list(A='a', C='c',E='e')
            .Object@regressions = formula()
            .Object@suffix = c('_1', '_2')
            .Object@latentSuffix = c('_1', '_2')
            callNextMethod(.Object)
          })

setGeneric('getMeasure', function (object) standardGeneric('getMeasure'))
setMethod('getMeasure', signature('Univariate'),
          function (object) {
            object@measure
          })

### Public methods
#' Get latent factors (ACDE)
#' 
#' @param object An object of class twinModel
#' @export
setGeneric('getLatentFactors', function (object) standardGeneric('getLatentFactors'))
setMethod('getLatentFactors', signature('Univariate'),
          function (object) {
            names(object@factors)
          })
#' Set latent factors (ACDE)
#' 
#' @param ... Each factor specified as factor=label, eg A='a'
#' @export
setGeneric('setLatentFactors', function (object, ...) standardGeneric('setLatentFactors'))
setMethod('setLatentFactors', signature('Univariate'),
          function (object, ...) {
            object@factors <- list(...)
            object
          })

setGeneric('setLatentSuffix', function (object, suffix) standardGeneric('setLatentSuffix'))
setMethod('setLatentSuffix', signature('Univariate', 'vector'),
          function (object, suffix) {
            object@latentSuffix <- suffix
            object
          })

setGeneric('suffixedLatent', function (object, factor) standardGeneric('suffixedLatent'))
setMethod('suffixedLatent', signature('Univariate', 'character'),
          function (object, factor) {
            paste0(factor, object@latentSuffix)
          })
#' Get parameter labels
#' 
#' @param object An object of class twinModel
setGeneric('getParameterLabels', function (object) standardGeneric('getParameterLabels'))
setMethod('getParameterLabels', signature('Univariate'),
          function (object) {
            unlist(object@factors, use.names = F)
          })
setGeneric('getLatentParameterLabel', function (object, fac) standardGeneric('getLatentParameterLabel'))
setMethod('getLatentParameterLabel', signature('Univariate', 'character'),
          function (object, fac) {
            if (!fac %in% getLatentFactors(object)) stop(paste0('Unknonw latent factor ', fac))
            object@factors[[fac]]
          })
#' Set parameter labels
#' 
#' 
setGeneric('setParameterLabels', function (object, ...) standardGeneric('setParameterLabels'))
setMethod('setParameterLabels', signature('Univariate'),
          function (object, ...) {
            lst <- list(...)
            for (el in names(lst)) {
              if (!el %in% names(object@factors)) stop(paste0(el, ' is not among the specified latent factors!'))
              object@factors[[el]] <- lst[[el]]
            }
            object
          })

#' Get loading on a latent variable
setGeneric('getLoading', function (object, fac) standardGeneric('getLoading'))
setMethod('getLoading', signature('Univariate', 'character'),
          function (object, fac) {
            ms <- suffixedMeasures(object)
            latParm <- getLatentParameterLabel(object, fac)
            labs <- paste0('c(', latParm, ', ', latParm,')*')
            lats <- suffixedLatent(object, fac)
            paste0(lats, rep(' =~ ', 2), rep(labs, 2), ms)
          })

setGeneric('getResiduals', function (object) standardGeneric('getResiduals'))
setMethod('getResiduals', signature('Univariate'),
          function (object) {
            ms <- suffixedMeasures(object)
            latParm <- getLatentParameterLabel(object, 'E')
            labs <- paste0(' ~~ c(', latParm,', ',latParm,')*')
            paste0(ms, labs, ms)
          })

setMethod('latentFactors', signature('Univariate'),
          function (object) {
            facs <- getLatentFactors(object)
            op <- ''
            for (fac in facs) {
              if (fac != 'E') {
                lats <- suffixedLatent(object, fac)
                ld <- paste(getLoading(object, fac), collapse='\n')
                op <- paste(op, ld, getCovariance(object, fac, lats), collapse='\n', sep='\n')
              }
            }
            if ('E' %in% facs) {
              res <- paste0(getResiduals(object), collapse='\n')
              op <- paste0(op, '\n', res)
            }
            op
          })


setMethod('getDefinitions', signature('Univariate'),
          function (object) {
            defs <- ''
            facs <- getLatentFactors(object)
            f2Lab <- object@factors
            totVar <- paste0('(', paste(unlist(f2Lab, use.names=F), collapse = '^2 + '), ')')
            defs <- paste0(defs, 'V_', getMeasure(object), ' := ', totVar, '\n')
            m <- getMeasure(object)
            for (fac in c('A', 'C', 'D')) {
              if (fac %in% facs) defs <- paste0(defs, fac, '_', m, '_share := ', getLatentParameterLabel(object, fac), '^2/V_', m, '\n')
            }
            if ('E' %in% facs) {
              defs <- paste0(defs, 'E_', m,'_share := ', getLatentParameterLabel(object, 'E'),'/V_', m)
            }
            defs
          })

setMethod('getCovariance', signature('Univariate', 'character'),
          function (object, factor, latents) {
            if (!factor %in% getLatentFactors(object)) stop('Factor must be one of ACDE');
            if (factor == 'A') cv <- paste0(latents[1], ' ~~ c(1, 0.5)*', latents[2]);
            if (factor == 'C') cv <- paste0(latents[1], ' ~~ c(1, 1)*', latents[2]);
            cv
          })

setGeneric('getRegressions', function (object) standardGeneric('getRegressions'))
setMethod('getRegressions', signature('Univariate'),
          function (object) {
            regs <- suffixFormula(object)
            paste(regs, collapse="\n")
          })

setMethod('objectToChar', signature('Univariate'),
          function (object) {
            regs <- ''
            if (length(object@regressions) > 0) {
              regs <- getRegressions(object)
            }
            latents <- paste0(latentFactors(object), '\n')
            defs <- getDefinitions(object)
            paste0(regs, latents, defs, sep='\n')
          })

#' Add regression to the twin model
#' @param reg A regression formula for a measurement without suffix for measurements in twin and co-twin.
#' @export
setGeneric('regress', function (object, reg) standardGeneric('regress'))
setMethod('regress', signature('twinModel', 'formula'),
          function (object, reg) {
            object@regressions <- reg
            object@measure <- as.character(reg[[2]])
            object
          })
setGeneric('removeRegressions', function (object) standardGeneric('removeRegressions'))
setMethod('removeRegressions', signature('Univariate'),
          function (object) {
            object@regressions <- formula()
            object
          })
setGeneric('suffixFormula', function (object) standardGeneric('suffixFormula'))
setMethod('suffixFormula', signature('Univariate'),
          function (object) {
            ms <- suffixedMeasures(object)
            regs <- object@regressions
            op <- c()
            tmp <- as.character(regs)
            for (m in ms) {
              op <- c(op, paste0(m, ' ~ ', tmp[[3]]))
            }
            op
          })