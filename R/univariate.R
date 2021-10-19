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
            object@regressions[[2]]
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
setGeneric('setLatentSuffix', function (object, suffix) standardGeneric('setLatentSuffix'))
setMethod('setLatentSuffix', signature('Univariate', 'vector'),
          function (object, suffix) {
            object@latentSuffix <- suffix
            object
          })

setGeneric('suffixedLatent', function (object, factor) standardGeneric('suffixedLatent'))
setMethod('suffixedLatent', signature('twinModel', 'character'),
          function (object, factor) {
            paste0(factor, object@latentSuffix)
          })

setMethod('latentFactors', signature('Univariate'),
          function (object) {
            ms <- suffixedMeasures(object)
            facs <- getLatentFactors(object)
            f2Lab <- object@factors
            defs <- ''
            totVar <- paste0('(', paste(unlist(f2Lab, use.names=F), collapse = '^2 + '), ')')
            if ('A' %in% facs) {
              labs <- paste0('c(', f2Lab[['A']],', ',f2Lab[['A']],')*')
              lats <- suffixedLatent(object, 'A')
              A <- paste(paste0(lats, rep(' =~ ', 2), rep(labs, 2), ms), collapse='\n')
              A <- paste0(A, '\n', getCovariance(object, 'A', lats), collapse='\n')
              defs <- paste0(defs, 'a2_share := ', f2Lab[['A']],'^2/', totVar, '\n')
            }
            if ('C' %in% facs) {
              labs <- paste0('c(', f2Lab[['C']],', ',f2Lab[['C']],')*')
              lats <- suffixedLatent(object, 'C')
              C <- paste(paste0(lats, rep(' =~ ', 2), rep(labs, 2), ms), collapse='\n')
              C <- paste(C, getCovariance(object, 'C', lats), collapse='\n', sep='\n')
              defs <- paste0(defs, 'c2_share := ', f2Lab[['C']],'^2/', totVar, '\n')
            }
            if ('E' %in% facs) {
              labs <- paste0(' ~~ c(', f2Lab[['E']],', ',f2Lab[['E']],')*')
              E <- paste(paste0(ms, rep(labs, 2), ms), collapse='\n')
              defs <- paste0(defs, 'e2_share := ', f2Lab[['E']],'/', totVar, '\n')
            }
            paste(A, C, E, defs, collapse='\n', sep='\n')
          })

setMethod('objectToChar', signature('Univariate'),
          function (object) {
            regs <- ''
            if (length(object@regressions) > 0) {
              regs <- suffixFormula(object)
              regs <- paste(regs, collapse="\n")
            }
            latents <- latentFactors(object)
            out <- paste0(regs, '\n', latents)
            out
          })