#' Create a twin model
#' @export
twinModel <- function () {
  return(new('twinModel'))
}

#' The twinModel class
#' 
#' @slot factors A, C, D, E factors in a twin model.
#' @slot regressions Regressions to perform on measurements
#' @slot suffix Suffix for the measurement in twin and co twin
#' @export
setClass('twinModel', prototype=list())

### Internal methods
setGeneric('suffixedMeasures', function (object) standardGeneric('suffixedMeasures'))
setMethod('suffixedMeasures', signature('twinModel'),
          function (object) {
            paste0(object@regressions[[2]], object@suffix)
          })
setGeneric('suffixFormula', function (object) standardGeneric('suffixFormula'))
setMethod('suffixFormula', signature('twinModel'),
          function (object) {
            ms <- suffixedMeasures(object)
            regs <- object@regressions
            r1 <- update(regs, paste0(ms[1], '~.'))
            r2 <- update(regs, paste0(ms[2], '~.'))
            c(r1, r2)
          })

# Get the covariance structure
setGeneric('getCovariance', function (object, factor, latents) standardGeneric('getCovariance'))
setMethod('getCovariance', signature('twinModel', 'character'),
          function (object, factor, latents) {
            if (!factor %in% getLatentFactors(object)) stop('Factor must be one of ACDE');
            if (factor == 'A') cv <- paste0(latents[1], ' ~~ c(1, 0.5)*', latents[2]);
            if (factor == 'C') cv <- paste0(latents[1], ' ~~ c(1, 1)*', latents[2]);
            cv
          })

### Methods that need to be defined among inheriting classes
# Method for ACDE factors
setGeneric('latentFactors', function (object) standardGeneric('latentFactors'))
# Method to convert object to character lavaan syntax
setGeneric('objectToChar', function (object) standardGeneric('objectToChar'))


### Public methods
#' Get latent factors (ACDE)
#' 
#' @param object An object of class twinModel
setGeneric('getLatentFactors', function (object) standardGeneric('getLatentFactors'))
setMethod('getLatentFactors', signature('twinModel'),
          function (object) {
            names(object@factors)
          })
#' Set latent factors (ACDE)
#' 
#' 
setGeneric('setLatentFactors', function (object, ...) standardGeneric('setLatentFactors'))
setMethod('setLatentFactors', signature('twinModel'),
          function (object, ...) {
            object@factors <- list(...)
            object
          })

#' Get parameter labels
#' 
#' @param object An object of class twinModel
setGeneric('getParameterLabels', function (object) standardGeneric('getParameterLabels'))
setMethod('getParameterLabels', signature('twinModel'),
          function (object) {
            unlist(object@factors, use.names = F)
          })


#' Add regression to the twin model
#' @param reg A regression formula for a measurement without suffix for measurements in twin and co-twin.
#' @export
setGeneric('regress', function (object, reg) standardGeneric('regress'))
setMethod('regress', signature('twinModel', 'formula'),
          function (object, reg) {
            object@regressions <- reg
            return(object)
          })
#' Set the suffix for measurments in twin and co-twin
#' @param t1 Measurement suffix for twin 1
#' @param t2 Measurement suffix for twin 2
#' @export
setGeneric('setSuffix', function (object, t1, t2) standardGeneric('setSuffix'))
setMethod('setSuffix', signature('twinModel'),
          function (object, t1, t2) {
            object@suffix <- c(t1, t2)
            return(object)
          })

#setMethod('show', signature('twinModel'),
#          function (object) {
#            op <- objectToChar(object)
#            cat(op)
#          })
setMethod('as.character', signature('twinModel'),
          function (x) {
            objectToChar(x);
          })
