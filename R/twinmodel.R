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
# Get the measure with suffixes
setGeneric('suffixedMeasures', function (object) standardGeneric('suffixedMeasures'))
setMethod('suffixedMeasures', signature('twinModel'),
          function (object) {
            paste0(object@measure, object@suffix)
          })

# Get the covariance structure
setGeneric('getCovariance', function (object, factor, latents) standardGeneric('getCovariance'))

### Public methods
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

setMethod('as.character', signature('twinModel'),
          function (x) {
            objectToChar(x);
          })

### Methods that need to be defined among inheriting classes
# Method for ACDE factors
setGeneric('latentFactors', function (object) standardGeneric('latentFactors'))
# Method to convert object to character lavaan syntax
setGeneric('objectToChar', function (object) standardGeneric('objectToChar'))
# Method to get parameter definitions
setGeneric('getDefinitions', function (object) standardGeneric('getDefinitions'))
