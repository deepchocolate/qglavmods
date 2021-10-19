# Bivariate twin model
# Univariate twin model
#' Create a twin model
#' @export
bivTwinModel <- function (mod1, mod2) {
  return(new('Bivariate', mod1, mod2))
}

validate <- function (object) {
  
}

setClass('Bivariate', contains='Univariate', 
         representation = list(
           mod1 = "Univariate",
           mod2 = "Univariate",
           corrParamPrefix = 'list'
           ))
# Constructor
setMethod('initialize', 'Bivariate',
          function (.Object, mod1, mod2) {
            .Object@mod1 = mod1
            .Object@mod2 = mod2
            .Object@corrParamPrefix = list(A='rA', C='rC', E='rE')
            callNextMethod(.Object)
          })
"### Cross-twin cross-trait correlations
A_{L1}_1 ~~ c(rA_{L1}_{L2}, rA_{L1}_{L2}_2)*A_{L2}_2
A_{L2}_1 ~~ c(rA_{L1}_{L2}, rA_{L1}_{L2}_2)*A_{L1}_2
# Constraints
rA_{L1}_{L2}_2 == 0.5*rA_{L1}_{L2}
### Cross-trait correlations
A_{L1}_1 ~~ c(rA_{L1}_{L2}, rA_{L1}_{L2})*A_{L2}_1
A_{L1}_2 ~~ c(rA_{L1}_{L2}, rA_{L1}_{L2})*A_{L2}_2"
setGeneric('getLatentCorrelations', function (object, objA, objB) standardGeneric('getLatentCorrelations'))
setMethod('getLatentCorrelations', signature('Bivariate', 'Univariate', 'Univariate'), 
          function (object, objA, objB) {
            facs1 <- getLatentFactors(objA)
            facs2 <- getLatentFactors(objB)
            # Genetics
            facA1 <- suffixedLatent(objA, 'A')
            facA2 <- suffixedLatent(objB, 'A')
            # Cross-twin cross trait
            meas1 <- getMeasure(objA)
            meas2 <- getMeasure(objB)
            paramLabel <- getCorrParamLabel(object, 'A', objA, objB)
            parm <- paste0('c(', paramLabel, ', ', paramLabel,'_2)*')
            parm2 <- paste0('c(', paramLabel, ', ', paramLabel,')*')
            out <- paste0(facA1[1],' ~~ ', parm, facA2[2], '\n')
            out <- paste0(out, facA1[2],' ~~ ',parm, facA2[1], '\n')
            #out <- paste0(out, 'rA_m1_m2_2 == 0.5*rA_m1_m2\n')
            out <- paste0(out, paramLabel, '_2 == 0.5*', paramLabel, '\n')
            # Cross-trait
            out <- paste0(out, facA1[1],' ~~ ', parm2, facA2[1], '\n')
            out <- paste0(out, facA1[2],' ~~ ', parm2, facA2[2], '\n')
            
            # Shared environment
            # Cross-twin cross-trait correlations
            facC1 <- suffixedLatent(objA, 'C')
            facC2 <- suffixedLatent(objB, 'C')
            paramLabel <- getCorrParamLabel(object, 'C', objA, objB)
            parm <- paste0('c(', paramLabel, ', ', paramLabel,')*')
            # Cross twin cross trait
            out <- paste0(out, facC1[1], ' ~~ ', parm, facC2[2], '\n')
            out <- paste0(out, facC1[2], ' ~~ ', parm, facC2[1], '\n')
            # Cross trait
            out <- paste0(out, facC1[1], ' ~~ ', parm, facC2[1], '\n')
            out <- paste0(out, facC1[2], ' ~~ ', parm, facC2[2], '\n')
            # Unique environment
            facE1 <- suffixedLatent(objA, 'E')
            facE2 <- suffixedLatent(objB, 'E')
            paramLabel <- getCorrParamLabel(object, 'E', objA, objB)
            parm <- paste0('c(', paramLabel, ', ', paramLabel, ')*')
            # Cross trait
            meas1 <- suffixedMeasures(objA)
            meas2 <- suffixedMeasures(objB)
            out <- paste0(out, meas1[1], ' ~~ ', parm, meas2[2], '\n')
            out <- paste0(out, meas1[2], ' ~~ ', parm, meas2[1], '\n')
            out
          })

setGeneric('getCorrParamPrefix', function (object, factor) standardGeneric('getCorrParamPrefix'))
setMethod('getCorrParamPrefix', signature('Bivariate', 'character'),
          function (object, factor) {
            object@corrParamPrefix[[factor]]
          })
setGeneric('getCorrParamLabel', function (object, factor, ...) standardGeneric('getCorrParamLabel'))
setMethod('getCorrParamLabel', signature('Bivariate', 'character'),
          function (object, factor, ...) {
            objs <- list(...)
            nms <- c()
            for (obj in objs) {
              nms <- c(nms, getMeasure(obj))
            }
            nms <- paste(nms, collapse='_')
            prm <- paste0(getCorrParamPrefix(object, factor), '_', nms)
            prm
          })

setMethod('objectToChar', signature('Bivariate'),
          function (object) {
            m1meas <- getMeasure(object@mod1)
            m2meas <- getMeasure(object@mod2)
            m1 <- setParameterLabels(object@mod1, A=paste0('a_',m1meas), C=paste0('c_', m1meas), E=paste0('e_', m1meas))
            suf <- paste0(m1meas, object@latentSuffix)
            m1 <- setLatentSuffix(m1, suf)
            m2 <- setParameterLabels(object@mod2, A=paste0('a_',m2meas), C=paste0('c_', m2meas), E=paste0('e_', m2meas))
            suf <- paste0(m2meas, object@latentSuffix)
            m2 <- setLatentSuffix(m2, suf)
            crs <- getLatentCorrelations(object, m1, m2)
            m1 <- objectToChar(m1)
            m2 <- objectToChar(m2)
            out <- paste0(m1, '\n', m2, '\n', crs)
            out
          })