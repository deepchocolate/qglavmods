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
           corrParamPrefix = 'list',
           constraints = 'list'
           ))
# Constructor
setMethod('initialize', 'Bivariate',
          function (.Object, mod1, mod2) {
            .Object@corrParamPrefix = list(A='rA', C='rC', E='rE')
            .Object@constraints = list()
            m1meas <- getMeasure(mod1)
            m2meas <- getMeasure(mod2)
            mod1 <- setParameterLabels(mod1, A=paste0('a_',m1meas), C=paste0('c_', m1meas), E=paste0('e_', m1meas))
            suf <- paste0(m1meas, mod1@latentSuffix)
            mod1 <- setLatentSuffix(mod1, suf)
            mod2 <- setParameterLabels(mod2, A=paste0('a_',m2meas), C=paste0('c_', m2meas), E=paste0('e_', m2meas))
            suf <- paste0(m2meas, mod2@latentSuffix)
            mod2 <- setLatentSuffix(mod2, suf)
            .Object@mod1 = mod1
            .Object@mod2 = mod2
            callNextMethod(.Object)
          })

#' Add a constraint to the model
#' @export
setGeneric('addConstraint', function (object, param, rhs) standardGeneric('addConstraint'))
setMethod('addConstraint', signature('Bivariate', 'character', 'character'),
          function (object, param, rhs) {
            object@constraints[[param]] = rhs
            object
          })

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
            out <- paste0(out, meas1[1], ' ~~ ', parm, meas2[1], '\n')
            out <- paste0(out, meas1[2], ' ~~ ', parm, meas2[2], '\n')
            out
          })

setMethod('getDefinitions', signature('Bivariate'),
          function (object) {
            # Variance explained in each trait
            def <- ''
            rA <- getCorrParamLabel(object, 'A', object@mod1, object@mod2)
            rC <- getCorrParamLabel(object, 'C', object@mod1, object@mod2)
            rE <- getCorrParamLabel(object, 'E', object@mod1, object@mod2)
            # ACE correlations
            # With the exception of rE these are just defined to reduce the no of estimates
            def <- paste0(def, 'corr_A := ', rA, '\n')
            def <- paste0(def, 'corr_C := ', rC, '\n')
            # Within this framework rE is covariance and has to be converted to a correlation
            def <- paste0(def, 'corr_E := ', rE, '/(sqrt(', object@mod1@factors[['E']], ')*sqrt(', object@mod2@factors[['E']],'))\n')
            # Contributions to phenotypic correlations
            m1 <- getMeasure(object@mod1)
            m2 <- getMeasure(object@mod2)
            def <- paste0(def, 'contr_A := corr_A*sqrt(A_', m1, '_share', ')*sqrt(A_', m2, '_share)\n')
            def <- paste0(def, 'contr_C := corr_C*sqrt(A_', m1, '_share', ')*sqrt(A_', m2, '_share)\n')
            def <- paste0(def, 'contr_E := corr_E*sqrt(A_', m1, '_share', ')*sqrt(A_', m2, '_share)\n')
            def <- paste0(def, 'pheno := contr_A + contr_C + contr_E\n')
            def <- paste0(def, 'contr_A_share := contr_A/pheno\n')
            def <- paste0(def, 'contr_C_share := contr_C/pheno\n')
            def <- paste0(def, 'contr_E_share := contr_E/pheno\n')
            def
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
            #m1meas <- getMeasure(object@mod1)
            #m2meas <- getMeasure(object@mod2)
            #m1 <- setParameterLabels(object@mod1, A=paste0('a_',m1meas), C=paste0('c_', m1meas), E=paste0('e_', m1meas))
            #suf <- paste0(m1meas, object@latentSuffix)
            #m1 <- setLatentSuffix(m1, suf)
            #m2 <- setParameterLabels(object@mod2, A=paste0('a_',m2meas), C=paste0('c_', m2meas), E=paste0('e_', m2meas))
            #suf <- paste0(m2meas, object@latentSuffix)
            #m2 <- setLatentSuffix(m2, suf)
            crs <- getLatentCorrelations(object, object@mod1, object@mod2)
            r1 <- getRegression(object@mod1)
            r2 <- getRegression(object@mod2)
            lat1 <- latentFactors(object@mod1)
            lat2 <- latentFactors(object@mod2)
            defs <- getDefinitions(object)
            def1 <- getDefinitions(object@mod1)
            def2 <- getDefinitions(object@mod2)
            cnstrnts <- ''
            for (cnstr in names(object@constraints)) {
              cnstrnts <- paste0(cnstrnts, cnstr, object@constraints[[cnstr]], '\n')
            }
            out <- paste0(r1, r2, paste0(lat1,'\n'), paste0(lat2,'\n'), crs, def1, def2, defs, cnstrnts, collapse='\n\n')
            out
          })