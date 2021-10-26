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

#' Get loadings.
#' This is used for the cholesky model.
setMethod('getLoading', signature('Bivariate', 'character'),
          function (object, fac) {
            l1 <- paste0(getLoading(object@mod1, fac), collapse='\n')
            l2 <- getLoading(object@mod2, fac)
            lbl <- getCorrParamLabel(object, fac, object@mod1, object@mod2)
            if (fac == 'A') parm <- paste0('c(', lbl, ', ', lbl, ')*')
            else parm <- paste0('c(', lbl, ', ', lbl, '_2)*')
            m1 <- suffixedMeasures(object@mod1)
            l2 <- paste0(l2, ' + ', parm, m1)
            lats1 <- suffixedLatent(object@mod1, fac)
            lats2 <- suffixedLatent(object@mod2, fac)
            cv1 <- getCovariance(object@mod1, fac, lats1)
            cv2 <- getCovariance(object@mod2, fac, lats2)
            out <- paste(l1, paste0(l2, collapse='\n'), cv1, cv2, sep='\n')
            if (fac == 'A') out <- paste0('\n', lbl, '_2 == 0.5*', lbl, '\n')
            out
          })


setMethod('getResiduals', signature('Bivariate'),
          function (object) {
            l1 <- getResiduals(object@mod1)
            m1 <- suffixedMeasures(object@mod1)
            parm <- getCorrParamLabel(object, 'E', object@mod1, object@mod2)
            l2 <- getResiduals(object@mod2)
            paste0(l1, '\n', paste0(l2, ' + c(', parm, ', ', parm, ')*'), m1)
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

#' Get parameter definitions for a bivariate cholesky model.
#' Contrary to the correlated factors solution, the definitions in the univariate
#' models cannot be used here.
setGeneric('getCholeskyDefinitions', function (object) standardGeneric('getCholeskyDefinitions'))
setMethod('getCholeskyDefinitions', signature('Bivariate'),
          function (object) {
            defs <- ''
            lt1 <- getParameterLabels(object@mod1)
            m1 <- getMeasure(object@mod1)
            V1A <- paste0('V_A_', m1, ' := ', getLatentParameterLabel(object@mod1, 'A'), '^2')
            V1C <- paste0('V_C_', m1, ' := ', getLatentParameterLabel(object@mod1, 'C'), '^2')
            V1E <- paste0('V_E_', m1, ' := ', getLatentParameterLabel(object@mod1, 'E'))
            defs <- paste0(defs, V1A, '\n', V1C, '\n', V1E, '\n')
            defs <- paste0(defs, '\nV_', m1, ' := V_A_', m1, ' + V_C_', m1, ' + V_E_', m1, '\n')
            defs <- paste0(defs, 'A_', m1, '_share := V_A_', m1, '/V_', m1, '\n')
            defs <- paste0(defs, 'C_', m1, '_share := V_C_', m1, '/V_', m1, '\n')
            defs <- paste0(defs, 'E_', m1, '_share := V_E_', m1, '/V_', m1, '\n')
            m2 <- getMeasure(object@mod2)
            V2A <- paste0('V_A_', m2, ' := ', getLatentParameterLabel(object@mod2, 'A'), '^2 + ', getCorrParamLabel(object, 'A', object@mod1, object@mod2),'^2')
            V2C <- paste0('V_C_', m2, ' := ', getLatentParameterLabel(object@mod2, 'C'), '^2 + ', getCorrParamLabel(object, 'C', object@mod1, object@mod2),'^2')
            V2E <- paste0('V_E_', m2, ' := ', getLatentParameterLabel(object@mod2, 'E'), ' + ', getCorrParamLabel(object, 'E', object@mod1, object@mod2))
            defs <- paste0(defs, V2A, '\n', V2C, '\n', V2E, '\n')
            defs <- paste0(defs, 'V_', m2, ' := V_A_', m2, ' + V_C_', m2, ' + V_E_', m2, '\n')
            defs <- paste0(defs, 'A_', m2, '_share := V_A_', m2, '/V_', m2, '\n')
            defs <- paste0(defs, 'C_', m2, '_share := V_C_', m2, '/V_', m2, '\n')
            defs <- paste0(defs, 'E_', m2, '_share := V_E_', m2, '/V_', m2, '\n')
            defs
          })

#' Get parameter definitions for a bivariate model.
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
            eLab1 <- getLatentParameterLabel(object@mod1, 'E')
            eLab2 <- getLatentParameterLabel(object@mod2, 'E')
            def <- paste0(def, 'corr_E := ', rE, '/(sqrt(', eLab1, ')*sqrt(', eLab2,'))\n')
            # Contributions to phenotypic correlations
            m1 <- getMeasure(object@mod1)
            m2 <- getMeasure(object@mod2)
            def <- paste0(def, 'contr_A := corr_A*sqrt(A_', m1, '_share', ')*sqrt(A_', m2, '_share)\n')
            def <- paste0(def, 'contr_C := corr_C*sqrt(C_', m1, '_share', ')*sqrt(C_', m2, '_share)\n')
            def <- paste0(def, 'contr_E := corr_E*sqrt(E_', m1, '_share', ')*sqrt(E_', m2, '_share)\n')
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

#' Formulate the correlated factors model
#' @param object The specified bivariate model.
#' @export
setGeneric("correlatedFactors", function (object) standardGeneric('correlatedFactors'))
setMethod('correlatedFactors', signature('Bivariate'),
          function (object) {
            crs <- getLatentCorrelations(object, object@mod1, object@mod2)
            r1 <- getRegressions(object@mod1)
            r2 <- getRegressions(object@mod2)
            lat1 <- latentFactors(object@mod1)
            lat2 <- latentFactors(object@mod2)
            defs <- getDefinitions(object)
            def1 <- getDefinitions(object@mod1)
            def2 <- getDefinitions(object@mod2)
            cnstrnts <- ''
            for (cnstr in names(object@constraints)) {
              cnstrnts <- paste0(cnstrnts, cnstr, object@constraints[[cnstr]], '\n')
            }
            out <- paste0(r1, '\n', r2, lat1, lat2, '\n', crs, def1, '\n', def2, '\n', defs, cnstrnts)
            out
          })

#' Formulate the cholesky model
#' @param object The specified bivariate model.
#' @export
setGeneric('cholesky', function (object) standardGeneric('cholesky'))
setMethod('cholesky', signature('Bivariate'),
          function (object) {
            r1 <- getRegressions(object@mod1)
            r2 <- getRegressions(object@mod2)
            l1 <- getLoading(object, 'A')
            l2 <- getLoading(object, 'C')
            res <- paste0(getResiduals(object), collapse='\n')
            cdefs <- getCholeskyDefinitions(object)
            defs <- getDefinitions(object)
            cnstrnts <- ''
            for (cnstr in names(object@constraints)) {
              cnstrnts <- paste0(cnstrnts, cnstr, object@constraints[[cnstr]], '\n')
            }
            paste(r1, r2, l1, l2, res, cdefs, defs, cnstrnts, sep='\n')
          })