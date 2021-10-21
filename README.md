# Quantitative Genetic LAVaan Models
`qglavmods` is an R-package to generate syntax for quantitative genetic
twin models.

Currently only ACE forms of univariate and bivariate models using the correlated factors solusion is supported
Current limitations:
- Binary/categorical outcomes are treated as continuous
- Only supporting ACE versions
- Identification only through standardization of latent variables
## Usage
```
# Install package
devtools::install_github('deepchocolate/qglavmods')
library(qglavmods)

# Create univariate model
modA <- univTwinModel()
# Add regression (intercept + female sex)
modA <- regress(modA, measure ~ 1 + female)
fitA < lavaan(as.character(modA), data=dataInWideFromat, gropu='zygocityColumn', group.label=c('MZ','DZ'), std.lv=T)
# Create a univariate model for a different trait
modB <- regress(modA, measure ~ 1 + female + somethingElse)
fitB < lavaan(as.character(modB), data=dataInWideFromat, gropu='zygocityColumn', group.label=c('MZ','DZ'), std.lv=T)

# Now do a bivariate version
modAB <- bivTwinModel(modA, modB)
fitAB < lavaan(as.character(modAB), data=dataInWideFromat, gropu='zygocityColumn', group.label=c('MZ','DZ'), std.lv=T)
```