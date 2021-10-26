# MRlap 0.0.1.1 (2021-10-26)

## Bug fixes
- Excluding IVs more strongly associated with the outcome than with the exposure
There was a typo leading to an underestimation of the t-statistics, and IVs more strongly associated with the outcome than with the exposure were not correctly identified.

# MRlap 0.0.1.0 (2021-10-15)

## Changes
- Case-control GWAS    

It is now possible to use case-control GWAS for the exposure and/or the outcome datasets.
If (at least) one of the datasets is coming from a case-control GWAS, the Sample size column should correspond to the total sample size (not the effective sample size).   
The analysis is performed on the observed scale (explaining why heritability estimates will differ from the ones usually obtained on the liability scale) and the approach is the same as the one used for continuous traits.    

- optimal number of simulations to estimate SE/COV

The number of simulation used to estimate the standard error of the corrected effect and the covariance between corrected and observed effect is no longer a parameter that needs to be chosen before running the analysis. This could be a problem if the value was too small, the SE/COV estimates were too noisy to correcly test for the difference between observed and corrected effects. The optimal value is now automatically determined when running the analysis (using the coefficients of variation for the variance of the corrected effect and covariance between corrected and observed effect). 

- p-value estimation for observed and corrected effects / number of IVs

The p-values are now reported for both observed and corrected effects, as well as the number of instruments used.

# MRlap 0.0.0.9 (2021-03-04)

initial version of the package, as described in the preprint

<!--- 
## Bug fixes

## New functions

## Documentation

## Error messages

## Performance


--->  
