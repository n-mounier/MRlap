# MRlap 0.0.1.0 (2021-10-05)

## Changes
- Case-control GWAS    

It is now possible to use case-control GWAS for the exposure and/or the outcome datasets.
If (at least) one of the datasets is coming from a case-control GWAS, the Sample size column should correspond to the effective sample size (not the total sample size). `Ncases` (number of cases) and `Ncontrols` can also be provided (instead or in addition to the effective sample size).   
If the data has been analyzed using a linear model, there are two options:    
... if the sample prevalence and the population prevalence are provided, the analysis will be performed on the liability scale,    
...if the sample prevalence and the population prevalence are not provided, the analysis will be perfomed on the observed scale (sililarly to what is done for continuous traits).   
If the data has been analyzed using a logisitic model, the sample prevalence and the population prevalence need to be provided and the analysis will be performed on the liability scale. Additionally, in this case, if the Z-statistics is not present, the effect size column can either correspond to the odds ratio (`OR`) or to the log odds ratio  (`b` or `beta`).

- optimal number of simulations to estimate SE/COV

The number of simulation used to estimate the standard error of the corrected effect and the covariance between corrected and observed effect is no longer a parameter that needs to be chosen before running the analysis. This could be a problem if the value was too small, the SE/COV estimates were too noisy to correcly test for the difference between observed and corrected effects. The optimal value is now automatically determined when running the analysis. 

- p-value estimation for observed and corrected effects

The p-values are now reported for both observed and corrected effects.

# MRlap 0.0.0.9 (2021-03-04)

initial version of the package, as described in the preprint

<!--- 
## Bug fixes

## New functions

## Documentation

## Error messages

## Performance


--->  
