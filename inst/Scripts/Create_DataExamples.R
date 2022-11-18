# to run on laptop #

set.seed(33)

library(tidyverse)

Exposure_Data <- data.table::fread("~/Documents/SGG/Projects/SampleOverlap/Data/Simulations/noY/GWAS_X.tsv")
Outcome_Data <- data.table::fread("~/Documents/SGG/Projects/SampleOverlap/Data/Simulations/noY/GWAS_Y100.tsv")

# take 750,000 SNPs
SNPs <- sample(x = 1:nrow(Exposure_Data), 750000)

Exposure_Data %>%
  slice(SNPs) -> SmallExposure_Data

Outcome_Data %>%
  slice(SNPs) -> SmallOutcome_Data

save(SmallExposure_Data, file="~/Documents/SGG/Projects/MRlap/data/SmallExposure_Data.rda", compress='xz')
save(SmallOutcome_Data, file="~/Documents/SGG/Projects/MRlap/data/SmallOutcome_Data.rda", compress='xz')




## Rdata for tests
# last update, 2022/11/18 (low h2 check)
library(MRlap)

# we use ~100K samples for BMI/SBP, with 0% of sample overlap
# (only weak-instrument bias and Winner's curse)
BMI <- system.file("data/", "BMI_Data.tsv.gz", package="MRlap")
SBP <- system.file("data/", "SBP_Data.tsv.gz", package="MRlap")

A = MRlap(exposure = BMI,
          exposure_name = "BMI_100Ksample",
          outcome = SBP,
          outcome_name = "SBP_100Ksample",
          ld = "~/eur_w_ld_chr",
          hm3 = "~/w_hm3.noMHC.snplist")

saveRDS(A, file="~/Documents/SGG/Projects/MRlap/inst/Data/A.RDS")



# we use simulated data (standard settings scenario), with 100% of sample overlap

data("SmallExposure_Data")
data("SmallOutcome_Data")

B = MRlap(exposure = SmallExposure_Data,
          exposure_name = "simulated_exposure",
          outcome = SmallOutcome_Data,
          outcome_name = "simulated_outcome",
          ld = "~/eur_w_ld_chr",
          hm3 = "~/w_hm3.noMHC.snplist",
          MR_threshold = 5e-10,
          MR_pruning_LD = 0.05)


saveRDS(B, file="~/Documents/SGG/Projects/MRlap/inst/Data/B.RDS")




