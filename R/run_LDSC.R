###### Function to run LDSC ######



# #' Run LDSC
# #'
# #' Use GenomicSEM to perform cross-trait LDSC analysis and returns
# #' the heritability of the exposure (SE) and the cross-trait intercept (SE)
# #'
# #' @param exposure_data xx
# #' @param outcome_data xx
# #'
# #' @inheritParams MRlap
# # #' @export
# NOT EXPORTED

run_LDSC <- function(exposure_data,
                     exposure_name,
                     K_exposure,
                     P_exposure,
                     outcome_data,
                     outcome_name,
                     K_outcome,
                     P_outcome,
                     ld, hm3, save_logfiles, verbose){

  # write down exposure/outcome data
  utils::write.table(exposure_data, paste0(exposure_name, ".tsv"), sep="\t", quote=F, row.names=F)
  utils::write.table(outcome_data, paste0(outcome_name, ".tsv"), sep="\t", quote=F, row.names=F)

  if(verbose) cat("> Munging exposure data... \n")
  invisible(utils::capture.output(GenomicSEM::munge( paste0(exposure_name, ".tsv"),
                                              hm3,
                                              exposure_name)))

  if(verbose & save_logfiles) cat("  Please check the log file", paste0(exposure_name, "_munge.log"),  "to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files\n")


  if(verbose) cat("> Munging outcome data... \n")
  invisible(utils::capture.output(GenomicSEM::munge( paste0(outcome_name, ".tsv"),
                                              hm3,
                                              outcome_name)))
  if(verbose & save_logfiles) cat("  Please check the log file", paste0(outcome_name, "_munge.log"),  "to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files\n")


  traits <- c(paste0(exposure_name, ".sumstats.gz"),
              paste0(outcome_name, ".sumstats.gz"))
  sample.prev <- c(P_exposure,P_outcome) # continuous traits
  population.prev <- c(K_exposure,K_outcome) # continuous traits

  trait.names<-c(exposure_name, outcome_name)


  if(verbose) cat("> Running cross-trait LDSC... \n")

  invisible(utils::capture.output(LDSCoutput <- GenomicSEM::ldsc(traits,
                                                          sample.prev,
                                                          population.prev,
                                                          ld,
                                                          trait.names)))

  if(verbose & save_logfiles) cat("  Please check the log file", paste0(c(traits, "ldsc.log"), collapse="_"),  "for detailed results of the cross-trait LDSC analysis\n")
  if(verbose & !save_logfiles) cat("  Please consider saving the log files and checking them to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files\n")


  # we need : h2 exposure & se + cross-trait intercept & se
  h2_exp = as.numeric(LDSCoutput$S[1,1])
  h2_exp_SE = sqrt(LDSCoutput$V[1,1])

  # if case-control, we need to rescale the intercept: ... * sqrt(Ntot/Neff)
  gcov_int = as.numeric(LDSCoutput$I[1,2])
  if(is.na(K_exposure) & !is.na(K_outcome)){
    N_tot_outcome = median(outcome_data$N)
    N_eff_outcome = 4 * N_tot_outcome * K_outcome * N_tot_outcome * (1 - K_outcome) / N_tot_outcome
    gcov_int = gcov_int * sqrt(N_tot_outcome / N_eff_outcome)
  } else if(!is.na(K_exposure) & is.na(K_outcome)){
    N_tot_exposure = median(exposure_data$N)
    N_eff_exposure = 4 * N_tot_exposure * K_exposure * N_tot_exposure * (1 - K_exposure) / N_tot_exposure
    gcov_int = gcov_int * sqrt(N_tot_exposure / N_eff_exposure)
  } else if(!is.na(K_exposure) & !is.na(K_outcome)){
    N_tot_exposure = median(exposure_data$N)
    N_eff_exposure = 4 * N_tot_exposure * K_exposure * N_tot_exposure * (1 - K_exposure) / N_tot_exposure
    N_tot_outcome = median(outcome_data$N)
    N_eff_outcome = 4 * N_tot_outcome * K_outcome * N_tot_outcome * (1 - K_outcome) / N_tot_outcome
    gcov_int = gcov_int * sqrt(N_tot_exposure / N_eff_exposure) * sqrt(N_tot_outcome / N_eff_outcome)
  }

  output = readLines(paste0(c(traits, "ldsc.log"), collapse="_"))
  gcov_int_SE = as.numeric(stringr::str_replace(stringr::str_split(output[stringr::str_detect(output, "Cross trait Intercept")], "\\(" )[[1]][2], "\\)", ""))

  # other parameters estimates to report
  int_exp = LDSCoutput$I[1,1]
  int_out = LDSCoutput$I[2,2]

  h2_out = as.numeric(LDSCoutput$S[2,2])
  h2_out_SE = sqrt(LDSCoutput$V[3,3])

  rgcov = as.numeric(LDSCoutput$S[1,2])
  # what is reported by GenomicSEM is sigma_g : the covariance, need to be rescaled to get the correlation
  rg = as.numeric(LDSCoutput$S[1,2]/sqrt(LDSCoutput$S[1,1]*LDSCoutput$S[2,2]))
  # what is reported by GenomicSEM is se(sigma_g)^2 with sigma_g being the covariance, not the correlation!
  rgcov_SE= sqrt(LDSCoutput$V[2,2])



  if(verbose) cat("> Cleaning temporary files... \n")
  # exposure.tsv
  # exposure_munge.log
  # exposure.sumstats.gz
  if(save_logfiles){
    file.remove(paste0(exposure_name, c(".tsv", ".sumstats.gz")))
  } else {
    file.remove(paste0(exposure_name, c(".tsv", "_munge.log", ".sumstats.gz")))
  }
  # outcome.tsv
  # outcome_munge.log
  # outcome.sumstats.gz
  if(save_logfiles){
    file.remove(paste0(outcome_name, c(".tsv", ".sumstats.gz")))
  } else {
    file.remove(paste0(outcome_name, c(".tsv", "_munge.log", ".sumstats.gz")))
  }
  # exposure.sumstats.gz_outcome.sumstats.gz_ldsc.log
  if(!save_logfiles) file.remove(paste0(c(traits, "ldsc.log"), collapse="_"))

  return(list("h2_LDSC" = h2_exp,
           "h2_LDSC_se" = h2_exp_SE,
           "lambda" = gcov_int,
           "lambda_se" = gcov_int_SE,
           "int_exp" = int_exp,
           "int_out" = int_out,
           "h2_out" = h2_out,
           "h2_out_se" = h2_out_SE,
           "rgcov" = rgcov,
           "rgcov_se" = rgcov_SE,
           "rg" = rg))
}
