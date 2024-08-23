###### Function to run MR ######



# #' Run MR
# #'
# #' Use TwoSampleMR to perform IVW-MR and returns causal effect estimate (SE)
# #' but also M (number of instruments), mean sample size for exposure / outcome
# #' and the set of IVs used
# #'
# #' @param exposure_data xx
# #' @param outcome_data xx
# #'
# #' @inheritParams MRlap
# #' @export
# NOT EXPORTED


run_MR <- function(exposure_data,
                   outcome_data,
                   do_pruning = TRUE,
                   user_SNPsToKeep = "",
                   MR_threshold = 5e-8,
                   MR_pruning_dist = 500,
                   MR_pruning_LD = 0,
                   MR_reverse = NULL,
                   MR_plink = NULL,
                   MR_bfile = NULL,
                   verbose = TRUE){

  # here we need to join exposure and outcome data (merge only based on rsid - chr/pos are optional)
  data <- dplyr::inner_join(exposure_data, outcome_data, by=c("rsid"), suffix=c(".exp", ".out"))
  # and we need to make sure alleles are aligned!
  data %>%
    dplyr::mutate(std_beta.out = dplyr::case_when(
      alt.exp == alt.out & ref.exp == ref.out ~ std_beta.out, # aligned
      ref.exp == alt.out & alt.exp == ref.out ~ -std_beta.out, # swapped
      TRUE ~ NA_real_ # otherwise exclude SNP
    ) ) %>%
    dplyr::filter(!is.na(std_beta.out)) -> data




  # doing pruning? Or using user-provided SNPsToKeep?
  if (!do_pruning) {  # use user-provided SNPsToKeep
    if(verbose) cat("> Using IVs provided by user... \n")
    SNPsToKeep = user_SNPsToKeep
    # add message if some SNPs missing in exposure / outcome data
    if(verbose & any(!SNPsToKeep %in% data$rids)) cat("some user-provided IVs are missing from the merged exposure-outcome dataset and will not be used\n")
    # get "pruned" data
    data %>%
      dplyr::filter(.data$rsid %in% SNPsToKeep) -> data_pruned

    if(verbose) cat("   ", format(nrow(data_pruned), big.mark=","), "IVs left \n")

  } else {
    if(verbose) cat("> Identifying IVs... \n")


    data %>%
      dplyr::filter(.data$p.exp<MR_threshold) -> data_thresholded
    if(verbose) cat("   ", format(nrow(data_thresholded), big.mark=","), "IVs with p <", format(MR_threshold, scientific = T), "\n")

    # if no IVs, stop()
    if(nrow(data_thresholded)==0) stop("no IV left after thresholding")

    # here remove the ones with
    if(!is.null(MR_reverse)){
      reverse_t_threshold  =  stats::qnorm(MR_reverse)
      data_thresholded %>%
        dplyr::filter( ( abs(.data$std_beta.exp) - abs(.data$std_beta.out)) /
                         sqrt(.data$std_SE.exp^2 + .data$std_SE.out^2) > reverse_t_threshold) -> data_thresholded_filtered
      if(verbose) cat("   ", nrow(data_thresholded)-nrow(data_thresholded_filtered), "IVs excluded - more strongly associated with the outcome than with the exposure, p <", format(MR_reverse, scientific = T), "\n")
      data_thresholded = data_thresholded_filtered
      rm(data_thresholded_filtered)
    }

    # if no IVs, stop()
    if(nrow(data_thresholded)==0) stop("no IV left after excluding IVs more strongly associated with the outcome than with the exposure")

    if(MR_pruning_LD>0){# LD-pruning
      data_thresholded %>%
        dplyr::transmute(SNP = .data$rsid,
                         chr_name = .data$chr.exp,
                         chr_start = .data$pos.exp,
                         pval.exposure = .data$p.exp) -> ToPrune

      if(verbose) cat("   Pruning : distance : ", MR_pruning_dist, "Kb", " - LD threshold : ", MR_pruning_LD, "\n")
      # Do pruning
      SNPsToKeep = c()
      # use remote clumping?
      if(is.null(MR_plink)){
        for(chr in unique(ToPrune$chr_name)){
          SNPsToKeep = c(SNPsToKeep, suppressMessages(TwoSampleMR::clump_data(ToPrune[ToPrune$chr_name==chr,], clump_kb = MR_pruning_dist, clump_r2 = MR_pruning_LD)$SNP))
        }
      } else { # use local clumping
        if(verbose) cat("Using ieugwasr::ld_clump with local PLINK")
        # need slightly different headers
        ToPrune = ToPrune %>% dplyr::select(rsid = SNP, #chr = chr_name, pos = chr_start,
                                            pval = pval.exposure)
        SNPsToKeep = ieugwasr::ld_clump(ToPrune,
                                        clump_kb = MR_pruning_dist,
                                        clump_r2 = MR_pruning_LD,
                                        plink_bin = MR_plink,
                                        bfile = MR_bfile)$rsid
      }
    } else{ # distance pruning
      data_thresholded %>%
        dplyr::transmute(SNP = .data$rsid,
                         chr_name = .data$chr.exp,
                         chr_start = .data$pos.exp,
                         pval.exposure = .data$p.exp) -> ToPrune
      prune_byDistance <- function(data, prune.dist=100, byP=T) {
        # data should be : 1st column rs / 2nd column chr / 3rd column pos / 4th column stat
        # if byP = T : stat = p-value -> min is better
        # if byP = F : stat = Zstat, beta.. -> max is better

        if(byP){
          SNP_order = order(data %>% dplyr::pull(4))
        } else {
          SNP_order = order(data %>% dplyr::pull(4), decreasing = T)
        }
        data = data[SNP_order,]
        snp=0
        while(T){
          snp=snp+1
          ToRemove=which(data$chr_name==data$chr_name[snp] & abs(data$chr_start - data$chr_start[snp])<prune.dist*1000)
          if(length(ToRemove)>1){
            ToRemove = ToRemove[-1]
            data = data[-ToRemove,]
          }
          if(snp==nrow(data)) break
        }

        return(unlist(data[,1]))

      }
      if(verbose) cat("   Pruning : distance : ", MR_pruning_dist, "Kb \n")
      SNPsToKeep = prune_byDistance(ToPrune, prune.dist=MR_pruning_dist, byP=T)
    }


  # get "pruned" data
    data_thresholded %>%
      dplyr::filter(.data$rsid %in% SNPsToKeep) -> data_pruned

    if(verbose) cat("   ", format(nrow(data_pruned), big.mark=","), "IVs left after pruning \n")
  }

  ## check, if IVs < 2, fail?


  if(verbose) cat("> Performing MR \n")

  TwoSampleMR::mr_ivw(data_pruned$std_beta.exp, data_pruned$std_beta.out,
                      data_pruned$std_SE.exp, data_pruned$std_SE.out) -> res_MR_TwoSampleMR

  TwoSampleMR::mr_egger_regression(data_pruned$std_beta.exp, data_pruned$std_beta.out,
                                   data_pruned$std_SE.exp, data_pruned$std_SE.out) -> res_MR_TwoSampleMR_Egger

  if(verbose) cat("   ",  "IVW-MR observed effect:", format(res_MR_TwoSampleMR$b, digits = 3), "(", format(res_MR_TwoSampleMR$se, digits=3), ")\n")

  return(list("alpha_obs" = res_MR_TwoSampleMR$b,
              "alpha_obs_se" = res_MR_TwoSampleMR$se,
              "egger_b" = res_MR_TwoSampleMR_Egger$b,
              "egger_se" = res_MR_TwoSampleMR_Egger$se,
              "egger_intercept_p" = res_MR_TwoSampleMR_Egger$pval_i,
              "n_exp" = mean(data_pruned$N.exp),
              "n_out" = mean(data_pruned$N.out),
              "IVs" = data_pruned %>% dplyr::select(.data$std_beta.exp, .data$std_SE.exp),
              "IVs_rs" = data_pruned$rsid))

}
