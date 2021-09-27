###### Function to get correction ######



# #' Get correction
# #'
# #' Calculate corrected causal effect estimate (SE) and the covariance between
# #' observed and corrected effects, and test for the difference between the two
# #'
# #' @param IVs xx
# #' @param lambda xx
# #' @param lambda_se xx
# #' @param h2_LDSC xx
# #' @param h2_LDSC_se xx
# #' @param alpha_obs xx
# #' @param alpha_obs_se xx
# #' @param n_exp xx
# #' @param n_out xx
# #' @param M xx
# #'
# #' @inheritParams MRlap
# #' # @export
# NOT EXPORTED

get_correction <- function(IVs, lambda, lambda_se, h2_LDSC, h2_LDSC_se,
                           alpha_obs, alpha_obs_se, n_exp, n_out, MR_threshold, verbose, s=100, sthreshold, extracheck=F){

  M=1150000 # consider M is a constant! number of independent markers genome-wide
  Tr = -stats::qnorm(MR_threshold/2)
  lambdaPrime = lambda/sqrt(n_exp*n_out)

  if(verbose) cat("> Estimating genetic architecture parameters... \n")

  # function for optimisation
  get_pi <- function(my_pi, sumbeta2, Tr, n_exp, h2_LDSC, M){
    if(0>=my_pi) return(1e6)

    sigma = sqrt(h2_LDSC/(my_pi*M))

    denominator = (my_pi * (2*(sigma^2 + 1/n_exp) * stats::pnorm( - Tr / sqrt( 1 + n_exp * sigma^2) ) +
                              2 * Tr *(n_exp * sigma^4 + 2*sigma^2 + 1/n_exp) *exp(-Tr^2 / (2*(n_exp*sigma^2+1)))/( sqrt(2*pi) * ( 1 + n_exp *sigma^2)^(3/2))) +
                     (1-my_pi)*1/n_exp * (2*stats::pnorm(-Tr) + 2*Tr *stats::dnorm(Tr))) * M

    return(abs(denominator-sumbeta2))
  }

  # get genetic architecture
  get_geneticArchitecture<- function(theta, Nexp, M, Tr){
    # theta is (effects, h2_LDSC)
    h2_LDSC = theta[length(theta)]
    effects = theta[-length(theta)]
    sumBeta2 = sum(effects^2)
    nSP = 5 # seems enough! convergence is very good, for all starting points
    Res_SP = data.frame(SP = 1:nSP,
                        SP_pi = NA_real_,
                        diff = NA_real_,
                        pi = NA_real_)

    for(i in 1:nSP){
      theta = 3 * 10^(stats::runif(1, -7, -1))
      res_optim = stats::optimise(get_pi, interval=c(1e-7, 0.3),
                                  sumbeta2=sumBeta2, Tr=Tr, n_exp=n_exp, h2_LDSC=h2_LDSC, M=M, tol = 1e-6)

      Res_SP[i, 2:4] = c(theta, res_optim$objective, res_optim$minimum)
    }

    # take best one
    Res_SP %>% dplyr::arrange(diff) %>% dplyr::slice(1) %>% dplyr::pull(pi) -> pi_x
    sigma = sqrt(h2_LDSC/(M*pi_x))


    return(c(as.numeric(pi_x), as.numeric(sigma)))
  }
  theta = c(IVs$std_beta.exp,h2_LDSC)
  res_genA = get_geneticArchitecture(theta, n_exp, M, Tr)

  get_alpha <- function(n_exp, lambdaPrime, pi_x, sigma, alpha_obs, Tr){

    A = stats::pnorm( - Tr / sqrt(1 + n_exp * sigma^2) )
    B = 2 * Tr * exp(-Tr^2 / (2*(n_exp*sigma^2+1)))/( sqrt(2*pi)* ( 1 + n_exp *sigma^2)^(3/2))
    C = (stats::pnorm(-Tr) + Tr *stats::dnorm(Tr))

    a = (pi_x * (2*(sigma^2 + 1/n_exp) * A +
                   B *(n_exp * sigma^4 + 2*sigma^2 + 1/n_exp)))
    b =  (1-pi_x)*2/n_exp * C


    d = a+b

    alpha = (alpha_obs * d - (lambdaPrime * pi_x * (2*A + B + sigma^2 * n_exp*B) +
                                lambdaPrime * (1 - pi_x) * 2 * C)) / ( pi_x * sigma^2 * (2*A + B * n_exp * sigma^2 + B))
    return(alpha)
  }

  if(verbose) cat("> Estimating corrected effect... \n")

  alpha_corrected = get_alpha(n_exp, lambdaPrime, res_genA[1], res_genA[2], alpha_obs, Tr)



  ## get SE and covariance
  get_correctedSE <- function(IVs, lambda, lambda_se, h2_LDSC, h2_LDSC_se, alpha_obs, alpha_obs_se, n_exp, n_out, M, Tr, s=1000, sthreshold, extracheck){

    get_s <- function(s){
      effects = IVs$std_beta.exp
      effects_se = IVs$std_SE.exp

      # simulate 500 lambda
      L = matrix(stats::rnorm(s, lambda, lambda_se), ncol=s)/sqrt(n_exp*n_out)

      # simulate 500 "instruments sets" - each column = 1 simulation
      E =  matrix(stats::rnorm(nrow(IVs)*s, effects, effects_se), ncol= s)

      # simulate 500 h2_LDSC
      H = matrix(stats::rnorm(s, h2_LDSC, h2_LDSC_se), ncol=s)

      # effects + h2 are needed to get pi and therefore sigma
      D = rbind(E, H)

      pis = apply(D, 2, function(x) get_geneticArchitecture(x, n_exp, M, Tr))

      # simulate 500 alpha_obs
      B =  matrix(stats::rnorm(s, alpha_obs, alpha_obs_se), ncol= s)



      # get 500 alpha_corrected + sd
      all_params = data.frame(pi = pis[1,],
                              sigma = pis[2,],
                              alpha = B[1,],
                              lambda = L[1,])


      all_params$corrected = apply(all_params, 1, function(x) get_alpha(n_exp, x[4],  x[1], x[2], x[3], Tr))
      return(all_params)
    }

    res = get_s(s)

    num_groups=10
    res %>%
      dplyr::group_by((row_number()-1) %/% (n()/num_groups)) %>%
      tidyr::nest() %>% pull(data) -> res_subsets
    subsets_se = unlist(lapply(res_subsets, function(x) stats::sd(x$corrected)))
    subsets_cov = unlist(lapply(res_subsets,  function(x) stats::cov(x$corrected, x$alpha)))

    needmore = F
    if(stats::sd(subsets_se) / base::mean(subsets_se) > sthreshold) needmore=T
    if(stats::sd(subsets_cov) / base::mean(subsets_cov) > sthreshold) needmore=T
    if(extracheck & (alpha_obs_se^2 + se_cov[1]^2 - 2* se_cov[2])<0) needmore=T

    temp = data.frame(s=nrow(res),
                      SE =  stats::sd(res$corrected),
                      ratioSE = stats::sd(subsets_se) / base::mean(subsets_se),
                      COV = stats::cov(res$corrected, res$alpha),
                      ratioCOV = stats::sd(subsets_cov) / base::mean(subsets_cov))

    while(needmore){
      res = rbind(res,get_s(s))
      res %>%
        dplyr::group_by((row_number()-1) %/% (n()/num_groups)) %>%
        tidyr::nest() %>% pull(data) -> res_subsets
      subsets_se = unlist(lapply(res_subsets, function(x) stats::sd(x$corrected)))
      subsets_cov = unlist(lapply(res_subsets,  function(x) stats::cov(x$corrected, x$alpha)))

      needmore = F
      if(stats::sd(subsets_se) / base::mean(subsets_se) > sthreshold) needmore=T
      if(stats::sd(subsets_cov) / base::mean(subsets_cov) > sthreshold) needmore=T
      if(extracheck & (alpha_obs_se^2 + se_cov[1]^2 - 2* se_cov[2])<0) needmore=T


      temp = rbind(temp,
                   data.frame(s=nrow(res),
                              SE =  stats::sd(res$corrected),
                              ratioSE = stats::sd(subsets_se) / base::mean(subsets_se),
                              COV = stats::cov(res$corrected, res$alpha),
                              ratioCOV = stats::sd(subsets_cov) / base::mean(subsets_cov)))
    }

    all_res = list(
      "res" = c(stats::sd(res$corrected), stats::cov(res$alpha, res$corrected), nrow(res)),
      "extensive_res" = temp)
    # normally, remove all temp and return just "res"
    return(all_res)
  }

  # get SE corrected effects + COV
  se_cov = get_correctedSE(IVs, lambda, lambda_se, h2_LDSC, h2_LDSC_se, alpha_obs, alpha_obs_se, n_exp, n_out, M, Tr,s)
  tmp = se_cov$extensive_res
  se_cov = se_cov$res

  if(verbose) cat("   ",  "corrected effect:", format(alpha_corrected, digits = 3), "(", format(se_cov[1], digits=3), ")\n")
  if(verbose) cat("   ",  "covariance between observed and corrected effect:", format(se_cov[2], digits=3), "\n")
  if(verbose) cat("           ",  se_cov[3], " simulations were used to estimate the variance and the covariance.\n")

  if(verbose) cat("> Testing difference between observed and corrected effect... \n")

  test_diff = (alpha_obs - alpha_corrected) /
    sqrt(alpha_obs_se^2 + se_cov[1]^2 - 2* se_cov[2])
  p_diff = 2*stats::pnorm(-abs(test_diff), lower.tail=T)
  return(list("alpha_corrected"=alpha_corrected,
              "alpha_corrected_se" = se_cov[1],
              "cov_obs_corrected" = se_cov[2],
              "test_diff"=test_diff,
              "p_diff"=p_diff,
              "pi_x" = res_genA[1],
              "sigma2_x" = res_genA[2]^2,
              "temp" = temp))
}
