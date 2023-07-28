#context("Run scDesign3")
library(scDesign3)

test_that("Run scDesign3", {
  data(example_sce)

  my_data <- construct_data(
    sce = example_sce,
    assay_use = "counts",
    celltype = "cell_type",
    pseudotime = "pseudotime",
    spatial = NULL,
    other_covariates = NULL,
    corr_by = "cell_type"
  )

  my_data2 <- construct_data(
    sce = example_sce,
    assay_use = "counts",
    celltype = "cell_type",
    pseudotime = "pseudotime",
    spatial = NULL,
    other_covariates = NULL,
    corr_by = "pseudotime",
    ncell = 10000
  )

  my_marginal1 <- fit_marginal(
    data = my_data,
    mu_formula = "1",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 1,
    usebam = FALSE
  )

   my_marginal2 <- fit_marginal(
    data = my_data,
    mu_formula = "s(pseudotime, bs = 'cr', k = 10)",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 1,
    usebam = FALSE
  )

  my_fit1 <- lapply(my_marginal1, function(x)x$fit)
  my_fit2 <- lapply(my_marginal2, function(x)x$fit)

  my_pvalue <- perform_lrt(my_fit2, my_fit1)

  my_marginal3 <- fit_marginal(
    data = my_data,
    mu_formula = "s(pseudotime, bs = 'cr', k = 10)",
    sigma_formula = "s(pseudotime, bs = 'cr', k = 3)",
    family_use = c(rep("nb", 5), rep("zip", 5)),
    n_cores = 2,
    usebam = TRUE
  )

  my_copula <- fit_copula(
    sce = example_sce,
    assay_use = "counts",
    marginal_list = my_marginal3,
    family_use = c(rep("nb", 5), rep("zip", 5)),
    copula = "vine",
    n_cores = 2,
    input_data = my_data$dat
  )

  my_quantile_mat <- fit_copula(
    sce = example_sce,
    assay_use = "counts",
    empirical_quantile = TRUE,
    marginal_list = my_marginal3,
    family_use = c(rep("nb", 5), rep("zip", 5)),
    copula = "vine",
    n_cores = 2,
    input_data = my_data$dat
  )
  
  # my_copula1 <- fit_copula(
  #   sce = example_sce,
  #   assay_use = "counts",
  #   marginal_list = my_marginal3,
  #   family_use = c(rep("nb", 5), rep("zip", 5)),
  #   copula = "gaussian",
  #   n_cores = 1,
  #   new_covariate = NULL,
  #   input_data = my_data$dat
  # )

  my_para <- extract_para(
    sce = example_sce,
    marginal_list = my_marginal3,
    n_cores = 1,
    family_use = c(rep("nb", 5), rep("zip", 5)),
    new_covariate = my_data$new_covariate,
    data = my_data$dat
  )

  my_newcount <- simu_new(
    sce = example_sce,
    mean_mat = my_para$mean_mat,
    sigma_mat = my_para$sigma_mat,
    zero_mat = my_para$zero_mat,
    quantile_mat = NULL,
    copula_list = my_copula$copula_list,
    n_cores = 1,
    family_use = c(rep("nb", 5), rep("zip", 5)),
    input_data = my_data$dat,
    new_covariate = my_data$new_covariate,
    important_feature = my_copula$important_feature
  )
  
  my_newcount2 <- simu_new(
    sce = example_sce,
    mean_mat = my_para$mean_mat,
    sigma_mat = my_para$sigma_mat,
    zero_mat = my_para$zero_mat,
    quantile_mat = my_quantile_mat$quantile_mat,
    copula_list = NULL,
    n_cores = 1,
    family_use = c(rep("nb", 5), rep("zip", 5)),
    input_data = my_data$dat,
    new_covariate = my_data$new_covariate,
    important_feature = my_copula$important_feature
  )

  my_simu <- scdesign3(
    sce = example_sce,
    assay_use = "counts",
    celltype = "cell_type",
    pseudotime = "pseudotime",
    spatial = NULL,
    other_covariates = NULL,
    mu_formula = "s(pseudotime, bs = 'cr', k = 10)",
    sigma_formula = "s(pseudotime, bs = 'cr', k = 3)",
    family_use = c(rep("nb", 5), rep("zip", 5)),
    n_cores = 2,
    usebam = FALSE,
    corr_formula = "pseudotime",
    copula = "vine",
    DT = TRUE,
    pseudo_obs = FALSE,
    ncell = 1000,
    return_model = TRUE
  )

  # my_simu2 <- scdesign3(
  #   sce = example_sce,
  #   assay_use = "counts",
  #   celltype = "cell_type",
  #   pseudotime = "pseudotime",
  #   spatial = NULL,
  #   other_covariates = NULL,
  #   mu_formula = "s(pseudotime, bs = 'cr', k = 10)",
  #   sigma_formula = "s(pseudotime, bs = 'cr', k = 3)",
  #   family_use = c(rep("nb", 5), rep("zip", 5)),
  #   n_cores = 2,
  #   usebam = FALSE,
  #   corr_formula = "pseudotime",
  #   empirical_quantile = TRUE,
  #   copula = "vine",
  #   DT = TRUE,
  #   pseudo_obs = FALSE,
  #   return_model = FALSE
  # )
})

