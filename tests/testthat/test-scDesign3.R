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

  my_marginal1 <- fit_marginal(
    data = my_data,
    mu_formula = "1",
    sigma_formula = "cell_type",
    family_use = "nb",
    n_cores = 1,
    usebam = FALSE
  )

   my_marginal2 <- fit_marginal(
    data = my_data,
    mu_formula = "s(pseudotime, bs = 'cr', k = 10)",
    sigma_formula = "s(pseudotime, bs = 'cr', k = 3)",
    family_use = c(rep("nb", 5), rep("zip", 5)),
    n_cores = 1,
    usebam = FALSE
  )

  my_marginal3 <- fit_marginal(
    data = my_data,
    mu_formula = "s(pseudotime, bs = 'cr', k = 10)",
    sigma_formula = "s(pseudotime, bs = 'cr', k = 3)",
    family_use = c(rep("nb", 5), rep("zip", 5)),
    n_cores = 1,
    usebam = TRUE
  )

  my_copula <- fit_copula(
    sce = example_sce,
    assay_use = "counts",
    marginal_list = my_marginal3,
    family_use = c(rep("nb", 5), rep("zip", 5)),
    copula = "vine",
    n_cores = 1,
    new_covariate = NULL,
    input_data = my_data$dat
  )

  my_copula1 <- fit_copula(
    sce = example_sce,
    assay_use = "counts",
    marginal_list = my_marginal3,
    family_use = c(rep("nb", 5), rep("zip", 5)),
    copula = "gaussian",
    n_cores = 1,
    new_covariate = NULL,
    input_data = my_data$dat
  )

  my_para <- extract_para(
    sce = example_sce,
    marginal_list = my_marginal3,
    n_cores = 1,
    family_use = c(rep("nb", 5), rep("zip", 5)),
    new_covariate = NULL
  )

  my_newcount <- simu_new(
    sce = example_sce,
    mean_mat = my_para$mean_mat,
    sigma_mat = my_para$sigma_mat,
    zero_mat = my_para$zero_mat,
    quantile_mat = my_copula$new_mvu,
    n_cores = 1,
    family_use = c(rep("nb", 5), rep("zip", 5)),
    new_covariate = my_data$new_covariate
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
    n_cores = 1,
    usebam = FALSE,
    corr_formula = "pseudotime",
    copula = "vine",
    DT = TRUE,
    pseudo_obs = FALSE,
    ncell = 1000,
    return_model = TRUE
  )

})

