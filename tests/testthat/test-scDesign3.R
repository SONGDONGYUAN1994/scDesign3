#context("Run scDesign3")
library(scDesign3)

test_that("Run scDesign3", {
  data(example_sce)

  my_data <- construct_data(sce = example_sce,
                            assay_use = "counts",
                            covariate_use = "pseudotime",
                            celltype = "cell_type",
                            pseudotime = "pseudotime",
                            spatial = NULL,
                            other_covariates = NULL,
                            corr_by = "pseudotime")

  my_marginal <- fit_marginal(data = my_data,
                              mu_formula = "s(pseudotime, bs = 'cr', k = 10)",
                              sigma_formula = "s(pseudotime, bs = 'cr', k = 5)",
                              family = "nb",
                              n_cores = 1,
                              usebam = FALSE)

  my_copula <- fit_copula(sce = example_sce,
                          assay_use = "counts",
                          marginal_list = my_marginal,
                          family = "nb",
                          copula = "vine",
                          cor_formula = "pseudotime",
                          n_cores = 1,
                          new_covariate = NULL,
                          input_data = my_data$dat)

  my_newcount <- simu_new(sce = example_sce,
                          marginal_list = my_marginal,
                          copula_list = my_copula,
                          n_cores = 1,
                          family = "nb",
                          new_covariate = my_data$new_covariate)

  my_simu <- scdesign3(sce = example_sce,
                       assay_use = "counts",
                       covariate_use = "pseudotime",
                       celltype = "cell_type",
                       pseudotime = "pseudotime",
                       spatial = NULL,
                       other_covariates = NULL,
                       predictor = "gene",
                       mu_formula = "s(pseudotime, bs = 'cr', k = 10)",
                       sigma_formula = "s(pseudotime, bs = 'cr', k = 5)",
                       family = "nb",
                       n_cores = 1,
                       usebam = FALSE,
                       cor_formula = "pseudotime",
                       copula = "vine",
                       DT = TRUE,
                       pseudo_obs = FALSE,
                       ncell = 1000)

})

