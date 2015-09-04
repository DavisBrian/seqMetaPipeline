context("phenotype")

data(pheno1)

test_that("formula checking", {
  expect_error(phenotype(pheno1, formula = "y~AGE", family = "gaussian", 
                         id = "id"), regexp = "Error in parameter 'formula'")
  expect_is(phenotype(pheno1, formula = "y~1", family = "gaussian", id = "id"), "phenotype")
  expect_is(phenotype(pheno1, formula = "y~1.", family = "gaussian", id = "id"), "phenotype")
})