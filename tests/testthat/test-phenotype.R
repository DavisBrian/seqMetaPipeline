context("phenotype")

data(pheno1)

test_that("formula checking)", {
  expect_error(phenotype(pheno1, formula = "y ~ AGE", family = "gaussian", 
                         id = "id"), 
               regexp = "One or more formula variables not in data")
  expect_is(phenotype(pheno1, formula = "y~1", family = "gaussian", id = "id"), 
            "phenotype")
  expect_is(phenotype(pheno1, formula = "y~.", family = "gaussian", id = "id"), 
            "phenotype")
  expect_is(phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                      family = "gaussian", id = "id"), "phenotype")
})

test_that("family checking)", {
  expect_is(phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                      family = "gaussian", id = "id"), "phenotype")
  expect_is(phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                      family = "binomial", id = "id"), "phenotype")
  expect_is(phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                      family = "cox", id = "id"), "phenotype")
  ## check missing family
  expect_error(phenotype(pheno1, formula = "y ~ age + pc1 + pc2", id = "id"), 
               regexp = "'family' must be specified!")

  ## check not character
  expect_error(phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                      family = gaussian(), id = "id"), 
               regexp = "'family' must be of type character")
  
  ## check invaild type
  expect_error(phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                         family = "abc", id = "id"), 
               regexp = "Only family 'gaussian', 'binomial', or 'cox' currently supported")
  
  expect_error(phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                         family = character(0), id = "id"), 
               regexp = "'family' must be specified!")
  
  expect_error(phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                         family = NULL, id = "id"), 
               regexp = "'family' must be specified!")
  
  ## multiple family type specified
  expect_error(phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                         family = c("gaussian", "gaussian"), id = "id"), 
               regexp = "Family can be only one of: 'gaussian', 'binomial', or 'cox'")
})

test_that("family id)", {
  # id col specified
  expect_is(pheno <- phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                      family = "gaussian", id = "id"), "phenotype")
  expect_equal(attr(pheno, "idCol"), "id")
  
  # id col mis specified
  expect_error(phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                         family = "gaussian", id = "gwas_ids"), 
               regexp = "is not a column in data")  
  
  # id col not specified
  expect_warning(p <- phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                                family = "gaussian"), 
                 regexp = "Using row names")
  expect_equal(attr(p, "idCol"), ".id")
  expect_equal(p[ , ".id"], rownames(pheno1))
  
  # multiple columns given
  expect_error(phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                         family = "gaussian", id = c("id", "y")), 
               regexp = "Only one id column can be specified.") 
  
  
  p2 <- pheno1[ , -1]  
  # column not there
  expect_error(phenotype(p2, formula = "y ~ age + pc1 + pc2", 
                      family = "gaussian", id = "id"), 
               regexp = "is not a column in data")
 

  p3 <- pheno1
  p3[floor(nrow(p3)/2), "id"] <- p3[1, "id"]
  
  # non unique id's
  expect_error(phenotype(p3, formula = "y ~ age + pc1 + pc2", 
                         family = "gaussian", id = "id"), 
               regexp = "Duplicate phenotype ids found")
  
  # missing id's
  p3[floor(nrow(p3)/2), "id"] <- NA
  expect_error(phenotype(p3, formula = "y ~ age + pc1 + pc2", 
                         family = "gaussian", id = "id"), 
               regexp = "Missing phenotype ids")
  
  # convert factor to character
  p4 <- pheno1
  p4$id <- as.factor(p4$id)
  expect_warning(phenotype(p4, formula = "y ~ age + pc1 + pc2", 
                         family = "gaussian", id = "id"), 
               regexp = "Converting id column to character")  
})