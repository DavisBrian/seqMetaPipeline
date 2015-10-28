context("phenotype")

data(pheno1)

pheno2 <- pheno1
pheno2$y <- sample(c(0,1), size = nrow(pheno1), replace = TRUE)

test_that("formula checking", {
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

test_that("family checking", {
  expect_is(phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                      family = "gaussian", id = "id"), "phenotype")

  expect_is(phenotype(pheno2, formula = "y ~ age + pc1 + pc2", 
                      family = "binomial", id = "id"), "phenotype")
  
#   expect_is(phenotype(pheno3, formula = "y ~ age + pc1 + pc2", 
#                       family = "cox", id = "id"), "phenotype")   

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

test_that("id checking", {
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

test_that("gender checking", {
  # gender col missing
  expect_is(pheno <- phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                               family = "gaussian", id = "id"), "phenotype")
  expect_equal(attr(pheno, "genderCol"), NULL)

  # gender col mis-specified
  expect_error(pheno <- phenotype(pheno1, formula = "y ~ age + pc1 + pc2", 
                               family = "gaussian", id = "id", gender = "sex"), 
            regexp = "is not a column in data")  
    
  p <- pheno1
  p$sex <- sample(c("F", "M"), size = nrow(p), replace = TRUE)
  # gender col specified 
  expect_is(pheno <- phenotype(p, formula = "y ~ age + pc1 + pc2", 
                               family = "gaussian", id = "id", gender = "sex"), 
            "phenotype")
  expect_equal(attr(pheno, "genderCol"), "sex")
  
  # multiple gender cols specified
  expect_error(pheno <- phenotype(p, formula = "y ~ age + pc1 + pc2", 
                                  family = "gaussian", id = "id", gender = c("sex", "y")), 
               regexp = "Only one gender column can be specified")  
  
  # missing gender's
  p[floor(nrow(p)/2), "sex"] <- NA
  expect_error(pheno <- phenotype(p, formula = "y ~ age + pc1 + pc2", 
                                  family = "gaussian", id = "id", gender = "sex"), 
               regexp = "Missing gender data")  
  
})

test_that("include checking", { 
  ids <- c("id3557", "id305") 
  p <- phenotype(pheno1, formula = "y~1", family = "gaussian", id = "id", 
                 include = ids)
  
  expect_equal(get_included(p), ids)
  expect_equal(p[ , get_idCol(p)], ids)
  
  expect_warning(p2 <- phenotype(pheno1, formula = "y~1", family = "gaussian", 
                                 id = "id", 
                                 include = c("id3557", "id305", "SAMPLE1")))
  
  expect_equal(get_included(p2), ids)
  expect_equal(p2[ , get_idCol(p2)], ids)
  
})


test_that("exclude checking", { 
  ids <- c("id3557", "id305") 
  p <- phenotype(pheno1, formula = "y~1", family = "gaussian", id = "id", 
                 exclude = ids)
  
  expect_equal(get_excluded(p), ids)
  expect_equal(any(p[ , get_idCol(p)] %in% ids), FALSE)
  
  expect_warning(p2 <- phenotype(pheno1, formula = "y~1", family = "gaussian", 
                                 id = "id", 
                                 exclude = c("id3557", "id305", "SAMPLE1")))
  
  expect_equal(get_excluded(p2), ids)
  expect_equal(any(p2[ , get_idCol(p2)] %in% ids), FALSE)
  
})

test_that("include-exclude checking", { 
  ids_include <- c("id3557", "id305", "id2977", "id4292") 
  ids_exclude <- c("id305") 
  p <- phenotype(pheno1, formula = "y~1", family = "gaussian", id = "id", 
                 include = ids_include, exclude = ids_exclude)
  
  ids <- setdiff(ids_include, ids_exclude)
  
  expect_equal(sort(get_included(p)), sort(ids))
  expect_equal(sort(get_excluded(p)), sort(setdiff(pheno1$id, ids)))
  
})

  