## code to prepare `anes2008_complete_cases` dataset goes here
library(dataverse)

# Get original data file directly from Harvard dataverse for
# Kropko, Jonathan; Goodrich, Ben; Gelman, Andrew; Hill, Jennifer, 2014,
# "Replication data for: Multiple Imputation for Continuous and Categorical Data: Comparing Joint Multivariate Normal and Conditional Approaches",
# https://doi.org/10.7910/DVN/24672, Harvard Dataverse, V3, UNF:5:QuxE8nFhbW2JZT+OW9WzWw== [fileUNF]

anes2008_original <-
  dataverse::get_dataframe_by_name(
    filename    = "base2008_3.tab",
    dataset     = "10.7910/DVN/24672",
    .f          = foreign::read.dta,
    original    = TRUE,
    server      = "dataverse.harvard.edu"
  )

# Apply pre-processing as in Kropko et al. 2014

anes2008_original <-
  anes2008_original[, c(6:8, 11:13, 25, 42, 44:46)]
levels(anes2008_original$imp_enviro) <-
  c("Not important", "Not important", "Important", "Important")
levels(anes2008_original$religion) <-
  c("Protestant",
    "Catholic/Orthodox",
    "Atheist/Other",
    "Atheist/Other")
anes2008_original <-
  anes2008_original[!as.logical(rowSums(is.na(anes2008_original[, c(2, 3, 9, 11)]))),]

anes2008_original[["education"]] <-
  factor(anes2008_original[["education"]], ordered = TRUE)
anes2008_original[["jobs_r"]] <-
  factor(anes2008_original[["jobs_r"]], ordered = TRUE)
anes2008_original[["income"]] <-
  factor(anes2008_original[["income"]], ordered = TRUE)

# Subset to complete cases
anes2008_complete_cases <- na.omit(anes2008_original)

# Store data to package
usethis::use_data(anes2008_complete_cases, overwrite = TRUE)
