library(CCP)
library(tidyverse)
library(candisc)


# This code will residualize X and Y datasets for a CCA. Data to residualize by is read in as covariates and an LM residualization is performed. 
# Data is then fit into a prepackaged CCA functions in R with permutation testing. See documentation 


x_data     <- read_csv("path/to/x_data.csv")      # must include subject id + x dataset vars
y_data     <- read_csv("path/to/y_data.csv")      # subject id + y dataset vars
covariates <- read_csv("path/to/covariates_to_residualize_by.csv")  # subject_id + covariates to residualize by 


subject_id_var <- "subject_id" #harmonize naming 
covariate_vars <- c("list", "of", "covs")


residualize_dataframe <- function(data, covariate_data, id_col, cov_cols) {

  # cca variables to residualize
  meas_cols <- setdiff(names(data), id_col)
  
  # initialize residuals data frame with just the IDs
  resid_df <- data %>% select(all_of(id_col))
  
  # loop over cca variables and residualize,
  for (m in meas_cols) {
    tmp <- data %>%
      select(all_of(c(id_col, m))) %>%
      left_join(covariate_data, by = id_col)
    
    form <- reformulate(termlabels = cov_cols, response = m)
    fit  <- lm(form, data = tmp, na.action = na.exclude)
    
    resid_df[[m]] <- resid(fit)
  }
  
  resid_df
}

# Residualize with above function
resid_x <- residualize_dataframe(x_data, covariates,
                                 id_col = subject_id_var,
                                 cov_cols = covariate_vars)

resid_y <- residualize_dataframe(y_data, covariates,
                                 id_col = subject_id_var,
                                 cov_cols = covariate_vars)

# Drop subject_id
X <- resid_x %>% select(-all_of(subject_id_var))
Y <- resid_y %>% select(-all_of(subject_id_var))



### IF YOU DO NOT NEED TO RESIDUALIZE DATA THEN JUST RUN THE CODE BELOW AND MAKE SURE TO HAVE PROPER X AND Y DATASET NAMES 


# Run CCA 
cca_res <- candisc::cancor(x = X, y = Y,
                           row.names = resid_x[[subject_id_var]])

# Permutation of CCA statistics
perm <- p.perm(X, Y, nboot = 1000)
print(perm$stat)


# Below is code to extract individual canonical variate scores and variable loadings. 
# Only 2 here but should be 1 less loading than the total number of vars 

# Get variate scores
scores_df <- bind_cols(
  tibble(!!subject_id_var := resid_x[[subject_id_var]]),
  as_tibble(cca_res$scores$X, .name_repair = "unique"), # fixes duplicate names if any 
  as_tibble(cca_res$scores$Y, .name_repair = "unique"))

# Get Loadings
x_loadings <- tibble(
  Variable  = colnames(X),
  Loading1  = cca_res$structure$X.xscores[, 1],
  Loading2  = cca_res$structure$X.xscores[, 2])

y_loadings <- tibble(
  Variable  = colnames(Y),
  Loading1  = cca_res$structure$Y.yscores[, 1],
  Loading2  = cca_res$structure$Y.yscores[, 2])

