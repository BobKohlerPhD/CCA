library(CCP)
library(tidyverse)
library(candisc)


# Set seed for reproducibility
set.seed(123)

n <- 100
cca_covariates <- data.frame(
  age = sample(20:60, n, replace = TRUE),
  sex = sample(c("M", "F"), n, replace = TRUE),
  education = sample(c("HS", "Bachelor", "Master", "PhD"), n, replace = TRUE),
  income = sample(c("Low", "Medium", "High"), n, replace = TRUE))

# The first column is extraneous, the second is subject_id, followed by measurement columns.
x_data_cca <- data.frame(
  extraneous = 1:n,
  subject_id = paste0("subj", 1:n),
  measure1 = rnorm(n, mean = 50, sd = 10),
  measure2 = rnorm(n, mean = 100, sd = 15))

# -------------------------------
# Create fake y_data_cca data frame
# -------------------------------
y_data_cca <- data.frame(
  extraneous = 1:n,
  subject_id = paste0("subj", 1:n),
  measure1 = rnorm(n, mean = 5, sd = 2),
  measure2 = rnorm(n, mean = 10, sd = 3))

# -----------------------------------------------------------
# Define the generalized residualization function using GLMs
# -----------------------------------------------------------
residualize_dataframe <- function(data, covariate_formula, covariate_data, subject_id_col = "subject_id", family = gaussian()) {
 
   # Drop index or unwanted variable
  df <- data[, -1]  
  
  # Get subject IDs
  subject_ids <- df[[subject_id_col]]
  
  measurement_names <- setdiff(names(df), subject_id_col)
  
  # Residualize each measurement using GLM
  residuals_list <- lapply(measurement_names, function(colname) {
    # Extract the measurement vector from the data frame
    measurement_vector <- df[[colname]]
    
    # Combine the measurement vector with the covariate data
    temp_data <- cbind(covariate_data, measurement_vector)
    # Rename the new column to the measurement's name so that the formula finds it
    colnames(temp_data)[ncol(temp_data)] <- colname
    
    # Build the model formula for the current measurement
    formula_str <- paste(colname, deparse(covariate_formula))
    fit_formula <- as.formula(formula_str)
    
    # Fit the model and extract residuals
    fit <- glm(fit_formula, family = family, data = temp_data, na.action = na.exclude)
    residuals(fit)
  })
  
  # Combine residuals into a data frame and re-attach the subject IDs
  residuals_df <- as.data.frame(residuals_list)
  names(residuals_df) <- measurement_names
  residuals_df <- cbind(subject_id = subject_ids, residuals_df)
  return(residuals_df)
}


# Create model with variables to residualize
cov_model <- as.formula(" ~ age + factor(sex) + factor(education) + factor(income)")

# Residualize X & Y dataset
residuals_x <- residualize_dataframe(x_data_cca, 
                                     cov_model, 
                                     cca_covariates, 
                                     subject_id_col = "subject_id", 
                                     family = gaussian())

residuals_y <- residualize_dataframe(y_data_cca, 
                                     cov_model, 
                                     cca_covariates,
                                     subject_id_col = "subject_id", 
                                     family = gaussian())

subject_id <- x_data_cca$subject_id

# Remove subject id and run CCA + permutation testing
residual_cca_x <- residuals_x[ ,-1]
residual_cca_y <- residuals_y[ ,-1]
cca_result <- candisc::cancor(x = residual_cca_x,
                              y = residual_cca_y,
                              row.names=subject_id)

cca_permutation_test <- p.perm(residual_cca_x, 
                               residual_cca_y, 
                               nboot = 100, 
                               rhostart = 1)
print(cca_permutation_test$stat)

# Extract variate scores for each participant 
cca_scores <- data.frame(src_subject_id = subject_id, 
                         cca_result$scores$X, 
                         cca_result$scores$Y)

# Extract variable loadings for the first 3 modes
x_var_names <- colnames(residual_cca_x)
x_cca <- data.frame(Variable = x_var_names, 
                    Loadings_X_Mode1 = cca_result$structure$X.xscores[, 1], 
                    Loadings_X_Mode2 = cca_result$structure$X.xscores[, 2])

# If a dataframe is needed for each mode
x_mode1_df <- x_cca %>% select(Variable, Loadings_X_Mode1)
x_mode2_df <- x_cca %>% select(Variable, Loadings_X_Mode2)

y_var_names <- colnames(residual_cca_y)
y_cca <- data.frame(Variable = y_var_names, 
                    Loadings_Y_Mode1 = cca_result$structure$Y.yscores[, 1], 
                    Loadings_Y_Mode2 = cca_result$structure$Y.yscores[, 2])

# If a dataframe is needed for each mode
y_mode1_df <- y_cca %>% select(Variable, Loadings_Y_Mode1)
y_mode2_df <- y_cca %>% select(Variable, Loadings_Y_Mode2)

# Descending Order
x_mode1_df <- x_mode1_df  %>%
  arrange(desc(abs(Loadings_X_Mode1))) %>%
  mutate(Variable = reorder(Variable, Loadings_X_Mode1))
# Plot first 15 variable loadings for X dataset
ggplot(x_mode1_df, aes(x= Variable,y=Loadings_X_Mode1, fill = Loadings_X_Mode1)) +
  geom_bar(position = "dodge", stat="identity",width =.85, color = "black", linewidth = 1) +
  scale_fill_gradient2(low = "darkblue",high = "dodgerblue", midpoint = 0) +
  coord_flip()+
  labs(x = "", y = "Variable Loading", title = "Mode 1 X")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.line = element_line(linewidth = 1.25), 
        axis.ticks.length = unit(.25,"cm"),
        axis.text.y = element_text(size = 16,  color = "black"),
        axis.text.x = element_text(size = 14,  color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16, color = "black", face = "bold"),
        legend.key.size  = unit(1.5,"cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.position = "none")

# Descending Order
y_mode1_df <- y_mode1_df %>%
  arrange(desc(abs(Loadings_Y_Mode1))) %>%
  mutate(Variable = reorder(Variable, Loadings_Y_Mode1))
# Plot first 15 variable loadings for Y dataset 
ggplot(y_mode1_df, aes(x= Variable,y=Loadings_Y_Mode1, fill = Loadings_Y_Mode1)) +
  geom_bar(position = "dodge", stat="identity",width =.85, color = "black", linewidth = 1) +
  scale_fill_gradient2(low = "darkred",high = "salmon1", midpoint = 0) +
  coord_flip(ylim = c(-0.35, 0.00), expand = T)+
  labs(x = "", y = "Variable Loading", , title = "Mode 1 Y")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.line = element_line(linewidth = 1.25), 
        axis.ticks.length = unit(.25,"cm"),
        axis.text.y = element_text(size = 16,  color = "black"),
        axis.text.x = element_text(size = 14,  color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16, color = "black", face = "bold"),
        legend.key.size  = unit(1.5,"cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.position = "none")
