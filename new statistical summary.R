library(readr)
library(dplyr)
library(tidyr)

######### CONFLICT EXPOSED PROTECTED ###############

# Reshape, drop NAs, and add the new variables
df_longConfProt <- Africa_LPI_Conflict_exposed___Protected_pops_PopLambda %>%
  pivot_longer(
    cols = matches("^X[0-9]{4}$"),  
    names_to = "year",
    values_to = "value"
  ) %>%
  mutate(
    year = as.integer(sub("^X", "", year)),
    Conflict = 1,
    Climate = 0,
    Protected = 1
  ) %>%
  drop_na(value)


################## CONFLICT EXPOSED UNPROTECTED ##############

# Reshape, drop NAs, and add the new variables
df_longConfUnprot <- Africa_LPI_Conflict_exposed___Unprotected_pops_PopLambda %>%
  pivot_longer(
    cols = matches("^X[0-9]{4}$"), 
    names_to = "year",
    values_to = "value"
  ) %>%
  mutate(
    year = as.integer(sub("^X", "", year)), 
    Conflict = 1,
    Climate = 0,
    Protected = 0
  ) %>%
  drop_na(value)


# Combine with the previous dataframe
df_longCombined <- bind_rows(df_longConfProt, df_longConfUnprot)


################## CONFLICT EXPOSED CLIMATE EXPOSED PROTECTED ##############

# Reshape, drop NAs, and add the new variables
df_longConfClimProt <- Africa_LPI_Climate_and_conflict_exposed___Protected_pops_PopLambda %>%
  pivot_longer(
    cols = matches("^X[0-9]{4}$"),  
    names_to = "year",
    values_to = "value"
  ) %>%
  mutate(
    year = as.integer(sub("^X", "", year)),  
    Conflict = 1,
    Climate = 1,
    Protected = 1
  ) %>%
  drop_na(value)


# Combine with the previous dataframe
df_longCombined <- bind_rows(df_longCombined, df_longConfClimProt)


################## CONFLICT EXPOSED CLIMATE EXPOSED UNPROTECTED ##############

# Reshape, drop NAs, and add the new variables
df_longConfClimUnprot <- Africa_LPI_Climate_and_Conflict_exposed___Unprotected_pops_PopLambda %>%
  pivot_longer(
    cols = matches("^X[0-9]{4}$"),  # match X1971, X1972, etc.
    names_to = "year",
    values_to = "value"
  ) %>%
  mutate(
    year = as.integer(sub("^X", "", year)),  # remove the X before converting to number
    Conflict = 1,
    Climate = 1,
    Protected = 0
  ) %>%
  drop_na(value)


################## NO CONFLICT POPULATIONS ##############
df_longNoConf <- Africa_LPI_no_conflict_pops_PopLambda %>%
  pivot_longer(
    cols = matches("^X[0-9]{4}$"),
    names_to = "year",
    values_to = "value"
  ) %>%
  mutate(
    year = as.integer(sub("^X", "", year)),
    Conflict = 0,
    Climate = 0,
    Protected = 0  
  ) %>%
  drop_na(value)

###### CLIMATE ONLY POPULATIONS ########
df_longClimate <- Africa_LPI_Climate_only_pops_PopLambda %>%
  pivot_longer(
    cols = matches("^X[0-9]{4}$"),
    names_to = "year",
    values_to = "value"
  ) %>%
  mutate(
    year = as.integer(sub("^X", "", year)),
    Conflict = 0,
    Climate = 0,
    Protected = 0  
  ) %>%
  drop_na(value)

#Combining 
# Add ExposureType
df_longConfProt <- df_longConfProt %>%
  mutate(ExposureType = "ConflictOnly_Prot")

df_longConfUnprot <- df_longConfUnprot %>%
  mutate(ExposureType = "ConflictOnly_Unprot")

df_longConfClimProt <- df_longConfClimProt %>%
  mutate(ExposureType = "ConflictClimate_Prot")

df_longConfClimUnprot <- df_longConfClimUnprot %>%
  mutate(ExposureType = "ConflictClimate_Unprot")

df_longNoConf <- df_longNoConf %>%
  mutate(ExposureType = "NoConflict")

df_longClimate <- df_longClimate %>%
  mutate(ExposureType = "ClimateOnly")

df_longClimate <- df_longClimate %>%
  mutate(
    Conflict = 0,
    Climate = 1,   # set Climate = 1 for ClimateOnly
    Protected = 0
  )

#Combine all datasets
df_longCombined <- bind_rows(df_longConfProt, df_longConfUnprot,
                             df_longConfClimProt, df_longConfClimUnprot,
                             df_longNoConf, df_longClimate)

df_longCombined_dt <- as.data.table(df_longCombined)
df_longCombined$ExposureType <- factor(df_longCombined$ExposureType,
                                       levels = c("NoConflict", 
                                                  "ConflictOnly_Unprot", 
                                                  "ConflictOnly_Prot", 
                                                  "ConflictClimate_Unprot", 
                                                  "ConflictClimate_Prot", 
                                                  "ClimateOnly"))

#Convert to data.table and join species/site info
df_longCombined_dt <- as.data.table(df_longCombined)
pop_final_dt_unique <- unique(as.data.table(pop_final[, c("ID", "Binomial", "site")]), by = "ID")

# Join population info
df_longCombined_dt <- pop_final_dt_unique[df_longCombined_dt, on = c("ID" = "population_id")]

#Relevel factor so NoConflict is reference
df_longCombined_dt$ExposureType <- relevel(df_longCombined_dt$ExposureType, ref = "NoConflict")

#Fit Linear Mixed Model
lmm <- lmer(value ~ year + Protected + Climate*Conflict  + (1 | site) + (1 | Binomial),
            data = df_longCombined_dt)


# Libraries
library(broom.mixed)
library(dplyr)
library(kableExtra)

# Tidy up fixed effects
fixed_tidy <- tidy(lmm, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    Estimate = round(estimate, 3),
    `Std. Error` = round(std.error, 3),
    `95% CI` = paste0("(", round(conf.low, 3), ", ", round(conf.high, 3), ")"),
    `p-value` = signif(p.value, 3)
  ) %>%
  select(term, Estimate, `Std. Error`, `95% CI`, `p-value`) %>%
  mutate(effect = "Fixed")

# Rename predictors for readability
fixed_tidy$term <- recode(fixed_tidy$term,
                          "(Intercept)" = "Intercept",
                          "year" = "Year",
                          "Protected" = "Protected",
                          "Climate" = "Climate",
                          "Conflict" = "Conflict",
                          "Climate:Conflict" = "Climate × Conflict")

# Random effects variances
rand_var <- as.data.frame(VarCorr(lmm)) %>%
  select(grp, vcov) %>%
  rename(Term = grp, Variance = vcov) %>%
  mutate(`Std. Dev.` = round(sqrt(Variance), 3),
         effect = "Random",
         Estimate = NA, `Std. Error` = NA, `95% CI` = NA, `p-value` = NA)

# Combine fixed + random effects
full_table <- bind_rows(
  fixed_tidy %>% rename(Term = term),
  rand_var
)

# Print final table
full_table %>%
  kable(
    caption = "Linear Mixed Model: Effects of Conflict, Climate, and Protection on Population Values",
    booktabs = TRUE,
    align = "lcccccc"
  ) %>%
  kable_styling(full_width = FALSE, position = "center",
                bootstrap_options = c("striped", "hover"))



# Load necessary libraries
library(lme4)
library(ggeffects)
library(ggplot2)



# 1. Create predicted values for Protected = 0
pred_protected0 <- ggpredict(lmm, 
                             terms = c("year", "Climate", "Conflict"), 
                             condition = c(Protected = 0))

# Plot for Protected = 0
plot_protected0 <- ggplot(pred_protected0, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(x = "Year", y = "Predicted Value", color = "Climate-Conflict", fill = "Climate-Conflict",
       title = "Predictions for Protected = 0") +
  theme_minimal()

# 2. Create predicted values for Protected = 1
pred_protected1 <- ggpredict(lmm, 
                             terms = c("year", "Climate", "Conflict"), 
                             condition = c(Protected = 1))

# Plot for Protected = 1
plot_protected1 <- ggplot(pred_protected1, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(x = "Year", y = "Predicted Value", color = "Climate-Conflict", fill = "Climate-Conflict",
       title = "Predictions for Protected = 1") +
  theme_minimal()

# Display plots
plot_protected0
plot_protected1

















