#Additional results 

pop_final <- pop_final %>%
  mutate(
    Protected_status = case_when(
      tolower(Protected_status) == "yes" ~ "Protected",
      tolower(Protected_status) == "no" ~ "Unprotected",
      TRUE ~ NA_character_
    )
  )


#Conflict coverage by class
coverage_by_class <- pop_final%>%
  group_by(Class, conflict_exposed) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Class) %>%
  mutate(percentage = round(100 * n / sum(n), 1))

class_totals <- coverage_by_class %>%
  group_by(Class) %>%
  summarise(total_n = sum(n), .groups = "drop")
coverage_by_class <- left_join(coverage_by_class, class_totals, by = "Class") %>%
  mutate(Class_label = paste0(Class, " (n=", total_n, ")"))

ggplot(coverage_by_class, aes(x = Class_label, y = percentage, fill = conflict_exposed)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +  # black borders
  scale_fill_manual(values = c("#1b9e77", "#d95f02")) +  # custom colors for levels of conflict_exposed
  geom_text(aes(label = paste0(percentage, "%")),
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3.5) +
  labs(
    title = "Class Coverage by Conflict Exposure",
    x = "Class", 
    y = "Percentage of Populations"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#Monitoring gaps 
#Monitoring duration per population 
monitoring_duration_df <- pop_final %>%
  filter(!is.na(popvalue)) %>%  # consider only years with data
  group_by(ID, conflict_exposed.x) %>%
  summarise(
    first_year = min(year, na.rm = TRUE),
    last_year = max(year, na.rm = TRUE),
    duration = last_year - first_year + 1,
    n_years = n(),  
    .groups = "drop"
  )
#Summary of monitoring duration y conflict exposure
monitoring_duration_df %>%
  group_by(conflict_exposed.x) %>%
  summarise(
    mean_duration = mean(duration),
    median_duration = median(duration),
    mean_years_monitored = mean(n_years),
    median_years_monitored = median(n_years),
    n_populations = n()
  )
ggplot(monitoring_duration_df, aes(x = factor(conflict_exposed.x), y = duration)) +
  geom_boxplot(fill = c("#56B4E9", "#E69F00")) +
  labs(
    x = "Conflict Exposed",
    y = "Monitoring Duration (years)",
    title = "Estimated Monitoring Duration by Conflict Exposure"
  ) +
  scale_x_discrete(labels = c("No", "Yes")) +
  theme_minimal()


#Taxonomic class monitoring gaps 
pops_sf_test <- pops_sf %>%
  mutate(
    conflict_exposed.x = ifelse(year < 1997, FALSE, conflict_exposed)
  )

taxonomic_gaps_all <- pops_sf_test %>%
  st_drop_geometry() %>%
  filter(!is.na(popvalue)) %>%
  filter(Class != "Plantae") %>%
  mutate(
    Class = case_when(
      Class %in% c("Dipneusti", "Actinopteri") ~ "Fishes",
      TRUE ~ Class
    ),
    exposure_protection = case_when(
      conflict_exposed.x == TRUE & protection_status == "Protected" ~ "Conflict exposed - Protected",
      conflict_exposed.x == TRUE & protection_status == "Unprotected" ~ "Conflict exposed - Unprotected",
      conflict_exposed.x == FALSE & protection_status == "Protected" ~ "No exposure - Protected",
      conflict_exposed.x == FALSE & protection_status == "Unprotected" ~ "No exposure - Unprotected",
      TRUE ~ "Unknown protection status"
    ),
    exposure_protection = factor(
      exposure_protection,
      levels = c(
        "Conflict exposed - Protected",
        "Conflict exposed - Unprotected",
        "No exposure - Protected",
        "No exposure - Unprotected",
        "Unknown protection status"
      )
    )
  ) %>%
  group_by(exposure_protection, Class) %>%
  summarise(
    species_monitored = n_distinct(species),
    populations_monitored = n(),
    .groups = "drop"
  )

ggplot(taxonomic_gaps_all, 
       aes(x = reorder(Class, populations_monitored), 
           y = populations_monitored, 
           fill = exposure_protection)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = populations_monitored),
            position = position_dodge(width = 0.9),
            hjust = -0.2,
            size = 3) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Conflict exposed - Protected" = "chartreuse3",
    "Conflict exposed - Unprotected" = "red3", 
    "No exposure - Protected" = "cyan3",
    "No exposure - Unprotected" = "chocolate1",
    "Unknown protection status" = "grey50"
  )) +
  labs(
    title = "Taxonomic Monitoring Gaps by Exposure and Protection Group - Africa",
    y = "Number of Populations Monitored",
    x = "Class",
    fill = "Exposure + Protection Status"
  ) +
  theme_minimal()



# LPI by taxonomic group 
lpi_by_class_results <- list()

for (cls in classes) {
  message("Processing class: ", cls)
  
  pop_class <- pop_final %>%
    filter(Class == cls)
  
  pop_class_wide <- pop_class %>%
    select(ID, Binomial = species, Country, year, popvalue) %>%
    mutate(
      Binomial = str_trim(tolower(Binomial)),
      year = paste0("X", year)
    ) %>%
    pivot_wider(names_from = year, values_from = popvalue)
  
  class_name_clean <- str_replace_all(cls, "[^A-Za-z0-9]", "_")
  infile_name <- paste0("LPI_class_", class_name_clean)
  
  create_infile(
    pop_data_source = pop_class_wide,
    name = infile_name,
    start_col_name = "X1950",
    end_col_name = "X2022"
  )
  result <- LPIMain(paste0(infile_name, "_infile.txt"))
  
  result_clean <- result %>%
    mutate(
      Year = as.numeric(row.names(.)),
      Class = cls
    ) %>%
    filter(!is.na(Year) & is.finite(LPI_final))
  
  lpi_by_class_results[[cls]] <- result_clean
}
lpi_by_class_df <- bind_rows(lpi_by_class_results)
ggplot(lpi_by_class_df, aes(x = Year, y = LPI_final, color = Class)) +
  geom_line(size = 1) +
  labs(
    title = "LPI Curves by Taxonomic Class",
    x = "Year",
    y = "Living Planet Index"
  ) +
  theme_minimal() +
  theme(legend.position = "right")





#Role of protection status 
#1 data source : WDPA
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
pops_sf <- pop_final %>%
  filter(!is.na(Longitude) & !is.na(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
# 3. Spatial join: check which population points fall inside protected areas
pops_sf <- st_join(pops_sf, protected_areas, join = st_within)

pops_sf <- pops_sf %>%
  mutate(
    protection_status = ifelse(!is.na(NAME), "Protected", "Unprotected")
  )
protection_summary <- pops_sf %>%
  group_by(protection_status, conflict_exposed.x) %>%
  summarise(n_populations = n(), .groups = "drop") %>%
  group_by(protection_status) %>%
  mutate(percentage_within_group = round(100 * n_populations / sum(n_populations), 1))
#Visualize the coverage
ggplot(protection_summary, aes(x = protection_status, y = percentage_within_group, fill = conflict_exposed.x)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(percentage_within_group, "%")),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4) +
  labs(title = "Proportion of Conflict Exposure\nWithin Protected vs Unprotected Areas",
       x = "Protection Status",
       y = "Percentage of LPI Populations",
       fill = "Conflict Exposure") +
  scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "#1a9850")) +
  theme_minimal()

#2 data source : LPI
protection_summary <- pops_sf %>%
  group_by(Protected_status, conflict_exposed.x) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Protected_status) %>%
  mutate(total = sum(count),
         percentage_within_group = round((count / total) * 100, 1)) %>%
  ungroup()

#Standarization of values into three categories : protected, unprotected, unknown status
pops_sf <- pops_sf %>%
  mutate(protection_status = case_when(
    Protected_status %in% c("Yes", "Protected", "Both") ~ "Protected",
    Protected_status %in% c("No", "Unprotected", "No (area surrounding the PA)",
                             "No (large survey area)") ~ "Unprotected",
    TRUE ~ "Unknown"
  ))
protection_summary <- pops_sf %>%
  group_by(protection_status, conflict_exposed.x) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(protection_status) %>%
  mutate(total = sum(count),
         percentage_within_group = round((count / total) * 100, 1)) %>%
  ungroup()
ggplot(protection_summary, aes(x = protection_status,
                               y = percentage_within_group,
                               fill = as.factor(conflict_exposed.x))) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(percentage_within_group, "%")),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4) +
  labs(title = "Proportion of Conflict Exposure\nWithin Protected vs Unprotected Areas",
       x = "Protection Status",
       y = "Percentage of LPI Populations",
       fill = "Conflict Exposure") +
  scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "#1a9850")) +
  theme_minimal()


# Map of africa 
ggplot(pops_sf) +
  geom_sf(data = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf"), 
          fill = "gray95", color = "gray80") +
  geom_sf(aes(color = protection_flag, shape = protection_flag), 
          alpha = 0.7, size = 2) +
  scale_color_manual(
    name = "Protection & Conflict Status",
    values = c(
      "Protected & Conflict" = "deepskyblue",
      "Protected & No Conflict" = "chartreuse3",
      "Unprotected & Conflict" = "brown2",
      "Unprotected & No Conflict" = "darkgoldenrod1",
      "Unkown" = "grey"
    )
  ) +
  scale_shape_manual(
    name = "Protection & Conflict Status",
    values = c(
      "Protected & Conflict" = 16,     # dot
      "Protected & No Conflict" = 17,  # triangle
      "Unprotected & Conflict" = 16,   # dot
      "Unprotected & No Conflict" = 17# triangle
    )
  ) +
  coord_sf(xlim = c(-20, 55), ylim = c(-40, 40)) +
  labs(
    title = "Map of African populations according to their status",
    subtitle = "Protection and Conflict Exposure",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )



#LPI curve by expsoure and protection status
pops_sf <- pops_sf %>%
  mutate(
    protection_flag = case_when(
      protection_status == "Protected" & exposure_group == "Conflict only" ~ "Conflict exposed - Protected",
      protection_status == "Protected" & exposure_group == "Both conflict + climate" ~ "Climate and conflict exposed - Protected",
      protection_status == "Unprotected" & exposure_group == "Conflict only"  ~ "Conflict exposed - Unprotected",
      protection_status == "Unprotected" & exposure_group == "Both conflict + climate" ~ "Climate and Conflict exposed - Unprotected",
      TRUE ~ "Unknown"
      )
 ) 
run_lpi_by_protection_flag <- function(data, flag_col, flag_levels) {
  
  map_dfr(flag_levels, function(flag) {
    message("Processing: ", flag)
    
    # Subset and drop geometry if present
    pop_subset <- data %>%
      filter(!!sym(flag_col) == flag) %>%
      st_drop_geometry() %>%  # <-- DROP GEOMETRY COLUMN
      filter(!Class %in% "Plantae") %>%
      mutate(Class = case_when(
        Class %in% c("Dipneusti", "Actinopteri") ~ "Fishes",
        TRUE ~ Class
      ))
    
    if (nrow(pop_subset) == 0) {
      warning("No data for group: ", flag)
      return(NULL)
    }
    
    # Pivot to wide format
    pop_wide <- pop_subset %>%
      select(ID, Binomial = species, Country, year, popvalue) %>%
      mutate(
        Binomial = str_trim(tolower(Binomial)),
        year = paste0("X", year)
      ) %>%
      pivot_wider(names_from = year, values_from = popvalue)
    
    # Identify year columns
    year_cols <- names(pop_wide)[grepl("^X\\d{4}$", names(pop_wide))]
    if (length(year_cols) == 0) {
      warning("No year columns found for group: ", flag)
      return(NULL)
    }
    
    start_col <- min(year_cols)
    end_col <- max(year_cols)
    
    infile_name <- paste0("Africa_LPI_", str_replace_all(flag, "[^A-Za-z0-9]", "_"))
    create_infile(
      pop_data_source = pop_wide,
      name = infile_name,
      start_col_name = start_col,
      end_col_name = end_col
    )
    
    result <- LPIMain(paste0(infile_name, "_infile.txt"), LEV_FLAG = TRUE)

    result %>%
      mutate(
        Year = as.numeric(row.names(.)),
        Group = flag
      ) %>%
      filter(!is.na(Year), is.finite(LPI_final))
  })
}
groups_to_run <- c(
  "Conflict exposed - Protected",
  "Conflict exposed - Unprotected",
  "Climate and conflict exposed - Protected",
  "Climate and Conflict exposed - Unprotected" # Note: Check your capitalization!
)

lpi_protection_curves <- run_lpi_by_protection_flag(pops_sf, "protection_flag", groups_to_run)
lpi_plot_data <- lpi_protection_curves %>%
  filter(
    (Group == "Climate and conflict exposed - Protected" & Year <= 2016) |
      (Group != "Climate and conflict exposed - Protected")
  )

percent_changes <- lpi_plot_data %>%
  group_by(Group) %>%
  summarise(
    start_year = min(Year),
    end_year = max(Year),
    start = LPI_final[Year == start_year][1],
    end = LPI_final[Year == end_year][1],
    percent_change = round((end - start) / start * 100, 1)
  ) %>%
  mutate(label = paste0(percent_change, "%"))

lpi_plot_data <- lpi_plot_data %>%
  left_join(percent_changes %>% select(Group, label), by = "Group")
ggplot(lpi_plot_data, aes(x = Year, y = LPI_final, color = Group)) +
  geom_line(size = 1) +
  labs(
    title = "Living Planet Index by Protection and Exposure Status - Africa",
    y = "LPI (relative abundance)",
    x = "Year",
    color = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  ) +
  geom_text(
    data = lpi_plot_data %>% group_by(Group) %>% filter(Year == max(Year)),
    aes(label = label),
    hjust = 0, vjust = -0.5, size = 4, show.legend = FALSE
  )











