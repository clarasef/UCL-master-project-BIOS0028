#Filter data for Kenya
case_kenya <-pop_final %>%
  filter(Country == "Kenya")

glimpse(case_kenya)
summary(case_kenya)
n_distinct(case_kenya$species) 
case_kenya %>% count(species, sort = TRUE)
table(case_kenya$year)
table(case_kenya$category)



#Map of Kenya proteted vs unprotected and conflicts 
library(rnaturalearth)
library(rnaturalearthdata)
case_kenya <- case_kenya %>%
  mutate(
    Protected_status = case_when(
      tolower(Protected_status) == "yes" ~ "Protected",
      tolower(Protected_status) == "no" ~ "Unprotected",
      TRUE ~ NA_character_
    )
  )
case_kenya_sf <- case_kenya %>%
  filter(year >= 1997,!is.na(Protected_status), !is.na(Latitude), !is.na(Longitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

#Map 
kenya_boundary <- ne_countries(scale = "medium", country = "Kenya", returnclass = "sf")
data(world.cities)
kenya_cities <- world.cities %>%
  filter(country.etc == "Kenya", pop > 50000)  
kenya_cities_sf <- st_as_sf(kenya_cities, coords = c("long", "lat"), crs = 4326)

#Conflict in Kenya 
kenya_conflict <- acled_africa %>%
  filter(year >= 2020,year <= 2025, country == "Kenya") %>% #different years according to time-serie
  filter(!is.na(latitude), !is.na(longitude), event_type == "Battles") %>%
  mutate(
    date = ymd(event_date),
    year = year(date),
    intensity = replace_na(fatalities, 0),
    event_id = row_number()
  ) %>%
  select(year, latitude, longitude, intensity, event_id)

kenya_conflict_sf <- kenya_conflict %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

#Protected areas from WDPA
protected_areas <- st_read("WDPA_WDOECM_Jun2025_Public_AF_shp-polygons.shp")
protected_areas_kenya <- protected_areas %>%
  filter(ISO3 == "KEN")

conflicts_in_protected <- st_join(kenya_conflict_sf, protected_areas_kenya, join = st_within)
n_conflicts_protected <- nrow(conflicts_in_protected %>% filter(!is.na(NAME)))
conflicts_by_area <- conflicts_in_protected %>%
  filter(!is.na(NAME)) %>%
  count(NAME, STATUS_YR, sort = TRUE)

#Plotting Map 
gggplot() +
  geom_sf(data = kenya_boundary, fill = "gray95", color = "black") +
  geom_sf(data = protected_areas_kenya, fill = "palegreen3", color = "darkgreen", alpha = 0.4) +
  geom_sf(data = case_kenya_sf, aes(color = Protected_status), size = 1.5, alpha = 0.8) + 
  geom_sf(data = kenya_cities_sf, aes(fill = "Cities"), shape = 21, size = 1.5, alpha = 0.7) +
  geom_sf(data = kenya_conflict_sf, aes(shape = "Conflict"), color = "darkorange", size = 2, alpha = 0.9) +
  geom_text(
    data = kenya_cities_sf,
    aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], label = name),
    size = 2.5, check_overlap = TRUE, fontface = "italic", nudge_y = 0.12
  ) +
  scale_color_manual(values = c("Protected" = "cyan3", "Unprotected" = "red")) +
  scale_fill_manual(values = c("Cities" = "blue"), name = "") +
  scale_shape_manual(values = c("Conflict" = 4), name = "") + 
  labs(
    title = "Population Sites in Kenya by Protection Status, Conflict Areas & Cities (2011-2019)",
    color = "Protection Status"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )




#Calculation of the LPI (conflict only, both climate and conflict, total population)
kenya_grouped_lpi <- function(data, group_name) {
  # Subset data if not all populations
  if (group_name != "All populations") {
    data <- data %>% filter(exposure_group == group_name)
  }
  pop_subset <- data %>%
    select(ID, Binomial = species, Country, year, popvalue) %>%
    mutate(
      Binomial = str_trim(tolower(Binomial)),
      year = paste0("X", year)
    ) %>%
    pivot_wider(names_from = year, values_from = popvalue)
  year_cols <- names(pop_subset)[grepl("^X\\d{4}$", names(pop_subset))]
  if (length(year_cols) == 0) {
    message(paste("No valid year columns for group:", group_name))
    return(NULL)
  }
  start_col <- min(year_cols)
  end_col <- max(year_cols)
  safe_group_name <- str_replace_all(group_name, "[^A-Za-z0-9]+", "_")
  infile_name <- paste0("Kenya_LPI_", safe_group_name)
  create_infile(
    pop_data_source = pop_subset,
    name = infile_name,
    start_col_name = start_col,
    end_col_name = end_col
  )
  result <- LPIMain(paste0(infile_name, "_infile.txt"))
  result_clean <- result %>%
    mutate(Year = as.numeric(row.names(.))) %>%
    filter(!is.na(Year) & is.finite(LPI_final)) %>%
    mutate(Group = group_name)
  
  return(result_clean)
}

kenya_groups <- c("Conflict only", "Both conflict + climate", "All populations")
lpi_kenya <- bind_rows(lapply(kenya_groups, function(g) kenya_grouped_lpi(case_kenya, g)))

ggplot(lpi_kenya, aes(x = Year, y = LPI_final, color = Group)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "Living Planet Index in Kenya by Exposure Group",
    x = "Year",
    y = "LPI (Relative Abundance)",
    color = "Group"
  ) +
  scale_color_manual(values = c(
    "No exposure" = "chartreuse3",
    "Conflict only" = "cyan3",
    "Both conflict + climate" = "darkred",
    "All populations" = "chocolate1"
  )) +
  theme_minimal(base_size = 14)



# LPI Kenya (conflict/protected vs unprotected, climate and conflict/protected vs unprotecetd)
kenya_sf <- pops_sf %>%
  filter(Country == "Kenya") %>%
  mutate(
    protection_flag = case_when(
      protection_status == "Protected" & exposure_group == "Conflict only" ~ "Conflict exposed - Protected",
      protection_status == "Protected" & exposure_group == "Both conflict + climate" ~ "Climate and conflict exposed - Protected",
      protection_status == "Unprotected" & exposure_group == "Conflict only" ~ "Conflict exposed - Unprotected",
      protection_status == "Unprotected" & exposure_group == "Both conflict + climate" ~ "Climate and Conflict exposed - Unprotected",
      TRUE ~ "Unknown"
    )
  )

#Calculation
run_lpi_by_protection_flag <- function(data, flag_col, flag_levels) {
  map_dfr(flag_levels, function(flag) {
    message("Processing: ", flag)
    
    pop_subset <- data %>%
      filter(!!sym(flag_col) == flag) %>%
      st_drop_geometry() %>%
      filter(!Class %in% "Plantae") %>%
      mutate(Class = case_when(
        Class %in% c("Dipneusti", "Actinopteri") ~ "Fishes",
        TRUE ~ Class
      ))
    
    if (flag == "Conflict exposed - Protected") {
      print(head(pop_subset))
      print(summary(pop_subset))
      # Optionally save to csv for offline inspection
      write.csv(pop_subset, paste0("pop_subset_", flag, ".csv"), row.names = FALSE)
    }
    
    if (nrow(pop_subset) == 0) {
      warning("No data for group: ", flag)
      return(NULL)
    }
    
    pop_wide <- pop_subset %>%
      select(ID, Binomial = species, Country, year, popvalue) %>%
      mutate(
        Binomial = str_trim(tolower(Binomial)),
        year = paste0("X", year)
      ) %>%
      pivot_wider(names_from = year, values_from = popvalue)
    
    year_cols <- names(pop_wide)[grepl("^X\\d{4}$", names(pop_wide))]
    if (length(year_cols) == 0) {
      warning("No year columns found for group: ", flag)
      return(NULL)
    }
    
    start_col <- min(year_cols)
    end_col <- max(year_cols)
    
    infile_name <- paste0("Kenya_LPI_", str_replace_all(flag, "[^A-Za-z0-9]", "_"))
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
  "Climate and Conflict exposed - Unprotected"
)

lpi_protection_curves <- run_lpi_by_protection_flag(kenya_sf, "protection_flag", groups_to_run)
lpi_plot_data <- lpi_protection_curves %>%
  filter((Group == "Climate and conflict exposed - Protected" & Year <= 2016) |
           (Group != "Climate and conflict exposed - Protected"))

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
  left_join(percent_changes %>% select(Group, label), by = "Group") %>%
  filter(Group != "Conflict exposed - Protected")

ggplot(lpi_plot_data, aes(x = Year, y = LPI_final, color = Group)) +
  geom_line(size = 1) +
  labs(
    title = "Living Planet Index by Protection and Exposure Status - Kenya",
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



#Coverage by class group 
coverage_by_class_kenya <- case_kenya %>%
  filter(Class != "Plantae") %>% 
  group_by(Class, conflict_exposed.x) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Class) %>%
  mutate(percentage = round(100 * n / sum(n), 1))
ggplot(coverage_by_class_kenya, aes(x = Class, y = percentage, fill = conflict_exposed.x)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Class Coverage by Conflict Exposure in Kenya",
       x = "Class", y = "Percentage of Populations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#By ptotection status
protection_summary_kenya <- pops_sf %>%
  filter(Country == "Kenya") %>%
  group_by(protection_status, conflict_exposed.x) %>%
  summarise(n_populations = n(), .groups = "drop") %>%
  group_by(protection_status) %>%
  mutate(percentage_within_group = round(100 * n_populations / sum(n_populations), 1))
  
#Visualize the coverage
ggplot(protection_summary_kenya, aes(x = protection_status, y = percentage_within_group, fill = conflict_exposed.x)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(percentage_within_group, "%")),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4) +
  labs(title = "Proportion of Conflict Exposure\nWithin Protected vs Unprotected Areas in Kenya",
       x = "Protection Status",
       y = "Percentage of LPI Populations",
       fill = "Conflict Exposure") +
  scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "#1a9850")) +
  theme_minimal()


#Monitoring gaps in kenya by taxonomic class
taxonomic_gaps <- kenya_sf %>%
  st_drop_geometry() %>%
  filter(!is.na(popvalue)) %>%
  filter(Class!= "Plantae") %>%
  filter(Class!= "Actinopteri") %>%
  group_by(protection_flag, Class) %>%
  summarise(
    species_monitored = n_distinct(species),
    populations_monitored = n(),
    .groups = "drop"
  )

ggplot(taxonomic_gaps, aes(x = Class, y = populations_monitored, fill = protection_flag)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = populations_monitored),
            position = position_dodge(width = 0.9),
            hjust = -0.2, # moves text outside the bar, can adjust
            size = 3) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Conflict exposed - Protected" = "chartreuse3",
    "Climate and conflict exposed - Protected" = "cyan3",
    "Conflict exposed - Unprotected" = "chocolate1",
    "Climate and Conflict exposed - Unprotected" = "red3" ))+
  labs(
    title = "Taxonomic Monitoring Gaps by Exposure Group - Kenya",
    y = "Number of Populations Monitored",
    x = "Class",
    fill = "Protection and Exposure Status"
  ) +
  theme_minimal()













