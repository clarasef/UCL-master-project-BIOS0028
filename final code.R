library(lubridate)
library(sf)
library(units)
library(purrr)
library(stringr)
library(ggplot2)

# Reshape dataset of Africa to long format, rename "value" to "popvalue", convert to numeric
livingplanetindex <- lpi_pops_Africa_20250508 %>%
  filter(System != "Marine") %>%
  pivot_longer(cols = `1950`:`2022`, names_to = "year")
names(livingplanetindex)[names(livingplanetindex) == "value"] <- "popvalue"
livingplanetindex$year <- as.numeric(as.character(livingplanetindex$year))
livingplanetindex$popvalue <- as.numeric(as.character(livingplanetindex$popvalue))

# 1. CLEAN AND PREPARE CONFLICT DATA (ACLED)
acled <- acled_africa %>%
  filter(!is.na(latitude), !is.na(longitude), event_type == "Battles") %>%
  mutate(
    date = ymd(event_date),
    year = year(date),
    intensity = replace_na(fatalities, 0),
    event_id = row_number()
  ) %>%
  select(year, latitude, longitude, intensity, event_id)

# 2. ADD A TEMPORAL LAG (e.g., 3 years after each event)
lag_years <- 3
conflict_lagged <- acled %>%
  rowwise() %>%
  mutate(lagged_year = list((year + 1):(year + lag_years))) %>%
  unnest(lagged_year) %>%
  ungroup()

#Convert the lagged conflict events to spatial points (simple features), with WGS84 geographic CRS.
conflict_lagged_sf <- st_as_sf(conflict_lagged, coords = c("longitude", "latitude"), crs = 4326)
#Then reproject to metric CRS (EPSG 3395) for distance calculations.
conflict_proj <- st_transform(conflict_lagged_sf, 3395)

# 3. CLEAN POPULATION DATA
#Keep only pops with coordinates, convert to spatial points, 
pop_sf <- livingplanetindex %>%
  filter(!is.na(Longitude), !is.na(Latitude)) %>%
  mutate(species = tolower(str_trim(Binomial))) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

#Reproject to metric CRS for spatial analysis.
pop_sf_proj <- st_transform(pop_sf, 3395)

# 4. GET START & END YEARS
#For each unique population (ID and species), find the earliest and latest years with data
pop_series_years <- livingplanetindex %>%
  group_by(ID, Binomial) %>%
  summarise(
    start_year = min(year, na.rm = TRUE),
    end_year = max(year, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(species = tolower(str_trim(Binomial)))

# JOIN YEARS TO SPATIAL DATA
#Add start/end year info to the spatial points by joining on population ID and species
pop_locations <- pop_sf_proj %>%
  mutate(species = tolower(str_trim(Binomial))) %>%
  left_join(pop_series_years, by = c("ID", "species"))

# 5. CONFLICT EXPOSURE (SPATIAL + TEMPORAL)
#Create a 100 km buffer around each population location (to consider nearby conflicts)
pop_buffered <- st_buffer(pop_locations, dist = 100000) # 100 km
#Find which conflict points fall inside these buffers (spatial intersection)
intersections <- st_intersects(pop_buffered, conflict_proj, sparse = TRUE)

#For each population, check if any nearby conflict event (within 100 km buffer) occurred within the population’s observed time span (start to end year).
#Returns TRUE if conflict exposure exists, FALSE otherwise.
conflict_exposure_logical <- map_lgl(seq_along(intersections), function(i) {
  site_conflicts <- conflict_proj[intersections[[i]], ]
  if (nrow(site_conflicts) == 0) return(FALSE)
  
  start <- pop_buffered$start_year[i]
  end <- pop_buffered$end_year[i]
  
  any(site_conflicts$lagged_year >= start & site_conflicts$lagged_year <= end)
})

# 6. CONFLICT EXPOSURE SUMMARY PER POPULATION ID
#Summarize conflict exposure as a boolean flag per population ID
#If no info, treat as not exposed (FALSE)
pop_with_conflict_clean <- tibble(
  ID = pop_buffered$ID,
  conflict_exposed = conflict_exposure_logical
) %>%
  group_by(ID) %>%
  summarise(conflict_exposed = any(conflict_exposed), .groups = "drop") %>%
  mutate(conflict_exposed = replace_na(conflict_exposed, FALSE))

# 7. DISTANCE-WEIGHTED CONFLICT INTENSITY (PER YEAR/POPULATION)
# Batched spatial join to reduce memory usage 
batch_size <- 1
years <- sort(unique(livingplanetindex$year))
year_batches <- split(years, ceiling(seq_along(years) / batch_size))

conflicts_nearby_list <- list()

#For each year (in batches):
#Select populations present in that year
#Join with conflict events for the same lagged year within 100 km
#Calculate distance between population and conflict
#Assign weights inversely proportional to distance (closer conflicts count more)
#Calculate weighted conflict intensity (intensity * weight)
#Summarize total weighted conflict and number of conflict events near each population per year
#Store results in list and perform garbage collection to free memory
for (i in seq_along(year_batches)) {
  message("Processing batch ", i, " of ", length(year_batches))
  
  batch_result <- map(year_batches[[i]], function(this_year) {
    pops_year <- livingplanetindex %>%
      filter(year == this_year) %>%
      mutate(species = tolower(str_trim(Binomial))) %>%
      left_join(
        pop_sf_proj %>% mutate(species = tolower(str_trim(Binomial))) %>% st_drop_geometry(),
        by = c("ID", "species")
      ) %>%
      st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
      st_transform(3395)
    
    # Rename conflict geometry to preserve it after join
    conflicts_year <- conflict_proj %>%
      filter(lagged_year == this_year) %>%
      mutate(conflict_geom = geometry)  
    
    nearby <- st_join(pops_year, conflicts_year, join = st_is_within_distance, dist = 100000, left = FALSE)
    if (nrow(nearby) == 0) return(NULL)
    
    nearby <- nearby %>%
      mutate(
        dist_m = st_distance(geometry, conflict_geom, by_element = TRUE),  # Use renamed column
        dist_km = set_units(dist_m, "km") %>% drop_units(),
        weight = 1 / (1 + dist_km)
      ) %>%
      filter(dist_km <= 100)
    
    nearby %>%
      st_drop_geometry() %>%
      group_by(ID, species, year = this_year) %>%
      summarise(
        weighted_conflict = sum(intensity * weight, na.rm = TRUE),
        conflict_count = n(),
        .groups = "drop"
      )
  })
  
  conflicts_nearby_list[[i]] <- bind_rows(batch_result)
  gc()
}
#Combine all batches into one dataframe
conflicts_filtered <- bind_rows(conflicts_nearby_list)


# 8. MERGE INTO FINAL DATA
#Clean and prepare population data with species and IDs
#Join weighted conflict intensity data per population/year.
#Join overall conflict exposure flag per population.
#Create a new categorical variable "conflict_group" based on exposure.
livingplanetindex2 <- livingplanetindex %>%
  mutate(
    ID = as.integer(ID),
    species = tolower(str_trim(Binomial))
  )

pop_final <- livingplanetindex2 %>%
  left_join(conflicts_filtered, by = c("ID", "species", "year")) %>%
  left_join(pop_with_conflict_clean, by = "ID") %>%
  mutate(
    conflict_group = ifelse(conflict_exposed, "Conflict-exposed", "No conflict")
  )

setwd("/Users/clarasefzig/Documents/Project ZSL/new_master_proj")
# 9. LPI – NO CONFLICT
#Calculate LPI for populations not affected by conflict
pop_no_conflict <- pop_final %>%
  filter(conflict_group == "No conflict")

pop_no_conflict_wide <- pop_no_conflict %>%
  select(ID, Binomial = species, Country, year, popvalue) %>%
  mutate(
    Binomial = str_trim(tolower(Binomial)),
    year = paste0("X", year)
  ) %>%
  pivot_wider(names_from = year, values_from = popvalue)

create_infile(
  pop_data_source = pop_no_conflict_wide,
  name = "Africa_LPI_no_conflict",
  start_col_name = "X1950",
  end_col_name = "X2022"
)

result_no_conflict <- LPIMain("Africa_LPI_no_conflict_infile.txt")
result_no_conflict$Year <- as.numeric(row.names(result_no_conflict))
result_no_conflict_clean <- result_no_conflict %>%
  filter(!is.na(Year) & is.finite(Year) & is.finite(LPI_final))

# 10. LPI – CONFLICT-EXPOSED
#Calculate LPI for populations affected by conflict
pop_conflict_exposed <- pop_final %>%
  filter(conflict_group == "Conflict-exposed")

pop_conflict_exposed_wide <- pop_conflict_exposed %>%
  select(ID, Binomial = species, Country, year, popvalue) %>%
  mutate(
    Binomial = str_trim(tolower(Binomial)),
    year = paste0("X", year)
  ) %>%
  pivot_wider(names_from = year, values_from = popvalue)

create_infile(
  pop_data_source = pop_conflict_exposed_wide,
  name = "Africa_LPI_conflict_filtered",
  start_col_name = "X1950",
  end_col_name = "X2022"
)

result_conflict_filtered <- LPIMain("Africa_LPI_conflict_filtered_infile.txt")
result_conflict_filtered_clean <- result_conflict_filtered %>%
  mutate(Year = as.numeric(row.names(.))) %>%
  filter(!is.na(Year) & is.finite(LPI_final))

# 11. LPI – ALL POPULATIONS
#Calculate LPI for all populations 
pop_all_wide <- pop_final %>%
  select(ID, Binomial = species, Country, year, popvalue) %>%
  mutate(
    Binomial = str_trim(tolower(Binomial)),
    year = paste0("X", year)
  ) %>%
  pivot_wider(names_from = year, values_from = popvalue)

create_infile(
  pop_data_source = pop_all_wide,
  name = "Africa_LPI_all",
  start_col_name = "X1950",
  end_col_name = "X2022"
)

result_all <- LPIMain("Africa_LPI_all_infile.txt")
result_all_clean <- result_all %>%
  mutate(Year = as.numeric(row.names(.))) %>%
  filter(!is.na(Year) & is.finite(LPI_final))

# Add group labels
df_no_conflict <- result_no_conflict_clean %>%
  mutate(Group = "No conflict")

df_conflict <- result_conflict_filtered_clean %>%
  mutate(Group = "Conflict-exposed")

df_all <- result_all_clean %>%
  mutate(Group = "All populations")


#ADDING CLIMATE EXPOSURE
library(readr)
library(dplyr)


# Select relevant disasters
relevant_disasters <- c("Drought", "Flood", "Extreme temperature", "Storm", "Wildfire")

climate.data <- climate.data %>%
  filter(Disaster.Type %in% relevant_disasters) %>%
  mutate(
    start_year = as.numeric(Start.Year),
    end_year = as.numeric(End.Year),
    latlon_available = !is.na(Latitude) & !is.na(Longitude)
  )
# Convert to sf for spatial disasters
climate_spatial <- climate.data %>%
  filter(latlon_available) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(3395)
#Fallback for disaster with only country information
climate_country <- dplyr::select(climate.data, Country, Start.Year, End.Year)

#Spatial join for disaster with coordinates
intersections_climate <- st_intersects(pop_buffered, climate_spatial, sparse = TRUE)

climate_exposure_spatial <- map_lgl(seq_along(intersections_climate), function(i) {
  site_disasters <- climate_spatial[intersections_climate[[i]], ]
  if (nrow(site_disasters) == 0) return(FALSE)
  
  start <- pop_buffered$start_year[i]
  end <- pop_buffered$end_year[i]
  
  any(site_disasters$start_year <= end & site_disasters$end_year >= start)
})

#National join for disasters without coordinates 
climate_exposure_country <- map_lgl(seq_len(nrow(pop_buffered)), function(i) {
  country <- pop_buffered$Country[i]
  start <- pop_buffered$start_year[i]
  end <- pop_buffered$end_year[i]
  
  disasters <- climate_country %>% filter(str_to_lower(country) == str_to_lower(pop_buffered$Country[i]))
  if (nrow(disasters) == 0) return(FALSE)
  
  any(disasters$start_year <= end & disasters$end_year >= start)
})

#Combine exposures
climate_exposure_combined <- climate_exposure_spatial | climate_exposure_country

pop_climate_exposure <- tibble(
  ID = pop_buffered$ID,
  climate_exposed = climate_exposure_combined
)
#Merge with populations exposed to conflict
pop_exposure_combined <- pop_with_conflict_clean %>%
  left_join(pop_climate_exposure, by = "ID") %>%
  mutate(climate_exposed = ifelse(is.na(climate_exposed), FALSE, climate_exposed))

#Add exposure group to populations
pop_exposure_combined <- pop_exposure_combined %>%
  distinct(ID, .keep_all = TRUE)


pop_final <- pop_final %>%
  left_join(pop_exposure_combined, by = "ID") %>%
  mutate(
    exposure_group = case_when(
      conflict_exposed.y & climate_exposed ~ "Both conflict + climate",
      conflict_exposed.y ~ "Conflict only",
      TRUE ~ "No exposure"
    )
  )

#Plot LPI by exposure group
library(tidyr)
library(purrr)
library(stringr)
library(rlpi)


pop_grouped_lpi <- function(group_name) {
  pop_subset <- pop_final %>%
    filter(exposure_group == group_name) %>%
    select(ID, Binomial = species, Country, year, popvalue) %>%
    mutate(
      Binomial = str_trim(tolower(Binomial)),
      year = paste0("X", year)
    ) %>%
    pivot_wider(names_from = year, values_from = popvalue)
  
  create_infile(
    pop_data_source = pop_subset,
    name = paste0("Africa_LPI_", str_replace_all(group_name, " ", "_")),
    start_col_name = "X1950",
    end_col_name = "X2022"
  )
  
  result <- LPIMain(paste0("Africa_LPI_", str_replace_all(group_name, " ", "_"), "_infile.txt"))
  result_clean <- result %>%
    mutate(Year = as.numeric(row.names(.))) %>%
    filter(!is.na(Year) & is.finite(LPI_final)) %>%
    mutate(Group = group_name)
  
  return(result_clean)
}

groups <- unique(pop_final$exposure_group)
groups <- setdiff(groups, "Climate only")
# Remove "Climate only" from the exposure groups
groups <- setdiff(unique(pop_final$exposure_group), "Climate only")
# Run the pop_grouped_lpi() for each group and combine results
lpi_grouped_exposure <- bind_rows(map(groups, pop_grouped_lpi))

# Add the overall LPI ("All populations")
lpi_combined <- bind_rows(
  lpi_grouped_exposure,
  df_all
)

df_conflict_climate <- pop_grouped_lpi("Both conflict + climate") %>%
  mutate(Group = "Both conflict + climate")

df_conflict_climate <- df_conflict_climate %>%
  filter(Year < 2016)

lpi_combined <- bind_rows(df_no_conflict, df_conflict, df_all, df_conflict_climate)
#Add percentage labels 
evolution_summary <- lpi_combined %>%
  group_by(Group) %>%
  summarise(
    first_year = min(Year),
    last_year = max(Year),
    first_lpi = LPI_final[Year == first_year][1],
    last_lpi = LPI_final[Year == last_year][1],
    percent_change = 100 * (last_lpi - first_lpi) / first_lpi
  )
ggplot(lpi_combined, aes(x = Year, y = LPI_final, color = Group)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c(
    "No conflict" = "chartreuse3",
    "Conflict-exposed" = "cyan3",
    "All populations" = "chocolate1",
    "Both conflict + climate" = "darkred"
  )) +
  labs(
    title = "Living Planet Index by Conflict Exposure - Africa",
    x = "Year",
    y = "LPI (Relative Abundance)",
    color = "Group"
  ) +
  theme_minimal(base_size = 14) +
  ylim(0, NA)+
  theme_minimal(base_size = 13) +
  ylim(0, NA) +
  theme(
    plot.title = element_text(face = "bold")
  ) +
  # Add labels at the end of each line with % change
  geom_text(data = evolution_summary, aes(
    x = last_year + 0.5,
    y = last_lpi,
    label = paste0(round(percent_change, 1), "%"),
    color = Group
  ), hjust = 0) 

