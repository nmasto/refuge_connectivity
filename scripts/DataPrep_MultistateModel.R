
###############################################################################L

# This code preps the data for the multi-state model to estimate 
# movement between refuges for mallards in western Tennessee. The project is 
# part of a mallard project at Tennessee Tech University in collaboration with 
# Tennessee Wildlife Resources Agency. 

# Dependencies: "all_dat_oct2019_march2023.rds" & "spatial_sanctuary.shp"


# Allison C Keever & Nick Masto
# Tennessee Tech University
# akeever@tntech.edu
# GitHub: akeever2 & nmasto
# Created: 5/5/21
# Last updated: 8/21/23

# Last updated for reproducible data extraction 
# and accurate spatial sanctuaries 8/21/2023

###############################################################################L

# Libraries and data-------------

#library(raster)    # geodata
#library(terra)
library(sf)
#library(tmaptools)
library(tidyverse) # data mgt
#library(lubridate) 



# Clean dat isolate winters November-February

tn <- readRDS("data/all_dat_oct2019_march2023.rds")

tn2 <- tn %>% 
  distinct(Longitude, Latitude, Time, .keep_all = TRUE) %>% # keep distinct records (remove 2614 locs)
  group_by(BirdsID) %>%                  # group by id 
  filter(!is.na(sex)) %>%                # remove "TEST birds"
  mutate(month = month(date),            # add m/y/d columns separate
         year = year(date),              
         day = day(date),
         hour = hour(Time),
         ord = yday(date),
         yday = paste(year, sep = "_", ord)) %>%  # year/day combo
  filter(!between(month, 3, 10)) %>%     # limit to winter (months 11-2)
  filter(n() >= 50) %>%                  # remove birds w/ essentially no data
  ungroup() %>%
  mutate(study_year = case_when(date >= ymd("2019-10-31") & date < ymd("2020-03-01") ~ "yr1",  # create "study years"
                                date >= ymd("2020-10-31") & date < ymd("2021-03-01") ~ "yr2",  # not entirely necess 
                                date >= ymd("2021-10-31") & date < ymd("2022-03-01") ~ "yr3",  # b/c we do it later
                                date >= ymd("2022-10-31") & date < ymd("2023-03-01") ~ "yr4")) %>% # but oh well; part 
  mutate(BirdsID_season = # Separate IDs for multiple study year (not real year)               # of nick's normal workflow
           paste(BirdsID, sep = "_", study_year)) %>% 
  group_by(BirdsID_season) %>%
  #filter(study_year != "yr1") %>%
  mutate(no_months = n_distinct(month), # Number of distinct months 
         no_days = n_distinct(ord),     # Number of days
         obs = 1:n()                    # Number rows for each id
         ) %>%
  ungroup() %>% 
  arrange(BirdsID, date) %>% # Omit first 4 days per Methods
  group_by(BirdsID) %>%
  mutate(no_seq_days = data.table::rleid(date)) %>%
  filter(no_seq_days > 4) %>%
  ungroup()

# Read in Refuges and overlay locations------------------

Refs <- read_sf("geo_data/spatial_sanctuary.shp") %>% 
  st_transform("+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs +type=crs")

tn_sf <- st_as_sf(tn2, coords = c("Longitude", "Latitude"), crs = "epsg:4326") %>%
  mutate(lon = sf::st_coordinates(.)[,1],
         lat = sf::st_coordinates(.)[,2]) %>%
  st_transform("+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs +type=crs")

# coords <- st_coordinates(tn_sf)
# tn_sf <- mutate(tn_sf, utm.easting = coords[,"X"], utm.northing = coords[,"Y"])

# Let's use st_join and our bounding box. A little messier but allows us to track
# The number of birds we're removing. st_intersection would remove NAs and be
# slightly cleaner but can't track the birds removed.

tn_sf_refs <- tn_sf %>% st_join(Refs) 

# Tracking no of birds as we filter
#tn_sf_refs %>% group_by(BirdsID)        # 452
#tn_sf_refs %>% group_by(BirdsID_season) # 569

# Read in and intersect with bounding box--------------

bbox <- read_sf("geo_data/tn_bbox.shp") %>% 
  st_transform("+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs +type=crs")

tn_sf_refs_bbox <- tn_sf_refs %>% st_intersection(bbox)

# bbox removed 9 birds and 86 "return" birds
#tn_sf_refs_bbox %>% group_by(BirdsID)        # 452 - 443 = 9
#tn_sf_refs_bbox %>% group_by(BirdsID_season) # 569 - 483 = 86


# Change a few things for structure / formatting------------- 
# Filter out NAs because those locations were not on Refuge

datum <- tn_sf_refs_bbox %>% 
  dplyr::select(Time, BirdsID, BirdsID_season, site, sex, age, mass, tarsus,
                date, month, year, day, hour, study_year, lon, lat,
                Name, Area_sqkm, Type, geometry) %>%
  rename(Refuge = Name, trackID = BirdsID, timestamp = Time, yr = year) %>% 
  filter(!is.na(Refuge), trackID != c("TEST3")) #%>% 
  #datum$timestamp <- as_datetime(datum$timestamp)

#saveRDS(datum, "data/datum.rds")
  
# Last check on number of indiviudals
datum %>% group_by(trackID)        # 398
datum %>% group_by(BirdsID_season) # 421

# Removing NAs (i.e., no Refuge use) removed 45 individuals
# And 62 indiviudals if we include "return wintering birds".

# Create new occasion variable for each date and day to get daily movements. 
# First create a dataframe that holds the occassion number
occ_datum <- data.frame(occday = seq.Date(make_date(year = 2019, month = 11, day = 5), 
                                          make_date(year = 2020, month = 2, day = 29), 
                                          1), 
                        occ = 1:117)

# Then we are going to take only the first location for each date-day occasion
# from 1000 to 1900 from November to February. 

datum <- datum %>% 
  #mutate(yr = year(timestamp), 
  #       month = month(timestamp), 
  #       hour = hour(timestamp)) %>%
  filter(month %in% c(1, 2, 11, 12), hour %in% c(10:19)) %>%
  mutate(dateday = make_date(year = yr, month = month, day = day(timestamp)),
         occday = make_date(year = case_when(month %in% c(1, 2) ~ 2020, 
                                             TRUE ~ 2019),
                            month = month, day = day(timestamp)), 
         season = case_when(
           yr == 2019 ~ 1, 
           yr == 2020 & month %in% c(1, 2) ~ 1, 
           yr == 2020 & month %in% c(11, 12) ~ 2, 
           yr == 2021 & month %in% c(1, 2) ~ 2,
           yr == 2021 & month %in% c(11, 12) ~ 3, 
           yr == 2022 & month %in% c(1, 2) ~ 3,
           yr == 2022 & month %in% c(11, 12) ~ 4, 
           yr == 2023 & month %in% c(1, 2) ~ 4),
         Refuge = case_when(
           Refuge == "Lake Isom" ~ "LINWR",
           Refuge == "Bean Switch" ~ "BS", 
           Refuge == "HopIn Refuge" ~ "HI", 
           Refuge == "Maness Swamp" ~ "M", 
           Refuge == "White Lake Refuge" ~ "WL", 
           Refuge == "Chickasaw" ~ "CNWR",
           Refuge == "Lake Lauderdale" ~ "LL", 
           Refuge == "Reelfoot North" ~ "RLNWR_N",
           Refuge == "Reelfoot South" ~ "RLNWR_S", 
           Refuge == "Big Lake" ~ "BLNWR",
           Refuge == "Horns Bluff" ~ "HB",
           Refuge == "Black Bayou" ~ "BB",
           Refuge == "Phillipy" ~ "P",
           TRUE ~ Refuge)) %>%
  group_by(trackID, dateday) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(occ_datum, by = "occday") 

# Data need to be set up in matrix with nrows = # individuals and ncols = # occs. 
# The values will be the state (i.e., Refuge) for each individual at each 
# occasion. Replace Refuge name with number...

unique(datum$Refuge) # Unique states

# RLNWR_N = 1, RLNWR_S = 2, LINWR = 3, WL = 4, CNWR = 5, HI = 6, M = 7, P = 8, 
# BB = 9, BLNWR = 10, BS = 11, LL = 12, HB = 13, not seen = 14. 

enc_hist <- datum %>% as_tibble() %>% # It was an sf object following update
  arrange(occday) %>% 
  pivot_wider(id_cols = c(trackID, season, age, sex), names_from = occ, names_prefix = "occ_", 
              values_from = Refuge, values_fill = NA) %>%
  mutate(across(starts_with("occ_"), ~ case_when(
    . == "BB" ~ "1", 
    . == "LINWR" ~ "2", 
    . == "BS" ~ "3", 
    . == "HI" ~ "4", 
    . == "M" ~ "5", 
    . == "WL" ~ "6", 
    . == "CNWR" ~ "7", 
    . == "BLNWR" ~ "8", 
    . == "LL" ~ "9", 
    . == "RLNWR_N" ~ "10", 
    . == "HB" ~ "11", 
    . == "RLNWR_S" ~ "12",
    . == "P" ~ "13",
    is.na(.) ~ "14", 
    TRUE ~ .))) %>%
  mutate(across(starts_with("occ_"), as.numeric), 
         across(starts_with("occ"), ~ replace(., is.na(.), 14))) %>%
  arrange(season)

# Get occasion (date and hour) of first capture
get_first <- function(x) min(which(x != 14))
f <- apply(enc_hist %>% select(starts_with("occ_")), 1, get_first)


# Need to get distances and sizes of Refuges
Refs <- Refs  %>%
  mutate(Name = case_when(
    Name == "Lake Isom" ~ "LINWR",
    Name == "Bean Switch" ~ "BS", 
    Name == "HopIn Refuge" ~ "HI", 
    Name == "Maness Swamp" ~ "M", 
    Name == "White Lake Refuge" ~ "WL", 
    Name == "Chickasaw" ~ "CNWR",
    Name == "Lake Lauderdale" ~ "LL", 
    Name == "Reelfoot North" ~ "RLNWR_N",
    Name == "Reelfoot South" ~ "RLNWR_S", 
    Name == "Big Lake" ~ "BLNWR",
    Name == "Horns Bluff" ~ "HB",
    Name == "Black Bayou" ~ "BB",
    Name == "Phillipy" ~ "P",
    TRUE ~ Name))

# Delete little river CA, MO which was never used
unique(Refs$Name)

Refs <- Refs[-8,] # Removing Little River Refuge, MO

# Ordered Refuges same as "states"
Refsorder <- c("BB", "LINWR", "BS", "HI", "M", "WL", "CNWR", 
               "BLNWR", "LL", "RLNWR_N", "HB", "RLNWR_S", "P")
Refs <- Refs[order(match(Refs$Name, Refsorder)),]

# Get distances between 
distances <- st_distance(Refs, Refs)
colnames(distances) <- rownames(distances) <- Refs$Name
dist_datum <- units::drop_units(as.data.frame(distances)) %>%
  rownames_to_column("From") %>%
  pivot_longer(cols = 2:14, names_to = "To", values_to = "Distance") %>%
  mutate(Distance = Distance * 0.001) %>%
  filter(Distance != 0) %>% 
  pivot_wider(names_from = To, values_from = Distance) %>% 
  select(BB, everything(), -From)

dist_datum <- t(matrix(t(dist_datum)[which(!is.na(dist_datum))], nrow = 12, 
                       ncol = 13)) 
#str(dist_datum)
#saveRDS(dist_datum, "data/dist_datum.rds")


d <- data.frame(dist = c(11.2, 58.4, 39.8, 48.4, 34.6, 63.1, 79.4, 55.9, 4.98,
                         66.6, 1.35, 2.08, 55.5, 37.2, 43.9, 17.7, 45.7, 68.9,
                         38.5, 16.9, 52.1, 9.21, 15.7, 16.1, 8.70, 65.1, 76.9,
                         120, 61.0, 54.5, 44.9, 51.7, 58.6, 6.15, 49.1, 65.7,
                         104, 50.7, 36.4, 42.2, 33.3, 40.2, 53.0, 66.2, 108,
                         50.6, 45.5, 35.6, 42.2, 49.1, 24.2, 48.9, 24.0, 41.1,
                         45.8, 33.4, 39.2, 44.5, 14.6, 69.2, 38.9, 61.5, 67.7,
                         62.6, 87.2, 87.6, 81.2, 83.7, 60.7, 23.3, 53.1, 60.3,
                         68.7, 1.65, 1.80, 61.7, 69.7, 2.56))
summary(d$dist)
median(d$dist)
boxplot(d$dist)
