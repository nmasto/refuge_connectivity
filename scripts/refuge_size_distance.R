###################################################################l 
#
#  Size and distance to refuges
#
#  Calculate size and distance to and among refuges 
#
#  Ally Keever, modified by Nick Masto 7/21/2023
#
#  Last updated:
# 
###################################################################l

# Read in libraries

library(tidyverse)
library(terra)
library(sf)

# Read in Data
AR_refs <- read_sf(dsn = "data", layer = "PADUS3_0Combined_StateAR")
TN_refs <- read_sf(dsn = "data", layer = "PADUS3_0Combined_StateTN")

refs <- read_sf("C:/Users/nmmasto42/Documents/TennesseeTech/ALLY_KEEVER/spatial_sanctuary.shp")
plot(refs)
chick <- TN_refs %>% filter(Unit_Nm == "Chickasaw National Wildlife Refuge")
reelfoot <- TN_refs %>% filter(Unit_Nm == "Reelfoot National Wildlife Refuge") 
plot(reelfoot[,"Loc_Ds"])

# Filter to the ones we need
TN_refs <- TN_refs %>%
  filter(Unit_Nm %in% c("Black Bayou Refuge", 
                        "Lake Isom National Wildlife Refuge",
                        "Bean Switch Refuge", "Hop-In Refuge", 
                        "Maness Swamp Refuge", "White Lake Refuge", 
                        "Chickasaw National Wildlife Refuge",
                        "Lake Lauderdale Refuge", 
                        "Reelfoot National Wildlife Refuge",
                        "Horns Bluff Refuge"), 
         Category == "Fee") %>%
  filter(!(Unit_Nm == "Reelfoot National Wildlife Refuge" & Loc_Ds == "NWR"))

plot(TN_refs)

AR_refs <- AR_refs %>%
  filter(Unit_Nm %in% c("Big Lake National Wildlife Refuge"),
         Des_Tp == "NWR")

# Clean up the refuges
Refs <- bind_rows(TN_refs %>%
                    rename("Name" = "Unit_Nm") %>%
                    filter(Name != "Black Bayou Refuge"), 
                  AR_refs %>%
                    rename("Name" = "Unit_Nm")) %>% 
  dplyr::select(Name, geometry) %>%
  bind_rows(st_union(TN_refs %>% 
                       filter(Unit_Nm == "Black Bayou Refuge") %>% 
                       select(Unit_Nm, geometry)) %>% 
              st_sf() %>%
              mutate(Name = "BB")) %>%
  mutate(Area_sqkm = units::drop_units(st_area(.) * 0.000001), 
         Name = case_when(
           Name == "Lake Isom National Wildlife Refuge" ~ "LINWR",
           Name == "Bean Switch Refuge" ~ "BS", 
           Name == "Hop-In Refuge" ~ "HI", 
           Name == "Maness Swamp Refuge" ~ "M", 
           Name == "White Lake Refuge" ~ "WL", 
           Name == "Chickasaw National Wildlife Refuge" ~ "CNWR",
           Name == "Lake Lauderdale Refuge" ~ "LL", 
           Name == "Reelfoot National Wildlife Refuge" ~ "RLNWR",
           Name == "Horns Bluff Refuge" ~ "HB", 
           Name == "Big Lake National Wildlife Refuge" ~ "BLNWR",
           TRUE ~ Name))

# Order of the refuges
reforder <- c("BB", "LINWR", "BS", "HI", "M", "WL", "CNWR", "BLNWR", 
              "LL", "RLNWR", "HB")
Refs <- Refs[order(match(Refs$Name, reforder)),]

# Get distances between 
distances <- st_distance(Refs, Refs)
colnames(distances) <- rownames(distances) <- Refs$Name
dist_datum <- units::drop_units(as.data.frame(distances)) %>%
  rownames_to_column("From") %>%
  pivot_longer(cols = 2:11, names_to = "To", values_to = "Distance") %>%
  mutate(Distance = Distance * 0.001) %>%
  filter(Distance != 0) %>% 
  pivot_wider(names_from = To, values_from = Distance) %>% 
  select(BB, everything(), -From)
dist_datum <- t(matrix(t(dist_datum)[which(!is.na(dist_datum))], nrow = 10, 
                       ncol = 11))


#########################################lll
# Read in Data
refs <- read_sf("C:/Users/nmmasto42/Documents/TennesseeTech/ALLY_KEEVER/spatial_sanctuary.shp")

names(refs)
# Clean up the refuges
Refs <- refs %>%
  mutate(Name = case_when(
    Name == "Bean Switch" ~ "BS",
    Name == "Big Lake" ~ "BL",
    Name == "Black Bayou" ~ "BB", 
    Name == "Chickasaw" ~ "CNWR", 
    Name == "HopIn Refuge" ~ "HI", 
    Name == "Horns Bluff" ~ "HB", 
    Name == "Lake Isom" ~ "LI",
    Name == "Lake Lauderdale" ~ "LL", 
    Name == "Little River CA" ~ "LR",
    Name == "Maness Swamp" ~ "MS",
    Name == "Phillipy" ~ "PH",     # Phillipy on its own? Or separate from BB?
    Name == "Reelfoot" ~ "RLNWR",
    Name == "White Lake Refuge" ~ "WLR",
    TRUE ~ Name)) %>% dplyr::select(Name, Area_sqkm, geometry)

# Order of the refuges
reforder <- c("BB", "PH", "LINWR", "LR", "RLNWR", "HI", "BS", "MS", 
              "WL", "CNWR", "BLNWR", "LL", "HB")
Refs <- Refs[order(match(Refs$Name, reforder)),]

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

median(c(7.34, 22.7, 2.96, 2.64, 7.2, 2.64, 4.93, 4.94, 44.2, 5.56))
sd(c(7.34, 22.7, 2.96, 2.64, 7.2, 2.64, 4.93, 4.94, 44.2, 5.56))

min(c(2.08, 52.0, 1.35, 39.9, 58.5, 48.4, 63.1,  
      55.9, 66.7, 79.4, 11.2, 34.6, 55.9, 1.80, 
      40.3, 58.7, 49.1, 67.8, 60.3, 69.8, 83.7, 
      15.7, 39.3, 55.4, 86.3, 104,  91.7, 53.9, 
      62.6, 86.9, 31.5, 46.5, 37.8, 33.3, 51.7, 
      42.2, 61.5, 53.1, 61.7, 81.3, 9.21, 33.4,
      16.1, 6.16, 65.8, 50.7, 42.2, 104,  37.3,   
      49.1, 8.71, 77.0, 61.1, 44.9, 120,  55.6, 
      65.1, 66.2, 50.7, 35.6, 108,  43.9, 53.0,
      14.7, 38.9, 44.5, 45.7, 24.2, 23.3, 62.6, 
      38.5, 24.0, 87.7, 52.2, 45.8, 68.9, 48.9, 17.7))
sd(c(2.08, 52.0, 1.35, 39.9, 58.5, 48.4, 63.1,  
     55.9, 66.7, 79.4, 11.2, 34.6, 55.9, 1.80, 
     40.3, 58.7, 49.1, 67.8, 60.3, 69.8, 83.7, 
     15.7, 39.3, 55.4, 86.3, 104,  91.7, 53.9, 
     62.6, 86.9, 31.5, 46.5, 37.8, 33.3, 51.7, 
     42.2, 61.5, 53.1, 61.7, 81.3, 9.21, 33.4,
     16.1, 6.16, 65.8, 50.7, 42.2, 104,  37.3,   
     49.1, 8.71, 77.0, 61.1, 44.9, 120,  55.6, 
     65.1, 66.2, 50.7, 35.6, 108,  43.9, 53.0,
     14.7, 38.9, 44.5, 45.7, 24.2, 23.3, 62.6, 
     38.5, 24.0, 87.7, 52.2, 45.8, 68.9, 48.9, 17.7))