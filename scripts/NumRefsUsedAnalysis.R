#############################################################################L

# Analysis to look at time effect of using multiple refuges

##############################################################################L


# Load packages
library(tidyverse)
library(lubridate)
library(sf)
library(tidybayes)
library(glmmTMB)
library(DHARMa)

# Read in data and change a few things for structure / formatting. Filter out 
# locations above 36.5 latitude because they were not in our study area
datum <- readRDS("data/datum.rds") %>% 
  #rename(Refuge = refuge, timestamp = t_, trackID = trackId) %>% 
  filter(#!(trackID %in% c("TEST3", "WLR.b022021.M", "HNWR.111.M", "LINWR.b1232020.M")), 
         lat <= 36.55, lat >= 35.65, lon <= -88.7, lon >= -90.2) 

#datum$timestamp <- as_datetime(datum$timestamp)
str(datum)

# # Make SF object
# sf_datum <- st_as_sf(datum, coords = c("lon", "lat"))
# plot(sf_datum[,2])
# plot(tn, add = TRUE)

# Set up data to get the number of refuges used by individual each year, 
# the number of days spend in TN area, and the first date of arrival/capture. 
datum <- datum %>% 
  filter(month(timestamp) %in% c(1,2,11,12)) %>% 
  mutate(Year = year(timestamp), 
         SeasonYear = case_when(
           month(timestamp) %in% c(1,2,3) ~ Year - 1, 
           TRUE ~ Year)) %>% 
  group_split(trackID) %>%
  map_df(~ .x %>% group_by(SeasonYear) %>% 
           mutate(bird_yr = cur_group_id()) %>% 
           ungroup()) %>%
  group_by(trackID, SeasonYear) %>%
  mutate(first_arrival = min(timestamp), 
         last_location = max(timestamp),
         num_days = (max(timestamp) - min(timestamp)) / ddays(1)) %>%
  summarise(refs_used = length(unique(Refuge[!is.na(Refuge)])), 
            first_arrival = min(first_arrival), 
            last_location = max(last_location), 
            num_days = min(num_days), #sum(unique(num_days)), 
            first_month = factor(month(min(first_arrival)), levels = c("11", "12", "1", "2")),
            num_yrs = factor(first(bird_yr), levels = c("1", "2")),
            sex = first(sex),
            age = first(age))

table(datum$sex)
table(datum$age)
table(datum$num_yrs)
table(datum$refs_used)
mean(datum$num_days)

# Now run a model with # refuges used ~ number of days in area + first month here
m1 <- glmmTMB(refs_used ~ num_days + first_month + num_yrs, 
              family = "truncated_poisson", 
              data = datum)
summary(m1)

p1 <- ggplotify::as.ggplot(~plot(simulateResiduals(fittedModel = m1, plot = TRUE)))
ggsave(filename = "Figures/QQplots.png", dpi = 600, width = 7, height = 5)


data.frame(confint(m1)) %>% rename(low = "X2.5..",
                                   high = "X97.5..") %>%
  mutate(odds_est_30d = exp(Estimate * 30),
         odds_low_30d = exp(low * 30),
         odds_high_30d = exp(high * 30))

summary(datum$num_days)

as.data.frame(
  predict(object = m1,
              newdata = expand_grid(num_days = seq(0, 120, 5),
                                    first_month = factor(c(1, 2, 11, 12), levels = c("11", "12", "1", "2")), 
                                    num_yrs = factor(c(1,2))), 
              type = "response", 
              se.fit = TRUE)) %>% 
  bind_cols(expand_grid(num_days = seq(0, 120, 5), 
                        first_month = factor(c(1, 2,  11, 12), levels = c("11", "12", "1", "2")), 
                        num_yrs = factor(c(1,2), levels = c("1", "2")))) %>%
  ggplot(aes(x = num_days, y = fit, ymin = fit - 1.96 * se.fit, 
             ymax = fit + 1.96 * se.fit))+
  geom_lineribbon(alpha = .2) +
  geom_line(color = "black", size = 0.7) +
  facet_grid(num_yrs~first_month, 
             labeller = labeller(num_yrs = c('1' = "Capture year birds", 
                                             '2' = "Return birds"),
                                 first_month = c("11" = "Captured/arrived \nin November", 
                                                 "12" = "Captured/arrived \nin December",
                                                 "1" = "Captured/arrived \nin January",
                                                 "2" = "Captured/arrived \nin February")))+
  scale_x_continuous(name = "Number of days within sanctuary network", 
                     breaks = seq(0, 120, 20)) +
  scale_y_continuous(name = "Number of sanctuaries used", breaks = seq(0, 5, 1))+
  scale_fill_brewer(palette = "Greys") +
  coord_cartesian(ylim = c(0, 5), xlim = c(0, 120)) +
  theme_classic() + 
  theme(text = element_text(size = 13), 
        legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA))


ggsave(filename = "Figures/numRefs.png", dpi = 600, width = 7.5, height = 5)


ggplot(data = datum, aes(refs_used)) + 
  geom_bar(aes(fill = num_yrs)) +
  scale_x_continuous(name = "Number of refuges used", breaks = seq(0, 7, 1)) + 
  scale_y_continuous(name = "Count", breaks = seq(0, 500, 25))+
  scale_fill_brewer(palette = "Set2", name = "Number of seasons") +
  theme_classic() + 
  theme(text = element_text(size = 12), 
        legend.position = c(.8,.9), 
        legend.direction = "horizontal")

ggsave(filename = "Figures/barplot_Refs.png", dpi = 600)

ggplot(data = datum, aes(num_days)) + 
  geom_histogram(binwidth = 5)+ 
  facet_wrap(~refs_used)+
  scale_x_continuous(name = "Number of days in study area", breaks = seq(0, 120, 10)) +
  scale_y_continuous(name = "Count", breaks = seq(0, 100, 5))+
  # scale_fill_brewer(palette = "Set2", name = "Number of seasons") +
  theme_classic() + 
  theme(text = element_text(size = 12), 
        legend.position = c(.8,.9), 
        legend.direction = "horizontal")

ggsave(filename = "Figures/hist_days_v_Refs.png", dpi = 600)
