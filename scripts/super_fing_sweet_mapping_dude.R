###############################################################################L


# Mapping and plotting for multi-state model


##############################################################################L

# Load necessary packages
library(tidyverse)
library(tidybayes)
library(sf)
library(tidygraph)
library(ggspatial)
library(ggraph)
library(ggrepel)
library(sp)
library(patchwork)

# Run script to load and prep data
source("scripts/DataPrep_MultistateModel.R")


# Load data
data("us_states", package = "spData")

ref_points <- cbind(Refs, st_coordinates(st_transform(st_centroid(Refs$geometry), 
                                                      crs = st_crs(us_states))))
bbox <- Polygon(matrix(c(-90.25, 35.7,
                         -88.65, 35.7,
                         -88.65, 37, 
                         -90.25, 37,
                         -90.25, 35.7),
                       ncol = 2, byrow = TRUE))
bbox <- SpatialPolygons(list(Polygons(list(bbox), ID = "a")), 
                        proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))
bbox <- st_as_sf(bbox, plot = FALSE, fill = "gray92", width = 4)

us.p <- ggplot() +
  geom_sf(data = us_states, fill = "white", col = "black")+
  coord_sf(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")+
  #geom_sf(data = Refs, aes(fill = Name, color = Name))+
  geom_sf(data = bbox, fill = "darkgray", size = 6, col = "black") +
  # scale_fill_brewer(palette = "Set3")+
  # scale_color_brewer(palette = "Set3")+
  # coord_sf(xlim = c(-90.2,-88.7), ylim = c(35.65,36.55), expand = FALSE)+
  # geom_text_repel(data = ref_points, aes(x = X, y = Y, label = Name),
  #           fontface = "bold", nudge_y = -.03, nudge_x = .03, min.segment.length = 1)
  theme_void()
  



probs_datum <- read.csv("Results/ref_probs.csv") %>% 
  filter(sex == 1, age == 1) %>% 
  rename(Probability = psi) %>%
  dplyr::select(-X) %>%
  mutate(lab.prob = round(Probability, 2)) %>%
  mutate(lab.prob = ifelse(lab.prob >0, lab.prob, ""))


prob_graph <- tbl_graph(st_transform(st_centroid(Refs),
                                     crs = st_crs(us_states)), 
                        probs_datum, node_key = "Name")
prob_coords <- create_layout(graph = prob_graph, 
                             layout = "manual", 
                             x = ref_points$X, 
                             y = ref_points$Y)
# ggplot() + 
#   # geom_sf(data = prob_graph %>% activate(edges) %>% as_tibble() %>% st_as_sf())+
#   geom_sf(data = prob_graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf())

label_datum <- data.frame(x = c(-90.1, -89.2, -89.9), 
                          y = c(35.7, 35.7, 36.2), 
                          lab = c("AR", "TN", "MO"))

refs.p <- ggraph(prob_coords) +   
  geom_sf(data = us_states, fill = "gray92", col = "black", size = 1)+
  geom_sf(data = st_transform(Refs, crs = st_crs(us_states)), color = "black", fill = "darkorange") +
  coord_sf(xlim = c(-90.2,-88.7), ylim = c(35.65,36.65), expand = TRUE)+
  geom_edge_diagonal(aes(width = Probability, alpha = Probability,
                         color=Probability))+
  geom_sf(data = st_transform(Refs, crs = st_crs(us_states)), fill = "darkorange", color = "black") +
  coord_sf(xlim = c(-90.2,-88.7), ylim = c(35.65,36.65), expand = TRUE)+
  # geom_edge_diagonal(aes(width = Probability, alpha = Probability,
  #                        color=Probability, label = lab.prob), angle_calc = 'along',
  #                    label_dodge = unit(2.5, 'mm'))+
  # geom_edge_loop(aes(width = Probability, alpha = Probability, color=Probability,
  #                    strength = .15)) +
  scale_edge_color_gradientn(colors = scales::brewer_pal("seq", palette = "Blues", 
                                                         direction = 1)(9)[6:9]) +
  scale_edge_alpha_continuous(range = c(0.01, 1), guide = "none") +
  geom_node_text(aes(label = Name, fontface = "bold"), repel = TRUE, size = 4) +
  scale_edge_width_continuous(range = c(1, 4), guide = "none") +
  geom_text(data = label_datum, aes(label = lab, x = x, y = y), size = 6.5, color = "black")+
  guides(fill = "none", width = "none") +
  theme_bw(base_size = 16) + 
  theme(axis.text.x = element_text(angle = 35, hjust = 0.5, vjust = 0.5)) +
  labs(x = "Longitude", y = "Latitude")
  #theme(legend.position = "top")

layout <- c(
  area(t = 2, l = 2, b = 8, r = 8),
  area(t = 1, l = 1, b = 4, r = 4))
plot(layout)

refs.p + us.p + plot_layout(design = layout)


# layout <- "
# ##BBBB
# AABBBB
# AABBBB"
# 
# us.p / refs.p
# us.p + refs.p + plot_layout(design = layout)
ggsave("Figures/study_area.png", width = 7, height = 6.5, dpi = 600)



