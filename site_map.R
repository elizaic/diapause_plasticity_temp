library(readxl)
library(tidyverse)
library(maps)
library(maptools)
library(sf)
library(rnaturalearth)
library("rnaturalearthdata")
library(rnaturalearthhires)
library(viridis)
library(cowplot)
library(mapdata)


# ##read in collection site data
collections <- read_xlsx("2022collection_site_coords.xlsx")
#
# ##range expansion points
# expansion <- read_xlsx("collection_sites.xlsx",
#                        sheet = "Sheet2")
# expansion <- as.data.frame(expansion)
# expansion$site <- as.factor(expansion$site)
# expansion$year <- as.factor(expansion$year)
# expansion <- expansion %>% arrange(year)
# expansionSubset <- expansion %>% filter(year == 2010 | year == 2013 | year == 2016)

##river data
shapeData <- rgdal::readOGR("ne_10m_rivers_lake_centerlines",
                            layer = "ne_10m_rivers_lake_centerlines")
shapeData@data$id <- rownames(shapeData@data)
watershedPoints <- fortify(shapeData, region = "id")
watershedDF <- merge(watershedPoints, shapeData@data, by = "id")

### Focusing on N America river data
NA_shapeData <- rgdal::readOGR("ne_10m_rivers_north_america",
                               "ne_10m_rivers_north_america")
NA_shapeData@data$id <- rownames(NA_shapeData@data)
NA_watershedPoints <- fortify(NA_shapeData, region = "id")
NA_watershedDF <- merge(NA_watershedPoints, NA_shapeData@data, by = "id")

NA_watershedDF$scalerank <- as.numeric(NA_watershedDF$scalerank)

NAwater <- NA_watershedDF %>%
  dplyr::filter(long < -95 & long >-120 & lat > 25 & lat <44.8)
NA_water <- NA_watershedDF %>%
  dplyr::filter(long < -95 & long >-120 & lat > 25 & lat <44.8)

#land data
landData <- rgdal::readOGR("ne_10m_land",
                           layer = "ne_10m_land")
landData@data$id <- rownames(landData@data)
land_watershedPoints <- fortify(landData, region = "id")
land_watershedDF <- merge(land_watershedPoints, landData@data, by = "id")

##state outlines
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
# states$ID <- toTitleCase(states$ID)
world <- ne_countries(scale = "medium", returnclass = "sf")
state_prov <- rnaturalearth::ne_states(c("united states of america", "canada"))

##colors for collection sites
edge_col = viridis(3)[2]
core_col = viridis(3)[1]

map <- ggplot() +
  theme_bw(base_size = 20) +
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.grid = element_blank()) +
  # geom_rect(xmin = -127, xmax=-65, ymin=25, ymax = 50,
  #           fill = 'grey50', color = 'black', size = 2) +
  geom_sf(data = world, fill = "white") + #white land
  # geom_sf(data = states, fill = 'white') + # state boundaries
  geom_path(data=watershedDF, # rivers
            aes(x = long, y = lat, group = group),
            color = 'lightblue', size=as.numeric(watershedDF$scalerank*.1)) +
  geom_path(data=NA_water, # more rivers
            aes(x = long, y = lat, group = group),
            color = 'lightblue', size=as.numeric(NA_water$scalerank*.06)) +

  # geom_point(data = expansion, #range expansion points
  #            aes(x = longitude, y = latitude, fill = year), shape = 22, size = 2) +
  # scale_fill_brewer(palette = "Reds", name = "Range edge") +
  geom_point(data = collections, #collection points
             aes(x = longitude, y = latitude), size = 4) +
  # scale_color_manual(values = c(edge_col, core_col), name = "Collection sites \n2017-2018") +
  # scale_shape_manual(values = c(16, 17), name = "Collection sites \n2017-2018") +
  geom_text(data = collections, aes(x = longitude, y  = latitude, label = label),
            position = position_nudge(x=0, y = .38), hjust = 'center', size = 5) +
  # guides(shape = guide_legend(order = 1), col = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-118, -110), ylim = c(30, 42), expand = FALSE)+  #crop to NA
  scale_x_continuous(breaks = seq(-118, -110, by = 4))
map

inset <- ggplot(data = world) +
  theme_void() +
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.grid = element_blank()) +
  # geom_rect(xmin = -127, xmax=-65, ymin=25, ymax = 50,
  #           fill = 'lightgrey', color = NA) +
  geom_sf(fill = "white") + #white land
  geom_sf(data = states, fill = NA, size = 0.25) + # state boundaries
  coord_sf(xlim = c(-127, -65), ylim = c(25, 50), expand = FALSE) +
  # theme(legend.position = "none",
  #       axis.text=element_blank(),
  #       axis.ticks=element_blank(),
  #       panel.grid.major=element_blank()) +
  geom_rect(xmin = -118, xmax=-110, ymin=30, ymax = 42,
            fill = NA, color = 'red', size = 1) +
  geom_rect(xmin = -127, xmax=-65, ymin=25, ymax = 50,
            fill = 'NA', color = 'black', size = 1)
inset

full_map <- map %>%
  ggdraw() +
  draw_plot(inset,
            x = 0.38, y = .7, width = .4, height = .4)
full_map
##export: width = 450 x 596 px

# full_map %>%
#   ggsave(filename = "Figure_2.pdf", units = 'mm', width = 180, height = , device='pdf', dpi=1000)


##export width = 600, maintain aspect ratio


