setwd("~/Desktop/PanPan/mapchimps/")

library(ggplot2)  # ggplot() fortify()
library(dplyr)  # %>% select() filter() bind_rows()
library(rgdal)  # readOGR() spTransform()
library(raster)  # intersect()
library(ggsn)  # north2() scalebar()
library(rworldmap)  # getMap()
library(spData)
library(tmap)
library(leaflet)
library(cartogram)
require(maptools)
library(maps)
library("ggspatial")

shpdata_ppaniscus<-readShapeLines(fn = "ppaniscus/data_0.shp")
shpdata_pt<- readShapeLines(fn = "pt/data_0.shp")
#extract features to a new polygon
shpdata_pte = shpdata_pt[1,]
shpdata_ptv = shpdata_pt[2,]
shpdata_pts = shpdata_pt[3,]
shpdata_ptt = shpdata_pt[4,]

world = world
Africa= world %>% filter(is.na(geom) | continent == 'Africa')
#4c5d4c -> central african (Pantrog trog)
#9dced9 -> western African (Pantrog verus)
#fe604c -> nigeria cameroon (Pantrog elioti)
#ffb35a -> eastern (Pantrog schweinfurthii)

myplot<-ggplot(Africa) + geom_sf(color = "#ABB0B8", fill = "#6F7378", size = 0.1)  + xlab("") + ylab("") + 
  geom_polygon(data = shpdata_ptt,aes(x = long, y = lat, group = group), fill = "#4c5d4c", color="#4c5d4c", alpha = 0.8, size = 0.3, show.legend = TRUE) +
  geom_polygon(data = shpdata_ptv,aes(x = long, y = lat, group = group), fill = "#9dced9", color="#9dced9", alpha = 0.8, size = 0.3, show.legend = TRUE) +
  geom_polygon(data = shpdata_pte,aes(x = long, y = lat, group = group), fill = "#fe604c", color="#fe604c", alpha = 0.8, size = 0.3, show.legend = TRUE) +
  geom_polygon(data = shpdata_pts,aes(x = long, y = lat, group = group), fill = "#ffb35a", color="#ffb35a", alpha = 0.8, size = 0.3, show.legend = TRUE) +
  geom_polygon(data = shpdata_ppaniscus,aes(x = long, y = lat, group = group), fill = "#e9d7cb", color="#e9d7cb", alpha = 1, size = 0.3, show.legend = TRUE) +
  theme_void() + theme(legend.position = "bottom") 
ggsave("myplot.tiff", myplot, width = 5, height = 5, dpi=600)

                   


