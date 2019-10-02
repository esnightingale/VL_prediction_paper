
################################################################################
# Mapping predictions
################################################################################
library(grid)
library(cartography)

shp_df <- fortify(shp, region = "OBJECTID")
shp_df <- merge(agg1325, shp_df, by.x = "OBJECTID", by.y = "id")
## Reorder the new data file to prevent tearing the polygons
shp_df <- shp_df[order(shp_df$OBJECTID,shp_df$Year), ]


dist <- readOGR(dsn="C:/Users/phpuenig/Dropbox/VL/RiskMapper/example shapefiles/india", 
                layer="test_IND_district")
dist@data <- setNames(data.frame(dist, do.call(rbind, strsplit(as.character(dist$dot_name),split=":"))),c("hid_id","dot_name","continent","country","state","district"))
dist@data$vl <- 0
dist@data$vl[(dist$state=="BIHAR")|
               (dist$state=="JHARKHAND"&(dist$district %in% c("DUMKA","GODDA","PAKAUR","SAHIBGANJ")))] <- 1
# (dist$state=="WEST BENGAL" & !(dist$district%in%c("BANKURA","HOWRAH","JALPAIGURI","KOCH BIHAR","KOLKATA","PASCHIM MEDINIPUR", "PURBA MEDINIPUR","PURULIA")))] <- 1

pal <- cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
india <- tm_shape(dist) +
  tm_fill(col = "vl", palette = pal[c(4,2)]) + 
  # tm_borders(col=) +
  tm_layout(legend.show = F) 
#tm_compass(type = "4star", position = c("right","top")) +
#tm_scale_bar(breaks = c(0, 250,500), position = c("left","bottom")) 
india

vl <- tm_shape(shp) + 
  #  tm_compass(type = "8star", position = c("left", "bottom")) +
  tm_scale_bar(breaks = c(0, 25, 50), position = c("left","bottom")) +
  tm_fill(col = "inc", 
          title = "Incidence per 10,000", 
          style = "jenks",
          palette = "plasma") +
  tm_layout(frame = F, legend.outside = T) +
  tm_facets(along = "Year")

vl
print(india, vp = viewport(0.8, 0.27, width = 0.4, height = 0.4))
tmap_animation(vl, filename = "vl_anim.gif", delay = 25)

library(ggplot2)
library(maps)
library(ggthemes)
library(gganimate)
vlmap <-ggplot() + 
  geom_polygon(data = shp, 
               aes(x = long, y = lat, group = group), 
               colour = "black", fill = NA) +
  theme_map()
vlmap

shp_df <- broom::tidy(shp, region = "OBJECTID")
head(shp_df)

vlmap +
  geom_point(aes(x = lon, y = lat, fill = inc, 
                 frame = "Year"))


mymap <- ggplot(shp_df, aes(x = long, y = lat, group = group, fill = inc)) +
  coord_map() + 
  geom_polygon(color = "grey") +  
  transition_states(shp_df$Year, transition_length = 7, state_length = 7, 
                    wrap = TRUE) + 
  labs(title = "Annual predicted incidence per 10,000 population",
       subtitle = "Source: KAMIS + Surveillance modelling")
gganimate::animate(mymap)

# tmap_animation(tm, filename = "animation.gif", width = NA,
#                height = NA, dpi = NA, delay = 40, loop = TRUE,
#                restart.delay = 0)
#tm_style("col_blind")

for(y in 2013:2025){
  print(mapplot(shp,agg1325[agg1325$Year==y,],Cases,"Median"))
}

map_inc <- mapplot(shp,agg1325,inc,"Incidence")
gganimate(map_inc)
gganimate(map_inc, filename = 'pred_inc.gif')

map_median + transition_time(Year) +
  labs(title = "Year: {frame_time}")


# Surveillance animate fcn works:
animate(stsobj,
        timeplot = list(as.Date = TRUE,
                        scales = list(x = list(format = "%G/%V"))))

# Need to add predictions into object
cases2 <- rbind(cases, mdfc_blk)
cases2[is.na(cases2)] <- 0
popfrac2 <- pops/rowSums(pops)
stsobj_pred <- sts(observed= cases2, start = start.sts, frequency = freq, neighbourhood = nbOrd, map = shp, population = pops)
saveHTML(animate(stsobj_pred,
                 timeplot = list(as.Date = TRUE,
                                 scales = list(x = list(format = "%G/%V")))))



plot.map <- ggplot(data = shp.df, aes(long, lat, group = group)) 
plot.map <- plot.map + geom_polygon(aes(fill = value)) +
  plot.map <- plot.map + geom_path(colour = 'black')
plot.map <- plot.map + scale_fill_gradient2(name = legend_title,high = cols[1],mid = cols[2],low = cols[3]) 
plot.map <- plot.map + coord_equal() #+ theme_map() #+ geom_text(data=shp.df, aes(x=long, y=lat, label=Block, group=Block), size=0.5)

plot.map + 
  transition_states(Year,
                    transition_length = 2,
                    state_length = 1)

p <- ggplot(agg1325, aes(x = Cases, y = pop, group=OBJECTID )) + 
  geom_point() + theme(legend.position = "none")
p
p + 
  transition_states(Year,
                    transition_length = 2,
                    state_length = 1)



