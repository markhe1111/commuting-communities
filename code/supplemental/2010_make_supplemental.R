

#install.packages('tigris')
#install.packages('sf')
#install.packages('tidyverse')
#install.packages('tmap')


library(tigris)
library(sf)
library(tidyverse)
library(tmap)


load('/Users/markhe 1/Dropbox/regiondemarcation/code/git/community-analysis/Postprocessing/Results/communities_and_hubs_cvt2.Rdata')
load('/Users/markhe 1/Dropbox/regiondemarcation//results/community_results/comms_supplement.Rdata')


filtered_50$final_comms

Countyshp <- counties(cb=TRUE) %>% st_as_sf() %>%
    mutate(FIPS = paste(STATEFP, COUNTYFP, sep="")) %>%
    filter(!(STATEFP %in% c("02", "15", "72", '66', '60', '78', '69'))) %>%
    st_transform("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

#us.citiesg.5M <- us.cities[us.cities$pop >500000,] %>%  
#   st_as_sf(coords = c('long', 'lat'),crs = "+init=epsg:4326") %>%
#   st_transform("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
library(htmlwidgets)

library(tnet)
rm(list = ls())
#load('/Users/markhe/Dropbox/regiondemarcation//results/community_results/comm_results.Rdata')

#load('/Users/markhe/Dropbox/regiondemarcation//results/community_results/comms_2_20.Rdata')

load('/Users/markhe 1/Dropbox/regiondemarcation/code/git/community-analysis/Postprocessing/Results/communities_and_monocentric.Rdata')


load('/Users/markhe 1/Dropbox/regiondemarcation/results/community_results/comms_supplement.Rdata')

#load('/Users/markhe/Dropbox 1/regiondemarcation/community_results/comms_supplement.Rdata')


cvt  <-function(cluslist, key) { 
  
  cluslist.char = lapply(cluslist, function(x) as.character(key[x,3]))
  
  as.char.fill <- function(x)  if(nchar(x)==4) paste('0', x,sep='')  else x 
  char.fil.clus <-function(comm)  unname(sapply( comm,  function(i) (as.char.fill(i) ) ))
  
  gp.list = lapply(cluslist.char, function(x) char.fil.clus(x))
  gp.list
} 


library(leaflet)
library(magrittr)
library(htmltools)
library(tnet)

 
comms = cvt(filtered$final_comms,key )
comms50 = cvt(filtered_50$final_comms, key)

conserv_comms = cvt( filtered_conserv$final_comms, key)
conserv_comms50 = cvt(filtered_conserv_50$final_comms , key)

permissive_comms = cvt(filtered_permissive$final_comms , key)
permissive_comms50 = cvt( filtered_permissive_50$final_comms , key)



load('/Users/markhe 1/Dropbox/regiondemarcation/code/git/community-analysis/Postprocessing/Results/communities_and_hubs_cvt2.Rdata')
Countyshp <- counties(cb=TRUE) %>% st_as_sf() %>%
  mutate(FIPS = paste(STATEFP, COUNTYFP, sep="")) %>%
  filter(!(STATEFP %in% c("02", "15", "72", '66', '60', '78', '69'))) %>%
  st_transform("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

anycomm_shp <-function(COMM){
  C_Cvt_shp <- map(COMM, function(x){Countyshp %>% filter(FIPS %in% x) %>% st_union() %>% st_sf})
  C_Cvt_shp <- do.call(what = sf:::rbind.sf,args = C_Cvt_shp)
  C_Cvt_shp$ID <- 1:nrow(C_Cvt_shp) %>% as.character()
  comms_map <- tm_shape(Countyshp) + tm_borders(lwd=0.3, col='gray75', alpha=.5)+
    tm_shape(C_Cvt_shp) + tm_polygons(col='ID', alpha=.6, lwd=1, legend.show = FALSE) + 
    tm_layout(frame = FALSE)
  comms_map
  
}


comm_map = anycomm_shp(   comms50 )
conserv_map = anycomm_shp(   conserv_comms )
monad_map = anycomm_shp(   my_singletons )

my_tot = c(my_comms, my_hubs, my_singletons)

