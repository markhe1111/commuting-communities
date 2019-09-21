

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




Countyshp <- counties(cb=TRUE) %>% st_as_sf() %>%
    mutate(FIPS = paste(STATEFP, COUNTYFP, sep="")) %>%
    filter(!(STATEFP %in% c("02", "15", "72", '66', '60', '78', '69'))) %>%
    st_transform("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

#us.citiesg.5M <- us.cities[us.cities$pop >500000,] %>%  
#   st_as_sf(coords = c('long', 'lat'),crs = "+init=epsg:4326") %>%
#   st_transform("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

msa2010 <- core_based_statistical_areas(year='2010', cb=TRUE) %>% 
    st_as_sf() %>% 
    filter(LSAD == "Metro") %>%
    st_transform("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")


divisions2010 <- st_read("/Users/markhe/Dropbox/regiondemarcation/rawdata/cb_2017_us_division_500k.shp") %>% 
  st_transform("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

megaregions2010 <- st_read('/Users/markhe/Dropbox/regiondemarcation/rawdata/mega_region_smooth.shp') %>%
  st_transform("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

  

comms_cvt_shp <- map(comms_cvt, function(x){Countyshp %>% filter(FIPS %in% x) %>% st_union() %>% st_sf})
comms_cvt_shp <- do.call(what = sf:::rbind.sf,args = comms_cvt_shp)
comms_cvt_shp$ID <- 1:nrow(comms_cvt_shp) %>% as.character()

hubspoke_cvt_shp <- purrr::map(hubspoke_cvt, function(x){Countyshp %>% filter(FIPS %in% x) %>% st_union() %>% st_sf})
hubspoke_cvt_shp <- do.call(what = sf:::rbind.sf,args = hubspoke_cvt_shp)
hubspoke_cvt_shp$ID <- 1:nrow(hubspoke_cvt_shp) %>% as.character()

island_cvt_shp <- Countyshp %>% filter(FIPS %in% unlist(island_cvt))
island_cvt_shp$ID <- 1:nrow(island_cvt_shp) %>% as.character()

tmap_options(max.categories = 78)

comms_map <- tm_shape(Countyshp) + tm_borders(lwd=0.3, col='gray75', alpha=.5)+
  tm_shape(comms_cvt_shp) + tm_polygons(col='ID', alpha=.6, lwd=1, legend.show = FALSE) + 
  tm_layout(frame = FALSE)

hubs_map <- tm_shape(Countyshp) + tm_borders(lwd=0.3, col='gray75', alpha=.5)+
  tm_shape(hubspoke_cvt_shp) + tm_polygons(col='ID', alpha=.6, lwd=1, legend.show = FALSE) + 
  tm_layout(frame = FALSE)

islands_map <- tm_shape(Countyshp) + tm_borders(lwd=0.3, col='gray75', alpha=.5)+
  tm_shape(island_cvt_shp) + tm_polygons(col='red', alpha=.6, lwd=1, legend.show = FALSE) + 
  tm_layout(frame = FALSE)
 # tm_shape(msa2010) + tm_borders(lwd=2, lty = 4) 

msa_divison_megaregion_map <- tm_shape(Countyshp) + tm_borders(lwd=0.3, col='gray75', alpha=.5)+
  tm_shape(msa2010) + tm_polygons(col='red', alpha=.6, lwd=1) + 
  tm_shape(megaregions2010) + tm_polygons(col='green', alpha=.6, lwd=1) + 
  tm_shape(divisions2010) + tm_borders(lwd=2)+
  tm_layout(frame = FALSE) + 
  tm_add_legend(type="fill",  col=c('red', 'green'), labels=c("MSA", 'Megaregions'))+
  tm_add_legend(type="line",  lwd=2.5, labels=c("Census Divisions"))



allmaps <- tmap_arrange(comms_map, hubs_map, islands_map, msa_divison_megaregion_map, ncol =2)

tmap_save(allmaps, filename = '/Users/markhe/Dropbox/regiondemarcation/manuscript/figs/Results_allregions.jpeg')




cropped_maps<-function(){
  
  chicago_canceled<-function(){
    #midwest = c(17, 18, 55, 26)
    midwest = c(17, 18, 55)
    View(key)
    
    # get ones w chicago 
    
    Countyshp[Countyshp$NAME=='Cook',]
    chicago = '17031'
    comms_w_chicago = Filter(function(x) is.element( chicago,x),  comms_cvt)
    
    cook_comms_cvt_shp <- map(comms_w_chicago, function(x){Countyshp %>% filter(FIPS %in% x) %>% st_union() %>% st_sf})
    cook_comms_cvt_shp <- do.call(what = sf:::rbind.sf,args = cook_comms_cvt_shp)
    cook_comms_cvt_shp$ID <- 1:nrow(cook_comms_cvt_shp) %>% as.character()
    
    
    cty_MW = Countyshp[Countyshp$STATEFP %in% midwest ,]
    
    MW_map  <- tm_shape(cty_MW) + tm_borders(lwd=0.4, col='gray75', alpha=.5)+
      tm_shape(cook_comms_cvt_shp) + tm_polygons(col='ID', alpha=.4, lwd=2, legend.show = FALSE) + 
      tm_layout(frame = FALSE)
    MW_map
  }

  ##-----
  
  Countyshp[Countyshp$NAME=='Bexar',]
  
  dallas = '48113'
  harris = '48201' 
  travis = '48453'
  bexar = '48029'
  
  chicago = '17031'
  
  SELcty = '36061'
    
  SEL = c(34,36)
  
  sf_oak_metro = 169
  
  
  
  
  
  
  dir.create( '/Users/markhe/Dropbox/regiondemarcation/manuscript/figs/examples/')
  setwd( '/Users/markhe/Dropbox/regiondemarcation/manuscript/figs/examples/')
  
  
  
  comm_mapper<-function(SELcty, SEL, all=F){
    
    
    shapes_SEL = Countyshp[Countyshp$STATEFP %in% SEL ,]
    
    if(all==F){
      comms_w_SEL = Filter(function(x) is.element( SELcty ,x),  comms_cvt)
      
      listcomms = lapply( 1:length(SELcty) , function(i)  Filter(function(x)   SELcty[i] %in%x ,  comms_cvt) )   
      
      comms_w_SEL = unlist(listcomms, recursive =F)
      
      SEL_comms_cvt_shp <- map(comms_w_SEL, function(x){Countyshp %>% filter(FIPS %in% x) %>% st_union() %>% st_sf})
      SEL_comms_cvt_shp <- do.call(what = sf:::rbind.sf,args = SEL_comms_cvt_shp)
      SEL_comms_cvt_shp$ID <- 1:nrow(SEL_comms_cvt_shp) %>% as.character()
      
    }else{
      
      SEL_comms_cvt_shp <- map(comms_cvt, function(x){Countyshp %>% filter(FIPS %in% x) %>% st_union() %>% st_sf})
      SEL_comms_cvt_shp <- do.call(what = sf:::rbind.sf,args = SEL_comms_cvt_shp)
      SEL_comms_cvt_shp$ID <- 1:nrow(SEL_comms_cvt_shp) %>% as.character()
      
    }
    
    SEL_map  <- tm_shape(shapes_SEL) + tm_borders(lwd=1, col='gray75', alpha=.7)+
      tm_shape(SEL_comms_cvt_shp) + tm_polygons(col='ID', alpha=.4, lwd=1, legend.show = FALSE) + 
      tm_layout(frame = FALSE)
    SEL_map
  }
  msa_mapper<-function( msa2010, SEL, msa ){
 
    shapes_SEL = Countyshp[Countyshp$STATEFP %in% SEL ,]
    
    SEL_map  <- tm_shape(shapes_SEL) + tm_borders(lwd=1, col='gray75', alpha=.7)+
      tm_shape(msa2010[msa,]) + tm_polygons(col='ID', alpha=.4, lwd=1, legend.show = FALSE) + 
      tm_layout(frame = FALSE)
    SEL_map
  }
    
     
  ## testing:: 
  comm_mapper( '36061',  c(34,36) )
  msa_mapper(  msa2010, SEL=c(34,36), msa = 259)
  
  
  wake = comm_mapper('37183', c(37,45,51,24) )
  
  tmap_save(wake, 
            filename = 'NC_example.pdf')
  
  
  
  
  chicago_comm = comm_mapper('17031', c(17,18))
  chicago_msa = msa_mapper(  msa2010, SEL=c(17,18), msa = 101)
  
  ny_comm = comm_mapper( '36061',  c(34,36) )
  ny_msa = msa_mapper(  msa2010, SEL=c(34,36), msa = 259)
  
  
  tx_msas = c(13, 287, 349, 167)
  tx_fips = c('48201','48029','48113','48453')
  
  tx_comm = comm_mapper(  SELcty= tx_fips ,  SEL=48, all = F)
  tx_msa = msa_mapper(  msa2010, SEL=48 , msa =  tx_msas)
  
  minn_comm = comm_mapper(  SELcty= '27053',  SEL=c(27), all = F)
  minn_msa = msa_mapper(  msa2010, SEL=c(27) , msa =  58)
  
  
  ###########
  sf_comm = comm_mapper(  SELcty=  '06075' ,  SEL='06', all = F)
  sf_msa = msa_mapper(  msa2010, SEL= '06' , msa =  169)
  
  example_maps <- tmap_arrange(   ny_comm, ny_msa, 
                                tx_comm, tx_msa,
                                minn_comm, minn_msa, ncol =2)
  
  example_maps_3 <- tmap_arrange(     tx_comm, tx_msa,
    chicago_comm, chicago_msa, 
                               minn_comm, minn_msa, ncol =2)
  
  sf_maps <- tmap_arrange(sf_comm, sf_msa , ncol =2)
  
  
  #allmaps <- tmap_arrange(comms_map, hubs_map, islands_map, msa_divison_megaregion_map, ncol =2)
  tmap_save(example_maps, 
            filename = '/Users/markhe/Dropbox/regiondemarcation/manuscript/figs/comms_case_examples.pdf')
  
   
}


tmap_save(allmaps, filename = '/Users/kaza/Dropbox/regiondemarcation/manuscript/figs/Results_allregions.pdf')

regioncounts <- c(comms_cvt %>% unlist(), hubspoke_cvt %>%unlist(), island_cvt %>% unlist()) %>% 
  enframe(name = NULL) %>%
  rename(FIPS = value) %>%
  group_by(FIPS) %>% summarize(count = n())

regioncountsmap <- Countyshp %>% left_join(regioncounts) %>%
  tm_shape() +
  tm_polygons('count', palette = 'YlOrRd', style="cat", n=6, colorNA =NULL, textNA = "Not in a Region", title = "Number of Regions", border.alpha = .3, lwd=.5, legend.is.portrait = FALSE ) +
  tm_layout(frame = FALSE)

tmap_save(regioncountsmap, filename = '/Users/kaza/Dropbox/regiondemarcation/manuscript/figs/Results_regioncount.pdf')

    

msa_why<-function(){
  msa_cbsa = core_based_statistical_areas(year='2010', cb=TRUE)   %>%  st_as_sf() %>% filter(LSAD == "Metro") 
  msa_cbsa = core_based_statistical_areas(year='2010', cb=TRUE)  
  df_msa_cbsa = msa_cbsa@data
  hmm = data(msa_cbsa)
  
  
  micro2010 <- core_based_statistical_areas(year='2010', cb=TRUE) %>% 
    st_as_sf() %>% 
    filter(LSAD == "Metro") %>%
    st_transform("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  
  
  
  
}
    

