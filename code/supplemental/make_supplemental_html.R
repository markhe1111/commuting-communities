library(htmlwidgets)

library(tnet)
rm(list = ls())
#load('/Users/markhe/Dropbox/regiondemarcation//results/community_results/comm_results.Rdata')

load('/Users/markhe/Dropbox/regiondemarcation//results/community_results/comms_2_20.Rdata')

load('/Users/markhe/Dropbox/regiondemarcation/code/git/community-analysis/Postprocessing/Results/communities_and_monocentric.Rdata')
load('/Users/markhe/Dropbox/regiondemarcation//results/community_results/comms_supplement.Rdata')



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
  
  #user.directory = "c:/Users/spamb/" #where your dropbox file is 
  user.directory = "/Users/markhe/"
  setwd(paste(user.directory, "Dropbox/regiondemarcation/code/git/community-analysis/data", sep = ''))
  
  load("continental_counties.Rdata") #loads continental.counties which is a dataframe necessary for plotting counties in U.S
  load("2010_undirected_edgelist.Rdata")
  
  
  CreateMap = function(group.list, county = NULL, state = NULL, 
                       named.counties = NULL, named.county.labels = NULL, heat.map = F){
    # Creates interactable leaflet map for communities. More on leaflet maps can be found here: https://rstudio.github.io/leaflet/
    #
    # Args:
    #   group.list: list of groups of counties. Each county is represented by a string that is the concatanation of its state and county fips codes.
    #               each group represents a community
    #
    #   state: fips code of a state, if not null will only plot counties belonging to that state
    #
    #   named.counties: list of county fips to be plotted with permanent labels (they persist despite not having a cursor over them)
    #
    #   named.county.labels: list of labels to be included for each county in named.counties. named.county.labels[[i]] should be a list of info (as strings) to be
    #                        included for named.counties[[i]]. 
    #
    #   heat.map: if true will make the map a heatmap indicating the number of communities each county belongs to
    #   
    #
    # Returns:
    #   A leaflet map 
    
    
    #Helper function to extract counties from a single state 
    substate.cover <- function(cover, state.fips){
      only.state.counties <- sapply(cover, function(y) y[sapply(y, function(x)  as.numeric(substr(x,1,2)) == state.fips)])
      only.state.counties <- only.state.counties[sapply(only.state.counties, length) > 0]
      return(only.state.counties)
    }
    
    #Modifying cover and map if only a single state is selected 
    if(!is.null(state)){
      group.list <- substate.cover(group.list, state)
      continental.counties <- continental.counties[continental.counties$STATEFP == state,]
    }
    
    
    #Filtering out counties from group.list not found in continental.counties which is used for plotting
    counties.not.found = list()
    list.index = 1
    
    
    
    for(i in 1:length(group.list)){
      counties = as.numeric(group.list[[i]])
      counties.in.shapefile = counties %in% continental.counties$fips
      if(sum(!counties.in.shapefile) != 0){
        missing.counties = counties[!counties.in.shapefile] #counties not found in continental.counties
        missing.counties = missing.counties[!(missing.counties %in% counties.not.found)] #counties not found that haven't been accounted for 
        if(length(missing.counties) > 0){ #adding unaccounted for counties to counties.not.found list
          for(j in 1:length(missing.counties)){
            counties.not.found[[list.index]] = missing.counties[j]
            list.index = list.index + 1
          }
        }
        group.list[[i]] = as.character(counties[counties.in.shapefile]) #filtering out counties from group.list 
      }
    }
    
    group.list = group.list[sapply(group.list,length) > 0] #filtering out empty communities 
    
    if(length(counties.not.found > 0)){ #if any counties were not found 
      print("Counties in the group list were not found in the shape file and will not be plotted.")
      print(paste("Counties not found: ", paste(counties.not.found, collapse = ", "), sep = ""))
    }
    
    
    #orders group list by number of counties 
    ordered.group.list = group.list[order(sapply(group.list,length), decreasing = TRUE)]
    
    #Subsetting data for specific county 
    if(!is.null(county)){
      boolean.subset = as.numeric(lapply(ordered.group.list, function(x) county %in% as.numeric(x))) == 1
      ordered.group.list = subset(ordered.group.list,boolean.subset)
    }
    
    
    #creates color palette for mapping 
    colors = sample(rainbow(length(ordered.group.list)))
    
    
    #community.key is a dataframe that will be used to determine which groups belong to which layers 
    #first column is the community id 
    #second column is used to keep track of which community as already been accounted for in a layer
    #subsequent columns indicate if the community will be included in layer 1, 2, 3, etc. 
    community.key = data.frame(1:length(ordered.group.list))
    colnames(community.key) = c("id")
    community.key$plotted = FALSE 
    
    
    
    #map.layers is a dataframe that will be used to find which counties belong in what layers as well as their color 
    #the first column is the fips codes for each county
    #columns labeled "includeX" will indicate whether each county will be included in layer X
    #columns labeled "colorX" will give the color for the county in layer X 
    map.layers = data.frame(continental.counties$fips)
    
    
    #a method to obtain the column name for map.layers
    GetLabel = function(x, color = FALSE){
      if(color){
        return(paste("color",as.character(x),sep = ""))
      } else {
        return(paste("include",as.character(x),sep = ""))
      }
    }
    
    
    #creating disjoint layers of communities 
    column = 3
    while(TRUE){
      #initializing columns of community.key and map.layers 
      community.key[,column] = FALSE
      map.layers[,GetLabel(column-2)] = FALSE
      map.layers[,GetLabel(column-2, TRUE)] = "transparent"
      
      for(i in community.key[!(community.key$plotted),"id"]){ #for the community id's that have not been plotted
        if(sum(ordered.group.list[[i]] %in% unlist(ordered.group.list[community.key[,column]])) == 0){ #if this community id does not overlap with any others in the layer 
          community.key[i,column] = TRUE
          community.key[i, "plotted"] = TRUE
          counties = as.numeric(ordered.group.list[[i]])
          county.bool = map.layers[,1] %in% counties #a boolean to select rows of counties that are in the current community 
          map.layers[county.bool, GetLabel(column-2)] = TRUE
          map.layers[county.bool, GetLabel(column-2, TRUE)] = substr(colors[i],1,7)
        }
        
      }
      if(sum(community.key$plotted) == nrow(community.key)){ #if every community has been accounted for 
        break
      } else {
        column = column + 1 #otherwise, we add a new layer 
      }
    }
    
    
    
    #tags will contain the data for the labels of the counties 
    tags = data.frame(continental.counties$NAME, continental.counties$fips)
    tags[,3] = ""
    colnames(tags) = c("name","fips","communities")
    tags$community.count = 0
    
    for(i in 1:length(ordered.group.list)){ #interating over each community
      counties = as.numeric(ordered.group.list[[i]]) 
      for(j in 1:length(counties)){ #adding group number to each county tag 
        county = counties[j]
        if(tags[tags$fips == county,"communities"] == ""){
          tags[tags$fips == county,"communities"] = as.character(i)
          tags[tags$fips == county,"community.count"] = 1
        } else {
          tags[tags$fips == county,"communities"] = paste(tags[tags$fips == county,"communities"], as.character(i), sep = ", ")
          tags[tags$fips == county,"community.count"] = tags[tags$fips == county,"community.count"] + 1
        }
      }
    }
    
    #creating labels in html format from the dataframe tags 
    labels <- sprintf(
      "<strong>%s</strong><br/>Community: %s",
      tags$name, tags$communities
    ) %>% lapply(htmltools::HTML)
    
    
    #initializing leaflet map
    map <- leaflet(continental.counties) %>% addPolygons(color = "#000000", fillColor = "#FFFFFF", fillOpacity = 1, weight = 1 ,group = "main",layerId = continental.counties$fips)
    
    
    
    #making heatmap 
    if(heat.map){
      
      pal = colorBin("YlOrRd", domain = tags[tags$community.count > 0,"community.count"], bins = c(1,2,3,4,5,6,7, 8))
      map = map %>% addPolygons(data = continental.counties[tags$community.count > 0,], color = "#000000", fillColor = ~pal(tags[tags$community.count > 0,"community.count"]), fillOpacity = 1,
                                weight = 1, group = "heat map") %>%
        addLegend("bottomright", pal = pal, values = ~tags[tags$community.count > 0,"community.count"], opacity = 1, title = "Community Count")
      
      
    } else {
      #adding layers if heatmap option isn't selected 
      for(i in 1:as.integer(ncol(map.layers)/2)){
        if(!is.null(county)){
          county.rows = map.layers[,GetLabel(i)]
        } else {
          county.rows = map.layers[map.layers[,1] == county ,GetLabel(i)]
        }
        map = map %>% addPolygons(data = continental.counties[county.rows,], color = "#000000", fillColor = map.layers[county.rows, GetLabel(i,TRUE)], fillOpacity = 0.5,
                                  weight = 1, group = as.character(i))
        
      }
    }
    
    
    map = map %>% addPolygons(data = continental.counties, color = "transparent", fillColor = "transparent", fillOpacity = 1, weight = 1 ,group = "labels",
                              label = labels)
    
    permanent.labels.included = F
    
    #permanent labels
    if(!is.null(named.counties)){
      if(length(named.counties) != length(named.county.labels)){
        print("Length of named.county.labels does not equal the length of named.counties. No permanent labels will be included")
      } else {
        has.permanent.label = continental.counties$fips %in% sapply(named.counties, as.numeric)
        
        permanent.county.names = tags$name[has.permanent.label]
        
        permanent.county.labels = list()
        
        #Adding county info to each label
        for(i in 1:length(named.counties)){
          county.info = named.county.labels[[i]]
          county.label = sprintf("<strong>%s</strong>", permanent.county.names[[i]])
          
          for(j in 1:length(county.info)){
            info.line = county.info[[j]]
            county.label = paste(county.label, "<br/>", info.line, sep = "")
          }
          permanent.county.labels[i] = county.label
        }
        
        permanent.county.labels = lapply(permanent.county.labels, htmltools::HTML)
        
        map = map %>% addPolygons(data = continental.counties[has.permanent.label,], color = "transparent", fillColor = "transparent", fillOpacity = 1, weight = 1 ,group = "permanent labels",
                                  label = permanent.county.labels, labelOptions = labelOptions(noHide = T, textOnly = T))
        permanent.labels.included = T
      }
    }
    
    #adding layer control
    overlayGroups = c("labels")
    
    if(permanent.labels.included){
      overlayGroups = c(overlayGroups, "permanent labels")
    }
    
    if(heat.map){
      overlayGroups = c(overlayGroups, "heat map")
    } else {
      overlayGroups = c(overlayGroups, as.character(1:as.integer(ncol(map.layers)/2)))
    }
    
    map = map %>% addLayersControl(
      baseGroups = c("main"),
      overlayGroups = overlayGroups,
      options = layersControlOptions(collapsed = TRUE)
    )
    
    return(map)
  }
  
  
  
  comms = cvt(filtered$final_comms,key )
  comms50 = cvt(filtered_50$final_comms, key)
  
  conserv_comms = cvt( filtered_conserv$final_comms, key)
  conserv_comms50 = cvt(filtered_conserv_50$final_comms , key)
  
  permissive_comms = cvt(filtered_permissive$final_comms , key)
  permissive_comms50 = cvt( filtered_permissive_50$final_comms , key)
  
  
  CreateMap(comms)
  CreateMap(comms50)
  
  CreateMap(conserv_comms)
  CreateMap(conserv_comms50)
  
  CreateMap(permissive_comms)
  CreateMap(permissive_comms50)
  
  
  
  
  CreateMap(cvt(list(islands_05$node), key))
  CreateMap(cvt(list(islands_05_95$node), key))
  
  
  CreateMap(cvt(list(islands_01$node), key))
  CreateMap(cvt(list(islands_01_95$node), key))
  
  
  CreateMap(islands_01)
  CreateMap(islands_01_95)
   