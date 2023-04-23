#' This is some description of this function.
#' @title to map data to the world or other locations
#'
#' @description By using this package, you could use functions of GeoMapping to map category data onto the world map or other type of maps, as long as longitude and latitude information was provided.
#'
#' @details see above
#'
#' @param haploResults: haplotype analysis results, could be either folder stores list of gene list files or a single file of gene list;
#' @param mapData: geographic data in 'shp' or 'json' format. If 'mapData' is NULL, map data would be downloaded from rnaturalearth (http://www.naturalearth.com).
#' @param type: could be either coastline or countries, default is countries.
#' @param region: the region to display on the map, could be "world", "China", or other countries, default is 'China'.
#' @param mapResolution: value of 10, 50 and 110, or string of 'small', 'medium', or 'large', default is 50.
#' @param longLim: limitation of longitude, default is c(-180.00, 180.00).
#' @param latLim: limitation of latitude, default is c(-90.00, 90.00).
#' @param landColor: the color of land, default is lightyellow.
#' @param oceanColor: the color of ocean, default is white.
#' @param pieLinetype: line type inherited from ggplot2, default is 2.
#' @param pieLinesize: line size inherited from ggplot2, default is 0.
#' @param label_graticule: parameter of coord_sf(), which is a character vector indicating which graticule lines should be labeled where. It could be 'N', 'S', 'W', 'E', or combination of these characters. Default is 'EW'.
#' @param pdfWidth: numeric value to indicate the width of pdf file, default is 15.
#' @param pdfHeight: numeric value to indicate the height of pdf file, default is 10.
#' @param pieMax: numeric value to indicate the maximum size of pie plot on geographic map, default is 2.
#' @param dataFiltration: logic value to indicatee whether to filtrate data by 'China', or other keywords, default is NULL.
#' @param pieBorderColor: border color of pie plot, default is NA.
#' @param pieAlpha: numeric value to define transparency of pie plot, default is 0.8.
#' @param expand: logic value to make the axes not overlap with each other, default is "T".
#' @param plotMargin: adjust graph edge margins;
#' @param MappingDataType: data types of geographical mapping, default is "haplotype";
#' @param mapBorderSize: border size of geographic mapping, default is 0.1.

#' @return files extracted and stored in the folder created by default
#' @export GeoMapping
#' @examples GeoMapping(haploResults = "./Haplotype.Analysis/",mapData = "./map.source/China.json")


###### NOTES ######
# there were many other parameters from geom_sf(), coord_sf(), geom_sf_label(), geom_sf_text(), stat_sf(), which could be introduced into the function for fine adjustment.
###################



GeoMapping = function(haploResults = NULL,
                      mapData = NULL,
                      type = "countries",
                      region = "China",
                      mapResolution = 50,
                      longLim = NULL,
                      latLim = NULL,
                      landColor = "lightyellow",
                      oceanColor = "white",
                      pieMax = 2,
                      pdfWidth = 15,
                      pdfHeight = 10,
                      expand = T,
                      pieLinetype = 2,
                      pieLinesize = 0,
                      pieBorderColor = NA,
                      pieAlpha = 0.8,
                      label_graticule = "EW",
                      plotMargin = 2,
                      MappingDataType = "haplotype",
                      mapBorderSize = 0.1,
                      dataFiltration = NULL, ...){
  require(reshape2)
  library(maps)
  require(rworldmap)
  require(rworldxtra)
  library(RColorBrewer)
  library(rnaturalearthdata)
  library(rnaturalearth)
  library(ggplot2)
  library(tidyverse)
  library(scatterpie)
  library(rgeos)
  library(data.table)
  library(pacman)
  library(grDevices)
  library(randomcoloR)
  library(dplyr)
  p_load(sf)

  if(is.null(haploResults)) stop("'haploResults' is required!\n")
  if(is.null(mapData)) cat("'mapData' is obtained from internet!\n")
  if(!is.null(mapData)) cat(paste0("'mapData' is reading from: ", mapData, "\n", "In this case, you have to check whether this map data is the one to be used in this project.\n"))

  dirPath = getwd()
  if(!dir.exists(paste0(dirPath, "/", "GeoMappingResults"))){
    dir.create(paste0(dirPath, "/", "GeoMappingResults"))
  }

  if(dir.exists(haploResults)){
    mapFileDir = paste0(haploResults, list.files(haploResults))
    genelist = c()
    for(i in 1:length(mapFileDir)){
      temp = strsplit(mapFileDir[i], "/")[[1]][length(strsplit(mapFileDir[i], "/")[[1]])]
      genelist = c(genelist, temp)
    }
  }else{
    mapFileDir = haploResults
    genelist = strsplit(mapFileDir, "/")[[1]][length(strsplit(mapFileDir, "/")[[1]])]
  }

  # map data
  if(is.null(mapData)){
    if(region == "world" | region == "World" | region == "WORLD"){
      if(type == "coast line" | type == "coastline" | type == "Coast Line" | type == "CoastLine"){
        map = ne_coastline(scale = mapResolution, returnclass = c("sf"))
      }
      if(type == "countries" | type == "Countries" | type == "COUNTRIES"){
        map = ne_countries(scale = mapResolution, type = "countries", continent = NULL,
                           country = NULL, geounit = NULL, sovereignty = NULL,
                           returnclass = "sf")
      }
    }
    if(region == "China" | region == "CHINA"){
      if(type == "coast line" | type == "coastline" | type == "Coast Line" | type == "CoastLine"){
        map = ne_coastline(scale = mapResolution, returnclass = c("sf"))
      }
      if(type == "countries" | type == "Countries" | type == "COUNTRIES"){
        map = ne_countries(scale = mapResolution, type = "countries", continent = NULL,
                           country = NULL, geounit = NULL, sovereignty = NULL,
                           returnclass = "sf")
      }
      map = map[which(map$sovereignt == "China" | map$sovereignt == "Taiwan"),]
    }
    if(region != "world" & region != "China"){
      if(type == "coast line" | type == "coastline" | type == "Coast Line" | type == "CoastLine"){
        map = ne_coastline(scale = mapResolution, returnclass = c("sf"))
      }
      if(type == "countries" | type == "Countries" | type == "COUNTRIES"){
        map = ne_countries(scale = mapResolution, type = "countries", continent = NULL,
                           country = NULL, geounit = NULL, sovereignty = NULL,
                           returnclass = "sf")
      }
      map = map[which(map$sovereignt == region),]
    }
  }
  if(!is.null(mapData)){
    map = st_read(mapData)
    map = map$geometry
  }

  ### longLim and latLim definition
  if(is.null(longLim)){
    if(region == "China" | region == "CHINA"){
      longLim = c(70, 140)
    }
    if(region == "World" | region == "world"){
      longLim = c(-180.00, 180.00)
    }
    if(length(grep("World|world|china|China", region)) != 1){
      stop("'longLim' is required if the 'region' to map is not China or World.")
    }
  }
  if(is.null(latLim)){
    if(region == "China" | region == "CHINA"){
      latLim = c(3, 55)
    }
    if(region == "World" | region == "world"){
      latLim = c(-90.00, 90.00)
    }
    if(length(grep("World|world|china|China", region)) != 1){
      stop("'latLim' is required if the 'region' to map is not China or World.")
    }
  }

  # a = 1; i = 1
    for(a in 1:length(genelist)){
      mapFileDir.temp = mapFileDir[grep(genelist[a], mapFileDir)]
      if(dir.exists(mapFileDir.temp)){
        mapFileDir.temp = paste0(mapFileDir.temp, "/", list.files(mapFileDir.temp))
      }
      if(length(grep(".frequency|.pdf|.log|.annotation.csv|.category.csv", mapFileDir.temp))>0){
        mapFileDir.temp = mapFileDir.temp[-grep(".frequency|.pdf|.log|.annotation.csv|.category.csv", mapFileDir.temp)]
      }
      if(!dir.exists(paste0(dirPath, "/GeoMappingResults/", genelist[a]))){
        dir.create(paste0(dirPath, "/GeoMappingResults/", genelist[a]))
      }
      mapFileDir.temp = sort(mapFileDir.temp)
      for(i in 1:length(mapFileDir.temp)){
        fileName = strsplit(mapFileDir.temp[i], "/")[[1]][length(strsplit(mapFileDir.temp[i], "/")[[1]])]
        fileName = gsub(".map|.csv|.txt|.xlsx", "", fileName,ignore.case = T)
        data = fread(mapFileDir.temp[i], header = T,na.strings = "")
        if(length(grep("Longitude|Latitude|longitude|latitude", names(data))) != 2){
          stop(paste0("Mapping file of ", mapFileDir[i] ," should contain at least 'Longitude' and 'Latitude'"))
        }
        data$order = NA
        for(m in 1:nrow(data)){
          data$order[m] = as.numeric(gsub("Hap", "", data$haplotype[m]))
        }
        data = data[order(data$order),]

        if(!is.null(dataFiltration)){
          if(length(grep("Class|class", names(data))) == 1){
            names(data) = gsub("class", "Class", names(data))
            data = data[which(data$Class == dataFiltration),]
          }else{
            cat(paste0("Mapping file of ", fileName, " DOES NOT contain 'Class', which is important for data filtration!"))
            next
          }
        }

        if(MappingDataType == "Type"){
          if(length(grep("Type", names(data))) == 1){
            data = data[,c("Longitude", "Latitude", "Type")]
            data_reshape = data.frame(dcast(data, Sample_ID+haplotype+Latitude+Longitude~Type, ),check.names = F)
          }else{
            cat(paste0("Mapping file: ", mapFileDir[i], " DOES NOT contain column of 'Type'!"))
            next
          }
        }

        if(MappingDataType == "haplotype"){
          data_reshape = data.frame(dcast(data, Latitude+Longitude~order, ),check.names = F)
        if(grepl("China",region,ignore.case = T)){
          category = fread(mapFileDir.temp[grep("class",mapFileDir.temp,ignore.case = T)[1]],header = T)
          category = data.frame(category[,c("Latitude","Longitude","Class")])
          category = distinct(category)
          category$Loc = paste0(category$Latitude,"_",category$Longitude)
          category$Latitude = NULL;category$Longitude = NULL
          data_reshape$Loc = paste0(data_reshape$Latitude,"_",data_reshape$Longitude)
          data_reshape = merge(data_reshape,category,by = "Loc",all.y= T)
          data_reshape = na.omit(data_reshape)
          data_reshape = data_reshape[grep("China",data_reshape$Class,ignore.case = T),]
          data_reshape$Class = NULL;data_reshape$Loc <- NULL
          data_reshape <- data_reshape[, c(c("Latitude","Longitude"), sort(as.numeric(setdiff(names(data_reshape),c("Latitude","Longitude")))))]
          names(data_reshape)[c(3:ncol(data_reshape))] = paste0("Hap",names(data_reshape)[c(3:ncol(data_reshape))])
          }
        }
        
          for(j in 1:nrow(data_reshape)){
            for(k in 1:2){
              if(length(grep("E", data_reshape[j,k])) == 1){
                data_reshape[j,k] = gsub("E", "", data_reshape[j,k])
              }
              if(length(grep("W", data_reshape[j,k])) == 1){
                data_reshape[j,k] = gsub("W", "", data_reshape[j,k])
                data_reshape[j,k] = paste0("-", data_reshape[j,k])
              }
              if(length(grep("N", data_reshape[j,k])) == 1){
                data_reshape[j,k] = gsub("N", "", data_reshape[j,k])
              }
              if(length(grep("S", data_reshape[j,k])) == 1){
                data_reshape[j,k] = gsub("S", "", data_reshape[j,k])
                data_reshape[j,k] = paste0("-", data_reshape[j,k])
              }
            }
          }
          for(j in 1:ncol(data_reshape)){
            data_reshape[,j] = as.numeric(data_reshape[,j])
          }
          categoryNum = length(names(data_reshape)[3:ncol(data_reshape)])

          for(i in 1:nrow(data_reshape)){
            data_reshape$count[i] =sum(data_reshape[i, c(3:(ncol(data_reshape)-1))])
          }
          data_reshape$radius = data_reshape$count*(pieMax/max(data_reshape$count))

          ## plotting
          pdf(paste0(dirPath, "/GeoMappingResults/", genelist[a], "/", fileName, ".geographically.mapped.to.", region, ".", type, ".", pieMax, ifelse(MappingDataType == "haplotype", ".haplotype.pdf", ".Type.pdf")), width = pdfWidth, height = pdfHeight)
          p = ggplot()+
            geom_sf(data = map,fill = landColor, size = mapBorderSize) +
            geom_scatterpie(data = data_reshape,
                            aes(x = Longitude,
                                y = Latitude,
                                r = radius,colour = color_array),
                            cols = names(data_reshape)[3:(ncol(data_reshape)-2)],
                            linetype = pieLinetype,
                            size = pieLinesize,
                            color = pieBorderColor,
                            alpha = pieAlpha)  +
            coord_sf(xlim = longLim, ylim = latLim, expand = expand, label_graticule = label_graticule) +
            geom_scatterpie_legend(r = data_reshape$radius,
                                   x = longLim[1] + 10,
                                   y = latLim[1] + 10,
                                   n=5,
                                   labeller = function(x) {as.integer(x*max(data_reshape$count)/pieMax)}) +
            labs(x = "Longitude", y = "Latitude") +
            theme_bw()+
            theme(panel.grid = element_blank(),
                  panel.background = element_rect(fill = oceanColor))
          p = p + theme(plot.margin=unit(rep(plotMargin,4),'lines'))
          print(p)
          dev.off()
      }
  }
}

