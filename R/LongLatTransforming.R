#' This is some description of this function.
#' @title to transform longitude or latitude into different format, either from dms to dd, or from dd to dms
#'
#' @description By using this package, you could use functions of LongLatTransforming to transform longitude or latitude into different format, either from dms to dd, or from dd to dms
#'
#' @details In the LongLatTransforming, the input data could either be a regular format, or files containing longitude and latitude.
#'
#' @param longLat: input data to be transformed. It could be either value with or without 'E|N|W|S' to indicate direction, or file including longitude or latitude column. For example, 
#'
#' @return value or file, depends on the input data.
#' @export LongLatTransforming
#' 
#' @examples LongLatTransforming("-116°25'W")
#' 
# LongLatTransforming("-116°15'17.9999999999836''E")
# LongLatTransforming("-116.255S")

LongLatTransforming = function(longLat = NULL, ...){
  library(data.table)
  library(stringr)
  library(xlsx)
  library(magrittr)
  library(sp)
  dirPath = getwd()
  
  if(!dir.exists(paste0(dirPath, "/LongLatTransformation/"))){
    dir.create(paste0(dirPath, "/LongLatTransformation/"))
  }
  
  dms2dd <- function(dms){  
    orientation = str_extract_all(dms, "E|W|N|S")[[1]]
    orientation1 = "W"
    dms1 = gsub("E|W|N|S", "", dms)
    char.dms <- str_c(dms1,orientation1) %>%    
      str_replace(replacement = "d", pattern = "°") %>% # 这里用英文符号'和"替换了中文符号′与″
      str_replace(replacement = "'", pattern = "′") %>% # 注意符号是中文字符还是英文字符
      str_replace(replacement = "\"", pattern = "″") %>%  # "与″ '与′ 是不一样的。 DMS类中需要英文字符
      str_replace(replacement = "\"", pattern = "〃")
    dms2 <- char2dms(char.dms)
    d <- dms2@deg # 从dms类中获取度、分、秒的数值
    m <- dms2@min
    s <- dms2@sec
    if(length(orientation) != 0){
      if(d < 0){
        dd = abs(d) + m/60 + s/3600
        dd = paste0("-", dd)
      }
      if(d >= 0){
        dd = d + m/60 + s/3600
      }
      dd = paste0(dd, orientation)
    }
    if(length(orientation) == 0){
      if(d < 0){
        dd = abs(d) + m/60 + s/3600
        dd = paste0("-", dd)
      }
      if(d >= 0){
        dd = d + m/60 + s/3600
      }
    }
    return(dd)
  }
  
  ddTodms <- function(dd){
    orientation = str_extract_all(dd, "E|W|N|S")[[1]]
    dd = gsub("E|W|N|S", "", dd)
    dms = dd2dms(as.numeric(dd))
    d = dms@deg
    m = dms@min
    s = dms@sec
    if(as.numeric(dd) < 0){
      d = paste0("-", d)
    }
    if(s != 0){
      if(length(orientation) > 0){
        dms = paste0(d, "°", m, "'", s, "''", orientation)
      }else{
        dms = paste0(d, "°", m, "'", s, "''")
      }
    }else{
      if(length(orientation) > 0){
        dms = paste0(d, "°", m, "'", orientation)
      }else{
        dms = paste0(d, "°", m, "'")
      }
    }
    return(dms)
  }
  
  if(!file.exists(longLat) & !dir.exists(longLat)){
    if(str_extract_all(longLat, "°") == "°"){
      dd = dms2dd(longLat)
      print(dd)
    }
    if(str_extract_all(longLat, "°") != "°"){
      dms = ddTodms(longLat)
      print(dms)
    }
  }
  
  if(length(grep(".xlsx", longLat)) == 1){
    longlat = read.xlsx(longLat, sheetIndex = 1)
    longlatcol = grep("longitude|Longitude|latitude|Latitude|Long|long|Lat|lat", names(longlat))
    fileName = gsub(".xlsx", "", strsplit(longLat, "/")[[1]][length(strsplit(longLat, "/")[[1]])])
    longlat$Longitude_transformed = NA
    longlat$Latitude_transformed = NA
    for(i in 1:length(longlatcol)){
      for(j in 1:nrow(longlat)){
        if(length(str_extract_all(names(longlat)[longlatcol[i]], "longitude|Longitude|Long|long")[[1]]) == 1){
          if(length(str_extract_all(longlat[j,longlatcol[i]], "°")[[1]]) == 1){
            longlat$Longitude_transformed[j] = dms2dd(longlat[j,longlatcol[i]])
          }
          if(length(str_extract_all(longlat[j,longlatcol[i]], "°")[[1]]) != 1){
            longlat$Longitude_transformed[j] = ddTodms(longlat[j,longlatcol[i]])
          }
        }
        if(length(str_extract_all(names(longlat)[longlatcol[i]], "latitude|Latitude|Lat|lat")[[1]]) == 1){
          if(length(str_extract_all(longlat[j,longlatcol[i]], "°")[[1]]) == 1){
            longlat$Latitude_transformed[j] = dms2dd(longlat[j,longlatcol[i]])
          }
          if(length(str_extract_all(longlat[j,longlatcol[i]], "°")[[1]]) != 1){
            longlat$Latitude_transformed[j] = ddTodms(longlat[j,longlatcol[i]])
          }
        }
      }
    }
  }
  if(length(grep(".csv|.txt", longLat)) == 1){
    longlat = fread(longLat, header = T)
    longlatcol = grep("longitude|Longitude|latitude|Latitude|Long|long|Lat|lat", names(longlat))
    fileName = gsub(".csv|.txt", "", strsplit(longLat, "/")[[1]][length(strsplit(longLat, "/")[[1]])])
    longlat$Longitude_transformed = NA
    longlat$Latitude_transformed = NA
    for(i in 1:length(longlatcol)){
      for(j in 1:nrow(longlat)){
        if(length(str_extract_all(names(longlat)[longlatcol[i]], "longitude|Longitude|Long|long")[[1]]) == 1){
          if(length(str_extract_all(longlat[j,longlatcol[i]], "°")[[1]]) == 1){
            longlat$Longitude_transformed[j] = dms2dd(longlat[j,longlatcol[i]])
          }
          if(length(str_extract_all(longlat[j,longlatcol[i]], "°")[[1]]) != 1){
            longlat$Longitude_transformed[j] = ddTodms(longlat[j,longlatcol[i]])
          }
        }
        if(length(str_extract_all(names(longlat)[longlatcol[i]], "latitude|Latitude|Lat|lat")[[1]]) == 1){
          if(length(str_extract_all(longlat[j,longlatcol[i]], "°")[[1]]) == 1){
            longlat$Latitude_transformed[j] = dms2dd(longlat[j,longlatcol[i]])
          }
          if(length(str_extract_all(longlat[j,longlatcol[i]], "°")[[1]]) != 1){
            longlat$Latitude_transformed[j] = ddTodms(longlat[j,longlatcol[i]])
          }
        }
      }
    }
  }

  fwrite(longlat, file = paste0(dirPath, "/LongLatTransformation/", fileName, ".LongLat.transformed.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
}