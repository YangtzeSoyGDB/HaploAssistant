#' This is some description of this function.
#' @title generate .nex file stores DATA and TRAITS for haplotype network inference using PopART
#'
#' @description By using this function, you could generate .nex file stores DATA and TRAITS for haplotype network inference using PopART
#'
#' @details see above
#'
#' @param haploResults: haplotype category and frequency of haplotype analysis results, could be either a folder stores gene list files or a single file of gene list;
#' @param categoryFile: longitude and latitude of each location.

#' @return files extracted and stored in the folder created by default.
#' @export NexForPopART
#' @examples NexForPopART(haploResults = "./Haplotype.Analysis/", categoryFile = "./Trait.category.csv")

NexForPopART <- function(haploResults = NULL, categoryFile = "none", ...){
  library(data.table)
  library(stringr)
  dirPath <- getwd()
  if(is.null(haploResults)) stop("'haploResults' is required!")
  if(!dir.exists(paste0(dirPath, "/", "PopARTInput"))){
    dir.create(paste0(dirPath, "/", "PopARTInput"))
  }

  if(dir.exists(haploResults)){
    haploList = paste0(haploResults, list.files(haploResults))
  }else{
    haploList = haploResults
  }
  genelist = data.frame(str_split(as.list(haploList), "/"))[]
  genelist = as.character(genelist[nrow(genelist),])
  # a= 1
  if(categoryFile == "none"){
    for(a in 1:length(genelist)){
      if(!dir.exists(paste0(dirPath, "/", "PopARTInput/", genelist[a]))){
        dir.create(paste0(dirPath, "/", "PopARTInput/", genelist[a]))
      }
      haploList.temp = haploList[grep(genelist[a], haploList)]
      if(dir.exists(haploList.temp)){
        haploList.temp = paste0(haploList.temp, "/", list.files(haploList.temp))
      }
      haploList.temp = haploList.temp[grep("category", haploList.temp)]
      for(b in 1:length(haploList.temp)){
        fileName <- gsub(".haplotype.category.csv", "", strsplit(haploList.temp[b], "/")[[1]][length(strsplit(haploList.temp[b], "/")[[1]])])
        haplo <- fread(haploList.temp[b], header = T)
        haplo_col = names(haplo)
        haplo = data.frame(haplo)
        names(haplo) = haplo_col
        ref.seq <- strsplit(haplo$Pattern[1], "")[[1]]
        for(i in 1:nrow(haplo)){
          temp1 <- strsplit(haplo$Pattern[i], "")[[1]]
          temp1[grep("-", temp1)] <- ref.seq[grep("-", temp1)]
          temp2 <- temp1[1]
          for(j in 2:length(temp1)){
            temp2 <- paste0(temp2, temp1[j])
          }
          haplo$Pattern[i] <- temp2
        }
        traitlabels <- names(haplo)[4:ncol(haplo)]
        matrix.body <- data.frame(matrix(nrow = nrow(haplo), ncol = 2))
        colnames(matrix.body) <- c("Haplotype", "Frequency")
        matrix.body$Haplotype <- paste0("Hap", seq(1:nrow(haplo)))
        for(i in 1:nrow(haplo)){
          temp1 <- haplo[i,4]
          for(j in 5:ncol(haplo)){
            temp1 <- paste(temp1, haplo[i,j], sep = ",")
          }
          matrix.body$Frequency[i] <- temp1
        }
        taxa.1 <- c("BEGIN TAXA;")
        taxa.2 <- paste0("DIMENSIONS NTAX=", nrow(haplo), ";")
        taxa.3 <- c("TAXLABELS")
        taxa.head <- paste0(taxa.1, "\n", taxa.2, "\n", "\n", taxa.3)
        nex.1 <- paste0("#NEXUS", "\n")
        nex.1.1 <- paste0("\n", "BEGIN TRAITS;")
        nex.2 <- paste0("  ", "Dimensions NTRAITS=", length(traitlabels), ";")
        nex.3 <- paste0("  ", "Format labels=yes missing=? separator=Comma;")
        nex.4 <- c("TraitLabels")
        for(j in 1:length(traitlabels)){
          nex.4 <- paste(nex.4, traitlabels[j], sep = " ")
        }
        nex.4 <- paste0("  ", nex.4, ";")
        nex.5 <- paste0("  ", "Matrix")
        blank.row <- c("")
        nex.6 <- c(";")
        nex.7 <- c("END;")
        nex.head <- paste0(nex.1.1, "\n", nex.2, "\n", nex.3, "\n", nex.4)
        nex.end <- paste0(nex.6, "\n","\n", nex.7)
        matrix.body.lines <- paste(matrix.body$Haplotype, matrix.body$Frequency, sep = " ")

        hap.1 <- c("BEGIN DATA;")
        hap.2 <- paste0("  ", "DIMENSIONS NTAX=", nrow(haplo), " nchar=", nchar(haplo$Pattern[1]), ";")
        hap.3 <- paste0("  ", "Format datatype=DNA missing=N gap=-;")
        hap.head <- paste0(hap.1, "\n", hap.2, "\n", hap.3, "\n", nex.5)
        hap.body <-  paste(matrix.body$Haplotype, haplo$Pattern, sep = " ")

        sink(paste0(dirPath, "/PopARTInput/", genelist[a], "/", fileName, ".nex"), type = "output", append = TRUE)
        ### sequence
        writeLines(nex.1)
        writeLines(hap.head)
        writeLines(hap.body)
        writeLines(nex.end)
        writeLines(blank.row)
        ### trait
        writeLines(nex.head)
        #writeLines(blank.row)
        writeLines(nex.5)
        writeLines(matrix.body.lines)
        writeLines(nex.end)
        sink(type = "output")
      }
    }
  }
  if(categoryFile != "none"){
    longLat <- data.frame(fread(categoryFile, header = T))
    longLat = longLat[,grep("Longitude|longitude|Latitude|latitude|Loc|loc|location|Location|Origination|origination", names(longLat))]
    longLat = longLat[!duplicated(longLat),]
    names(longLat) = gsub("longitude", "Longitude", names(longLat))
    names(longLat) = gsub("latitude", "Latitude", names(longLat))
    names(longLat) = gsub("Loc|loc|location|Location|Origination|origination", "Loc", names(longLat))
    # i =1
    for(i in 1:nrow(longLat)){
      if(length(grep("W",longLat$Longitude[i])) != 0){
        longLat$Longitude[i] <- gsub("W", "", longLat$Longitude[i],ignore.case = T)
        longLat$Longitude[i] <- paste0("-", longLat$Longitude[i])
        longLat$Longitude[i] <- gsub("°", "", longLat$Longitude[i])
      }
      if(length(grep("E",longLat$Longitude[i])) != 0){
        longLat$Longitude[i] <- gsub("E", "", longLat$Longitude[i],ignore.case = T)
        longLat$Longitude[i] <- gsub("°", "", longLat$Longitude[i])
      }
      if(length(grep("S", longLat$Latitude[i])) != 0){
        longLat$Latitude[i] <- gsub("S", "", longLat$Latitude[i],ignore.case = T)
        longLat$Latitude[i] <- paste0("-", longLat$Latitude[i])
        longLat$Latitude[i] <- gsub("°", "", longLat$Latitude[i])
      }
      if(length(grep("N", longLat$Latitude[i])) != 0){
        longLat$Latitude[i] <- gsub("N", "", longLat$Latitude[i],ignore.case = T)
        longLat$Latitude[i] <- gsub("°", "", longLat$Latitude[i])
      }
    }

    # a = 1; b = 1
    for(a in 1:length(genelist)){
      if(!dir.exists(paste0(dirPath, "/PopARTInput/", genelist[a]))){
        dir.create(paste0(dirPath, "/PopARTInput/", genelist[a]))
      }
      haploList.temp = haploList[grep(genelist[a], haploList)]
      if(dir.exists(haploList.temp)){
        haploList.temp = paste0(haploList.temp, "/", list.files(haploList.temp))
      }
      haploList.temp = haploList.temp[grep("category", haploList.temp,ignore.case = T)]
      for(b in 1:length(haploList.temp)){
        fileName <- gsub(".haplotype.category.csv", "", strsplit(haploList.temp[b], "/")[[1]][length(strsplit(haploList.temp[b], "/")[[1]])])
        haplo <- fread(haploList.temp[b], header = T)
        haplo_col = names(haplo)
        haplo = data.frame(haplo)
        names(haplo) = haplo_col
        ref.seq <- strsplit(haplo$Pattern[1], "")[[1]]
        for(i in 1:nrow(haplo)){
          temp1 <- strsplit(haplo$Pattern[i], "")[[1]]
          temp1[grep("-", temp1)] <- ref.seq[grep("-", temp1)]
          temp2 <- temp1[1]
          for(j in 2:length(temp1)){
            temp2 <- paste0(temp2, temp1[j])
          }
          haplo$Pattern[i] <- temp2
        }
        traitlabels <- data.frame(names(haplo)[4:ncol(haplo)])
        colnames(traitlabels) <- c("Loc")
        trait.labels <- merge(x = traitlabels, y = longLat, all = TRUE)
        #### if category was Loc, Longitude and Latitude would be attached to each loc
        if(nrow(traitlabels) == nrow(trait.labels)){
          trait.row <- c()
          for(i in 1:nrow(traitlabels)){
            temp1 <- grep(traitlabels[i,1], trait.labels$Loc)
            trait.row <- c(trait.row, temp1)
          }
          trait.labels <- trait.labels[trait.row,]
          matrix.body <- data.frame(matrix(nrow = nrow(haplo), ncol = 2))
          colnames(matrix.body) <- c("Haplotype", "Frequency")
          matrix.body$Haplotype <- paste0("Hap", seq(1:nrow(haplo)))
          for(i in 1:nrow(haplo)){
            temp1 <- haplo[i, 4]
            for(j in 5:ncol(haplo)){
              temp1 <- paste(temp1, haplo[i,j], sep = ",")
            }
            matrix.body$Frequency[i] <- temp1
          }
          taxa.1 <- c("BEGIN TAXA;")
          taxa.2 <- paste0("DIMENSIONS NTAX=", nrow(haplo), ";")
          taxa.3 <- c("TAXLABELS")
          taxa.head <- paste0(taxa.1, "\n", taxa.2, "\n", "\n", taxa.3)
          nex.1 <- paste0("#NEXUS", "\n")
          nex.1.1 <- paste0("\n", "BEGIN TRAITS;")
          nex.2 <- paste0("  ", "Dimensions NTRAITS=", nrow(traitlabels), ";")
          nex.3 <- paste0("  ", "Format labels=yes missing=? separator=Comma;")
          nex.4 <- c("TraitLabels")
          nex.4lat <- c("TraitLatitude")
          nex.4long <- c("TraitLongitude")
          for(j in 1:nrow(trait.labels)){
            nex.4 <- paste(nex.4, traitlabels$Loc[j], sep = " ")
            nex.4lat <- paste(nex.4lat, trait.labels$Latitude[j], sep = " ")
            nex.4long <- paste(nex.4long, trait.labels$Longitude[j], sep = " ")
          }
          nex.4 <- paste0("  ", nex.4, ";")
          nex.4lat <- paste0("  ", nex.4lat, ";")
          nex.4long <- paste0("  ", nex.4long, ";")
          nex.5 <- paste0("  ", "Matrix")
          blank.row <- c("")
          nex.6 <- c(";")
          nex.7 <- c("END;")
          nex.head <- paste0(nex.1.1, "\n", nex.2, "\n", nex.3, "\n", nex.4lat, "\n", nex.4long, "\n", nex.4)
          nex.end <- paste0(nex.6, "\n","\n", nex.7)
          matrix.body.lines <- paste(matrix.body$Haplotype, matrix.body$Frequency, sep = " ")

          hap.1 <- c("BEGIN DATA;")
          hap.2 <- paste0("  ", "DIMENSIONS NTAX=", nrow(haplo), " nchar=", nchar(haplo$Pattern[1]), ";")
          hap.3 <- paste0("  ", "Format datatype=DNA missing=N gap=-;")
          hap.head <- paste0(hap.1, "\n", hap.2, "\n", hap.3, "\n", nex.5)
          hap.body <-  paste(matrix.body$Haplotype, haplo$Pattern, sep = " ")

          sink(paste0(dirPath,  "/PopARTInput/", genelist[a], "/", fileName, ".nex"), type = "output", append = TRUE)
          ### sequence
          writeLines(nex.1)
          writeLines(hap.head)
          writeLines(hap.body)
          writeLines(nex.end)
          writeLines(blank.row)
          ### trait
          writeLines(nex.head)
          #writeLines(blank.row)
          writeLines(nex.5)
          writeLines(matrix.body.lines)
          writeLines(nex.end)
          sink(type = "output")
        }else{
          matrix.body <- data.frame(matrix(nrow = nrow(haplo), ncol = 2))
          colnames(matrix.body) <- c("Haplotype", "Frequency")
          matrix.body$Haplotype <- paste0("Hap", seq(1:nrow(haplo)))
          for(i in 1:nrow(haplo)){
            temp1 <- haplo[i, 4]
            for(j in 5:ncol(haplo)){
              temp1 <- paste(temp1, haplo[i,j], sep = ",")
            }
            matrix.body$Frequency[i] <- temp1
          }
          taxa.1 <- c("BEGIN TAXA;")
          taxa.2 <- paste0("DIMENSIONS NTAX=", nrow(haplo), ";")
          taxa.3 <- c("TAXLABELS")
          taxa.head <- paste0(taxa.1, "\n", taxa.2, "\n", "\n", taxa.3)
          nex.1 <- paste0("#NEXUS", "\n")
          nex.1.1 <- paste0("\n", "BEGIN TRAITS;")
          nex.2 <- paste0("  ", "Dimensions NTRAITS=", nrow(traitlabels), ";")
          nex.3 <- paste0("  ", "Format labels=yes missing=? separator=Comma;")
          nex.4 <- c("TraitLabels")
          for(j in 1:nrow(traitlabels)){
            nex.4 <- paste(nex.4, traitlabels[j,1], sep = " ")
          }
          nex.4 <- paste0("  ", nex.4, ";")
          nex.5 <- paste0("  ", "Matrix")
          blank.row <- c("")
          nex.6 <- c(";")
          nex.7 <- c("END;")
          nex.head <- paste0(nex.1.1, "\n", nex.2, "\n", nex.3, "\n", nex.4)
          nex.end <- paste0(nex.6, "\n","\n", nex.7)
          matrix.body.lines <- paste(matrix.body$Haplotype, matrix.body$Frequency, sep = " ")

          hap.1 <- c("BEGIN DATA;")
          hap.2 <- paste0("  ", "DIMENSIONS NTAX=", nrow(haplo), " nchar=", nchar(haplo$Pattern[1]), ";")
          hap.3 <- paste0("  ", "Format datatype=DNA missing=N gap=-;")
          hap.head <- paste0(hap.1, "\n", hap.2, "\n", hap.3, "\n", nex.5)
          hap.body <-  paste(matrix.body$Haplotype, haplo$Pattern, sep = " ")

          sink(paste0(dirPath,  "/PopARTInput/", genelist[a], "/", fileName, ".nex"), type = "output", append = TRUE)
          ### sequence
          writeLines(nex.1)
          writeLines(hap.head)
          writeLines(hap.body)
          writeLines(nex.end)
          writeLines(blank.row)
          ### trait
          writeLines(nex.head)
          #writeLines(blank.row)
          writeLines(nex.5)
          writeLines(matrix.body.lines)
          writeLines(nex.end)
          sink(type = "output")
        }
      }
    }
  }
}

