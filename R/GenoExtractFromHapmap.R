#' This is some description of this function.
#' @title to extract samples from raw.hapmap data
#'
#' @description By using this package, you could use functions of GenoExtractFromHapmap to extract samples from raw.hapmap data
#'
#' @details see above
#'
#' @param trait: the trait you are working with;
#' @param categoryFile: file stores category of each sample;
#' @param genoFile: a folder stores the hapmap file of each chromosome, or a single hapmap file of specific chromosome.

#' @return files extracted and stored in the folder created by default.
#' @export GenoExtractFromHapmap
#' @examples GenoExtractFromHapmap(trait = "Trait",categoryFile = "./Trait.category.csv",genoFile = "./soy_hapmap/")

GenoExtractFromHapmap <- function(trait = NULL, categoryFile = NULL, genoFile = NULL, ...){
  library(data.table)
  dir.path <- getwd()
  if(is.null(trait)) trait = "unknown"
  if(is.null(categoryFile)) stop("'categoryFile' is required!")
  if(is.null(genoFile)) stop("'genoFile' is required!")
  if(!dir.exists(paste0(dir.path, "/", trait, "_geno.extracted.from.hapmap"))){
    dir.create(paste0(dir.path, "/", trait, "_geno.extracted.from.hapmap"))
  }

  sample <- as.data.frame(fread(categoryFile, header =T))
  sample <- data.frame(sample[,1])
  sample[,1] <- sample[!duplicated(sample[,1]),1]
  colnames(sample) <- c("sampleID")

  if(dir.exists(genoFile)){
    hapmapList = paste0(genoFile, list.files(genoFile))
  }else{
    hapmapList = genoFile
  }

  vector2string = function(vec){
    if(length(vec) > 1){
      str = vec[1]
      for(i in 2:length(vec)){
        str = paste0(str, ",", vec[i])
      }
    }else{
      str = vec
    }
    return(str)
  }

  # i = 1;j=1
  for(i in 1:length(hapmapList)){
    fileName <- gsub(".hapmap|.txt|.csv|.hmp", "", strsplit(hapmapList[i], "/")[[1]][length(strsplit(hapmapList[i], "/")[[1]])],ignore.case = T)
    hapmap <- fread(hapmapList[i], header = T)
    col_name <- colnames(hapmap)
    setkeyv(hapmap, names(hapmap))
    eSample = intersect(sample$sampleID, names(hapmap))
    unSample = setdiff(sample$sampleID, names(hapmap))
    unSample = vector2string(vec = unSample)
    exSample = vector2string(vec = eSample)
    names(hapmap)[1:11] <- c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panel","QCcode")
    selectCol = c(names(hapmap)[1:11],eSample)
    hapmap = data.frame(hapmap,check.names = F)
    hapmapSelected = hapmap[,selectCol]
    hapmapSelected$chrom = gsub("Gm0|Chr0|Gm|Chr","",hapmapSelected$chrom,ignore.case = T)
    Chr = unique(hapmapSelected$chrom)
    for(j in 1:length(Chr)){
      hapmapSelectedChr = hapmapSelected[which(hapmapSelected$chrom==Chr[j]),]
      fwrite(hapmapSelectedChr, paste0(dir.path, "/", trait, "_geno.extracted.from.hapmap/", trait, "_hapmap_Gm",Chr[j], ".hapmap"), row.names = F, col.names = T, quote = F, sep = "\t")  
      sink(paste0(dir.path, "/", trait, "_geno.extracted.from.hapmap/", trait, "_hapmap_Gm",Chr[j], ".geno.extraction.log"), type = "output", append = T)
      writeLines(paste0("\n###### Genotype extraction from ", fileName, " ######"))
      writeLines(paste0("Trait: ", trait))
      writeLines(paste0("SampleList: ", categoryFile))
      writeLines(paste0("genoFile: ", hapmapList[i]))
      writeLines(paste0("Date: ", date(), "\n"))
      writeLines(paste0("Not existed sample(s): ", unSample, "\n"))
      writeLines(paste0("Identified samples: ", exSample, "\n"))
      writeLines("### THE END ###\n")
      sink(type = "output")
       }
  }
}
