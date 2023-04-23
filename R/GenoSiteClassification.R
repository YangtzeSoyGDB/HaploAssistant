#' This is some description of this function.
#' @title to classify SNP sites according to relative position to a certain gene
#'
#' @description By using this package, you could use functions of GenoSiteClassification to classify SNP sites according to relative position to a certain gene
#'
#' @details see above
#'
#' @param annoFile: a SNP annotation folder from the fuction SNPAnnotation;
#' @param hapmapFile: regional hapmap file extracted from hapmap, could be either a file or a folder storing these type of files;
#'
#' @return files extracted and stored in the folder created by default.
#' @export GenoSiteClassification
#' @examples GenoSiteClassification(annoFile = "./SNP.Annotation/", hapmapFile = "./Hapmap.Regional.Extracted.Files/")

GenoSiteClassification <- function(annoFile = NULL, hapmapFile = NULL, ...){
  library(data.table)
  library(stringr)
  dirPath = getwd()
  if(is.null(annoFile)) stop("'annoFile' is required!")
  if(is.null(hapmapFile)) stop("'hapmapFile' is required!")

  if(!dir.exists(paste0(dirPath, "/", "Geno.Site.Classification"))){
    dir.create(paste0(dirPath, "/", "Geno.Site.Classification"))
  }

  if(dir.exists(hapmapFile)){
    genoList = paste0(hapmapFile, list.files(hapmapFile))
  }else{
    genoList = hapmapFile
  }
  if(length(grep(".log", genoList))){
    genoList = genoList[-grep(".log", genoList)]
  }

  if(dir.exists(annoFile)){
    annoList = paste0(annoFile, list.files(annoFile))
  }else{
    annoList = annoFile
  }
  type = c("upstream|downstream|intergenic|gene.related|intron.utr|in.total")
  annoList = annoList[grep(".csv", annoList)]
  if(length((grep(".location.csv", annoList)))>0){
    annoList = annoList[-grep(".location.csv", annoList)]
  }
  # i = 1; j = 3;k = 1
  for(i in 1:length(genoList)){
    geno = fread(genoList[i], header = T)
    geneName = gsub(".hmp|.hapmap", "", strsplit(genoList[i], "/")[[1]][length(strsplit(genoList[i], "/")[[1]])],ignore.case = T)
    geneName <- strsplit(unlist(geneName),"_", fixed = TRUE)[[1]][[1]]
    annoList.temp = annoList[grep(geneName, annoList)]

    for(j in 1:length(annoList.temp)){
      anno = fread(annoList.temp[j], header = T)
      if(nrow(anno)==0) next
      annoType = str_extract_all(annoList.temp[j], type)[[1]]
      anno$Pos = as.numeric(anno$Pos)
      names(anno) = gsub("Pos", "pos", names(anno),ignore.case = F)
      GeneID <- unique(anno$Gene_ID)
      if(length(GeneID)==1){
        #genos = geno[,][pos %in% anno$pos]
        genos = geno[which(geno$pos %in% anno$pos),]
        annogeno = merge(x = anno, y = genos, by = c("pos"), all.y = T)
        fwrite(annogeno, file = paste0(dirPath, "/", "Geno.Site.Classification", "/", GeneID, "_", annoType, ".csv"), row.names = F, col.names = T, quote = F, sep = "\t")
      }else{
        for(k in 1:length(GeneID)){
          annoTemp = anno[which(anno$Gene_ID==GeneID[k]),]
          genos = geno[,][pos %in% annoTemp$pos]
          annogeno = merge(x = anno, y = annoTemp, by = c("pos"), all.y = T)
          fwrite(annogeno, file = paste0(dirPath, "/", "Geno.Site.Classification", "/", GeneID[k], "_", annoType, ".csv"), row.names = F, col.names = T, quote = F, sep = "\t")
        }
      }



          }
  }
}



