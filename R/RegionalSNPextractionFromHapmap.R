#' This is some description of this function.
#' @title extract SNP genotype of of a certain regionFiles either to specific trait or not from hapmap files
#'
#' @description By using this function, you could extract SNP genotype of of a certain regionFiles either to specific trait or not from hapmap files
#'
#' @details see above
#'
#' @param sampleHapmapFile: hapmap file after extracting samples, could be either a single hapmap file or a folder stores hapmap files;
#' @param gffFile: generic feature format file of certain genome, and this file defines the start and end position of specific gene, which should contain at least 9 columns such as "seqid", "source", "type", "start", "end", "score", "strand", "phase" and "attributes".
#' @param geneModel: gene model you want to perform haplotype analysis, or a file that stores gene names with the column names "Gene_ID".

#' @return files extracted and stored in the folder created by default
#' @export RegionalSNPextractionFromHapmap
#' @examples RegionalSNPextractionFromHapmap(sampleHapmapFile = "./Trait_geno.extracted.from.hapmap/",
#'                                           gffFile = "./Wm82.a2.v1.gene.gff",
#'                                           geneModel = "Glyma.18G056600")

RegionalSNPextractionFromHapmap <- function(sampleHapmapFile = NULL, gffFile= NULL, geneModel = NULL, ...){
  library(stringr)
  library(data.table)
  
  if(is.null(sampleHapmapFile)) stop("'sampleHapmapFile' is required!")
  if(is.null(gffFile) & is.null(geneModel)) stop("'gffFile' and 'geneModel' is required!")
  dirPath <- getwd()
  
  if(!dir.exists(paste0(dirPath, "/", "Hapmap.Regional.Extracted.Files"))){
    dir.create(paste0(dirPath, "/", "Hapmap.Regional.Extracted.Files"))
  }
  
  if(dir.exists(sampleHapmapFile)){
    hapmapList = paste0(sampleHapmapFile, list.files(sampleHapmapFile))
    if(length(grep(".log", hapmapList))>0){
      hapmapList = sort(hapmapList[-grep(".log", hapmapList)])
    }
  }else{
    hapmapList = sampleHapmapFile
  }
  
  if(file.exists(geneModel)){
    geneList = fread(geneModel, header = T)
    names(geneList) = gsub("gene_id|geneID", "Gene_ID", names(geneList))
    geneList <- geneList[!duplicated(geneList),]
    geneList <- na.omit(geneList)
    if(length(grep("gene_id|geneID|Gene_ID", names(geneList))) == 0){
      stop("'geneModel' file should contain 'Gene_ID' (column)!")
    }
  }else{
    geneList = data.frame(geneModel)
    names(geneList) = "Gene_ID"
  }
  
  
  gff = data.frame(fread(gffFile, header = T))
  names(gff) = c("chr","source","type","start","end","score","strand","phase","attributes")
  gff.gene = merge(gff,geneList,by.x="attributes",by.y="Gene_ID",all.y=T)
  gff.gene = na.omit(gff.gene)
  if(nrow(gff.gene)==nrow(geneList)){
    # i =1
    gff.gene = data.frame(matrix(ncol=9,nrow=0))
    for(i in 1:nrow(geneList)){
      gff.temp = gff[grep(geneList$Gene_ID[i],gff$attributes),]
      gff.gene = rbind(gff.gene,gff.temp)
    }
    if(length(grep("gene",gff.gene$type,ignore.case = T))==nrow(geneList)){
      gff.gene <- gff.gene[grep("gene",gff.gene$type,ignore.case = T),]
    }else{
      gff.gene <- gff.gene[grep("mRNA",gff.gene$type,ignore.case = T),]
    }
  }
  block = data.frame(matrix(ncol = 5,nrow=nrow(gff.gene)))
  names(block) = c("Gene_ID","Chr","Start","End","Formated")
  block$Gene_ID = gff.gene$attributes;block$Chr = gff.gene$chr
  block$Start = gff.gene$start;block$End = gff.gene$end;
  block$Formated = paste0("Chr",block$Chr,":",block$Start,"..",block$End)
  
  sink(paste0(dirPath, "/","Hapmap.Regional.Extracted.Files/", "regional.geno.extraction.from.hapmap.log"), type = "output", append = T)
  writeLines(paste0("\n###### regional geno extraction from hapmap ######\n"))
  writeLines(paste0("sampleHapmapFile: ", sampleHapmapFile))
  if(!is.null(gffFile)) writeLines(paste0("gffFile: ", gffFile))
  if(!is.null(geneModel)) writeLines(paste0("geneModel: ", gffFile))
  writeLines(paste0("Date: ", date(), "\n"))
  sink(type = "output")
  
  block$Chr = gsub("Chr0|Gm0|Chr|Gm", "", block$Chr,ignore.case = T)
  chrList = unique(block$Chr)
  
  # i = 1;j = 1
  for(i in 1:length(chrList)){
    hapmap = fread(hapmapList[grep(paste0(chrList[i],"[.]"), hapmapList)], header = T)
    block.temp = block[which(block$Chr == chrList[i]),]
    for(j in 1:nrow(block.temp)){
      hapmap =  hapmap[which(hapmap$chrom == block.temp$Chr[j]),]
      hapmaps = hapmap[which((hapmap$pos >= as.numeric(block.temp$Start[j])) & (hapmap$pos <= as.numeric(block.temp$End[j]))),]
      #hapmaps = hapmap[,][chrom == block.temp$Chr[j]][pos >= as.numeric(block.temp$Start[j])][pos <= as.numeric(block.temp$End[j])]
      if(nrow(hapmaps) == 0){
        sink(paste0(dirPath, "/","Hapmap.Regional.Extracted.Files/", "regional.geno.extraction.from.hapmap.log"), type = "output", append = T)
        writeLines(paste0(block.temp$Gene_ID[j], "_Chr", block.temp$Chr[j], "_", block.temp$Start[j], "..", block.temp$End[j], ": No SNPs identified."))
        sink(type = "output")
      }
      if(nrow(hapmaps) >= 1){
        hapmaps$pos = as.numeric(hapmaps$pos)
        hapmaps<- hapmaps[order(hapmaps$pos),]
        fwrite(hapmaps, file = paste0(dirPath, "/","Hapmap.Regional.Extracted.Files/", block.temp$Gene_ID[j], "_Chr", block.temp$Chr[j], "_", block.temp$Start[j], "..", block.temp$End[j], ".hapmap"), row.names = F, col.names = T, quote = F, sep = "\t")
      }
    }
  }
}

