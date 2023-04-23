#' This is some description of this function.
#' @title to identify gene related SNP and annotation
#'
#' @description By using this package, you could use the function of GeneSpecificSNPannotation to identify gene related SNP and annotation.
#'
#' @details see above
#'
#' @param SNPAnnoFile: SNP annotation files for users self sequencing;
#' @param geneModel: gene model you want to perform haplotype analysis, or a file with a column named "Gene_ID";
#' @param gffFile: a generic feature format file of certain genome. This file defines the start and end position of specific gene, it should contains at least "seqid", "source", "type", "start", "end", "score", "strand", and "attributes";
#' 
#' @return files extracted and stored in the folder created by default
#' @export GeneSpecificSNPannotation
#' @examples GeneSpecificSNPannotation(SNPAnnoFile = "./SNP_variation/", geneModel = "Glyma.18G056600", gffFile = "./Wm82.a2.v1.gene.gff")

GeneSpecificSNPannotation <- function(SNPAnnoFile = NULL, geneModel = NULL, gffFile=NULL, ...){
  library(data.table)
  library(stringr)
  library(dplyr)
  dirPath <- getwd()
  if(is.null(SNPAnnoFile)) stop("'SNPAnnoFile' is required!")
  if(is.null(geneModel)) stop("'geneModel' is required!")
  if(is.null(gffFile)) stop("'gffFile' is required!")
  
  if(!dir.exists(paste0(dirPath, "/", "SNP.Annotation"))){
    dir.create(paste0(dirPath, "/", "SNP.Annotation"))
  }

  if(dir.exists(SNPAnnoFile)){
    annoList = paste0(SNPAnnoFile, list.files(SNPAnnoFile))
  }else{
    annoList = SNPAnnoFile
  }

  if(file.exists(geneModel)){
    geneList = fread(geneModel, header = T)
    names(geneList) = gsub("gene_id|geneID", "Gene_ID", names(geneList))
    if(length(grep("gene_id|geneID|Gene_ID", names(geneList))) == 0){
      stop("'geneModel' file should contain 'Gene_ID' (column)!")
    }
  }else{
    geneList = data.frame(geneModel)
    names(geneList) = "Gene_ID"
  }
  
  geneList <- distinct(geneList)
  if(file.exists(gffFile)){
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
    geneList <- gff.gene[,c(9,1)]
    names(geneList) <- c("Gene_ID","Chrom")
    geneList$Chrom <- gsub("Gm0|Chr0|CHR|GM", "", geneList$Chrom,ignore.case = T)
    geneList$Chrom <- paste0("Gm",geneList$Chrom)
    chrList = unique(geneList$Chrom)
    
  }
  
  # i = 1; j = 1
  for(i in 1:length(chrList)){
    anno = fread(annoList[grep(paste0(chrList[i], ".csv"), annoList)], header = T)
    gene.temp = geneList[which(geneList$Chrom == chrList[i]),]
    for(j in 1:nrow(gene.temp)){
      annos = anno[grep(gene.temp$Gene_ID[j], anno$Gene_ID),]
      annos$Pos <- as.numeric(annos$Pos)
      gene.loc <- data.frame(matrix(nrow = 1, ncol = 5))
      colnames(gene.loc) <- c("Gene_ID", "Chr", "Start", "End", "Formated")
      gene.loc$Gene_ID <- gene.temp$Gene_ID[j]; gene.loc$Chr <- gsub("Gm0|Chr0|CHR|GM", "", annos$Chrom[1],ignore.case = T)
      gene.loc$Start <- min(annos$Pos); gene.loc$End <- max(annos$Pos)
      gene.loc$Formated[1] <- paste0("Chr", gene.loc$Chr[1], ":", gene.loc$Start[1], "..", gene.loc$End[1])
      fwrite(gene.loc, file = paste0(dirPath, "/", "SNP.Annotation/", gene.temp$Gene_ID[j], ".location.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
      sink(paste0(dirPath, "/", "SNP.Annotation/", gene.temp$Gene_ID[j], ".location.txt"), type = "output")
      print(gene.loc)
      sink(type = "output")
      annos <- annos[order(annos$Pos),]
      fwrite(annos, file = paste0(dirPath, "/", "SNP.Annotation/", gene.temp$Gene_ID[j], ".SNP.annotation.in.total.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
      s1 <- grep("nonsynonymous", annos$Effect_type)
      s2 <- grep("synonymous", annos$Effect_type)
      s3 <- grep("splicing", annos$Effect_type)
      s4 <- grep("stopgain", annos$Effect_type)
      s5 <- grep("stoploss", annos$Effect_type)
      s6 <- grep("intergenic", annos$Effect_type)
      s7 <- grep("intronic", annos$Effect_type)
      s8 <- grep("UTR", annos$Effect_type)
      s9 <- grep("upstream", annos$Effect_type)
      s10 <- grep("downstream", annos$Effect_type)
      gene.related <- unique(c(s1,s2,s3,s4,s5))
      intergenic <- s6
      intron.utr <- c(s7, s8)
      upstream <- s9
      downstream <- s10
      annos$Chrom <- gsub("Gm0|Chr0|CHR|GM", "", annos$Chrom,ignore.case = T)
      SNP.anno.gene.related <- annos[gene.related,]
      SNP.anno.intergenic <- annos[intergenic,]
      SNP.anno.intron.utr <- annos[intron.utr,]
      SNP.anno.upstream <- annos[upstream,]
      SNP.anno.downstream <- annos[downstream,]
      SNP.anno.gene.related <- SNP.anno.gene.related[order(SNP.anno.gene.related$Pos),]
      SNP.anno.intergenic <- SNP.anno.intergenic[order(SNP.anno.intergenic$Pos),]
      SNP.anno.intron.utr <- SNP.anno.intron.utr[order(SNP.anno.intron.utr$Pos),]
      SNP.anno.upstream <- SNP.anno.upstream[order(SNP.anno.upstream$Pos),]
      SNP.anno.downstream <- SNP.anno.downstream[order(SNP.anno.downstream$Pos),]
      fwrite(SNP.anno.gene.related, file = paste0(dirPath, "/", "SNP.Annotation/", gene.temp$Gene_ID[j], ".SNP.annotation.gene.related.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
      fwrite(SNP.anno.intergenic, file = paste0(dirPath, "/", "SNP.Annotation/", gene.temp$Gene_ID[j], ".SNP.annotation.intergenic.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
      fwrite(SNP.anno.intron.utr, file = paste0(dirPath, "/", "SNP.Annotation/", gene.temp$Gene_ID[j], ".SNP.annotation.intron.utr.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
      fwrite(SNP.anno.upstream, file = paste0(dirPath, "/", "SNP.Annotation/", gene.temp$Gene_ID[j], ".SNP.annotation.upstream.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
      fwrite(SNP.anno.downstream, file = paste0(dirPath, "/", "SNP.Annotation/", gene.temp$Gene_ID[j], ".SNP.annotation.downstream.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
    }
  }
}
