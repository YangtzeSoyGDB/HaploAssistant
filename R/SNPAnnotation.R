#' This is some description of this function.
#' @title to annotate where a specific SNP is relative to the nearest gene or within certain region of a gene.
#'
#' @description By using this package, you could use functions of SNPAnnotation to annotate where a specific SNP is relative to the nearest gene or within certain region of a gene.
#'
#' @details see above
#'
#' @param genomeFasta: genomic sequence in fasta format, which is required;
#' @param gffFile: generic feature format file of certain genome, and this file defines the start and end position of specific gene, which should contain at least 9 columns such as "seqid", "source", "type", "start", "end", "score", "strand", "phase" and "attributes".
#' @param hapmapFile: sample segment hapmap file, could be either a file or a folder storing these type of files;
#' @param geneScope: a region to define the upstream and downstream of a gene, genomic position out of this scope would be defined as intergenic. Default value is 2000;
#' @param nthreads: to define how many cores of this computer would be used to process the task. Default is 1. Make sure not to excede the maximum cores of this computer;
#'
#' @return files extracted and stored in the folder created by default.
#' @export SNPAnnotation
#' @examples SNPAnnotation(genomeFasta = "./Gmax_275_v2.0.fa", gffFile = "./Wm82.a2.v1.gene.gff", hapmapFile = "./Hapmap.regional.extracted.files/")

SNPAnnotation = function(genomeFasta = NULL,gffFile = NULL, hapmapFile = NULL, geneScope = 2000, nthreads = 1, ...){
  library(parallel)
  library(Biostrings)
  library(data.table)
  library(stringr)
  library(tidyr)
  library(dplyr)
  dirPath = getwd()

  if(nthreads > 1){
    cores = detectCores()
    if(nthreads <= cores){
      cl.core = makeCluster(getOption("cl.cores", as.numeric(nthreads)))
    }
    if(nthreads > cores) stop(paste0("'nthreads' to be used is out of the maximum core number of this computer: ", cores, ". Please refine the 'nthreads' parameter."))
  }

  if(is.null(genomeFasta)) stop("'genomeFasta' is required!")
  if(is.null(hapmapFile)) stop("'hapmapFile' is required!")
  if(is.null(gffFile)) stop("'gffFile' is required!")


  if(!dir.exists(paste0(dirPath, "/SNP.Annotation/"))){
    dir.create(paste0(dirPath, "/SNP.Annotation/"))
  }

  if(dir.exists(hapmapFile)){
    SNPList = paste0(hapmapFile, list.files(hapmapFile))
  }else{
    SNPList = hapmapFile
  }
  if(length(grep(".log", SNPList,ignore.case = T))>0){
    SNPList = SNPList[-grep(".log", SNPList,ignore.case = T)]
  }

  ## genome file read in
  genome = readBStringSet(genomeFasta, format = "fasta")
  chr = names(genome) %>% as.data.frame()
  names(chr) = "chrName"
  if(length(grep("scaffold", chr$chrName, SNPList,ignore.case = T))>0){
    chr = chr[-grep("scaffold", chr$chrName, SNPList,ignore.case = T),]
  }
  chrPrefix = stringr::str_extract_all(chr, "Chr|Gm")[[1]]

  ## gff file read in
  gff = fread(gffFile, header = T,check.names=F)
  gff$chr = gsub("Chr0|Gm0|Chr|Gm", "", gff$chr,ignore.case = T)

  # i = 9
  snpannotation = function(SNPList, genome, gff, geneScope){
    # function: judge whether a data is hapmap format
    is.hapmap = function(data){
      if(length(grep("chrom|pos|strand|rs#|alleles|assembly#|center|protLSID|assayLSID|panel|QCcode", names(data))) == 11){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }

    # function: judge whether a data is vcf format
    is.vcf = function(data){
      if(length(grep("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT", names(data))) == 9){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }

    # function: identification of the nearest gene of a SNP
    nearestGeneIdentification = function(data = NULL, gff = gff, ...){
      if(is.null(data)) stop("'data' is required!")
      if(is.null(gff)) stop("'gff' is required!")
      gff.temp1 = gff[which(gff$chr == data$Chrom[1] & gff$start <= data$Pos[1]),]
      distance1 = min(abs(gff.temp1$start - data$Pos[1]))
      gff.temp2 = gff[which(gff$chr == data$Chrom[1] & gff$end >= data$Pos[1]),]
      distance2 = min(abs(gff.temp2$end - data$Pos[1]))
      if(distance1 >= distance2){
        gff.temp3 = gff.temp2[which(gff.temp2$end==(data$Pos[1]+distance2)),]
      }else{
        gff.temp3 = gff.temp1[which(gff.temp1$start==(data$Pos[1]-distance1)),]
      }
      return(gff.temp3)
    }

    ## function: to judage wheater a SNP is within a gene or not
    isInGene = function(data = NULL, gff = gff, geneScope = geneScope){
      if(is.null(data)) stop("'data' is required!")
      if(is.null(gff)) stop("'gff' is required!")
      gff.temp = nearestGeneIdentification(data = data, gff = gff)
      if(data$Pos[1] >= min(gff.temp$start) - geneScope & data$Pos[1] <= max(gff.temp$end) + geneScope){
        isingene = TRUE
      }else{
        isingene = FALSE
      }
      return(isingene)
    }

    # function: define the relative location of a SNP to a gene, and give
    AminoAcidChange = function(data = NULL, gffg = gffcds, ...){
      if(is.null(data)) stop("'data' is required!")
      if(is.null(gffg)) stop("'gffg' is required!")
      ref = data$Ref[1]; alt = data$Alt[1]; pos = data$Pos[1];chr = data$Chrom[1]
      gffg = gffg[order(gffg$start),]
      altpos = 0
      for(i in 1:nrow(gffg)){
        if(data$Pos[1] >= gffg$start[i] & data$Pos[1] <= gffg$end[i]){
          altpos = altpos + data$Pos[1] - gffg$start[i] + 1
        }
        if(data$Pos[1] > gffg$start[i] & data$Pos[1] > gffg$end[i]){
          altpos = altpos + gffg$end[i] - gffg$start[i] + 1
        }
      }
      aminopos = ceiling(altpos/3)
      refseq = c();altseq = c()
      for(i in 1:nrow(gffg)){
        temp = toString(Biostrings::subseq(genome[paste0(chrPrefix, ifelse(chr < 10, paste0("0", chr), chr))], start = gffg$start[i], end = gffg$end[i]))
        refseq = paste0(refseq, temp)
      }
      refseq = Biostrings::DNAString(refseq)
      altseq = refseq
      Biostrings::subseq(altseq, start = altpos, end = altpos) = Biostrings::DNAString(alt)
      if(gffg$strand[1] == "-"){
        refseq = Biostrings::reverseComplement(refseq)
        altseq = Biostrings::reverseComplement(altseq)
        altpos = length(refseq)-altpos + 1
        aminopos = ceiling(altpos/3)
      }
      aminoref = Biostrings::translate(Biostrings::RNAString(refseq))
      aminoalt = Biostrings::translate(Biostrings::RNAString(altseq))
      refamino = toString(aminoref[aminopos])
      altamino = toString(aminoalt[aminopos])
      if(refamino == altamino){
        type = "Synonymous SNV"
      }
      if(refamino == "*" & altamino != "*"){
        type = "Stoploss"
      }
      if(refamino != "*" & altamino == "*"){
        type = "Stopgain"
      }
      if(refamino != "*" & altamino != "*" & refamino != altamino){
        type = "Nonsynonymous SNV"
      }
      type = paste0(type, ";", ref, altpos, alt, ";", refamino, aminopos, altamino)
      return(type)
    }

    ## function: judge wheater a SNP is or not in the gff defined region
    isInGff = function(data = NULL, gffs = NULL, ...){
      if(is.null(data)) stop("'data' is required!")
      if(is.null(gffs)) stop("'gffs' is required!")
      isingff = "NA"
      for(i in 1:nrow(gffs)){
        if(data$Pos[1] >= gffs$start[i] & data$Pos[1] <= gffs$end[i]){
          isingff = TRUE
        }
      }
      if(isingff == "NA"){
        isingff = FALSE
      }
      return(isingff)
    }

    ## function:
    SNPDefinition = function(data = NULL, gff = gff, geneScope = geneScope,geneName=geneName){
      if(is.null(data)) stop("'data' is required!")
      if(is.null(gff)) stop("'gff' is required!")
      gff.near = nearestGeneIdentification(data = data, gff = gff)
      if(length(grep(geneName,gff.near$attributes))>0){
        gff.near = gff.near[grep(geneName,gff.near$attributes),]
        gff.near = gff[grep(geneName, gff$attributes),]
        gff.near = gff.near[which(gff.near$type != "gene"),]
        gff.near = gff.near[which(gff.near$type != "mRNA"),]
        geneList = unique(gff.near$attributes)
        # i = 1
        type1 = c()
        # i=1
        for(i in 1:length(geneList)){
          gff.temp = gff.near[which(gff.near$attributes == geneList[i]),]
          gff.temp1 = data.frame(matrix(nrow = (nrow(gff.temp)+2), ncol = ncol(gff.temp)),check.names=F)
          names(gff.temp1) = names(gff.temp)
          gff.temp1[1:nrow(gff.temp),] = gff.temp
          gff.temp1$chr = gff.temp$chr[1]
          gff.temp1$source = gff.temp$source[1]
          gff.temp1$score = "."
          gff.temp1$strand = gff.temp$strand[1]
          gff.temp1$phase[(nrow(gff.temp1)-1):nrow(gff.temp1)] = "."
          gff.temp1$type[nrow(gff.temp1)-1] = "upstream"
          gff.temp1$type[nrow(gff.temp1)] = "downstream"
          gff.temp1$type[nrow(gff.temp1)-1] = "upstream"
          gff.temp1$attributes = gff.temp$attributes[1]
          if(gff.temp1$strand[1] == "-"){
            gff.temp1$start[nrow(gff.temp1)-1] = max(gff.temp$end) + 1
            gff.temp1$end[nrow(gff.temp1)-1] = max(gff.temp$end) + geneScope
            gff.temp1$end[nrow(gff.temp1)] = min(gff.temp$start)  - 1
            gff.temp1$start[nrow(gff.temp1)] = min(gff.temp$start) - geneScope
          }else{
            gff.temp1$start[nrow(gff.temp1)] = max(gff.temp$end) + 1
            gff.temp1$end[nrow(gff.temp1)] = max(gff.temp$end) + geneScope
            gff.temp1$end[nrow(gff.temp1)-1] = min(gff.temp$start)  - 1
            gff.temp1$start[nrow(gff.temp1)-1] = min(gff.temp$start) - geneScope
          }
          gffcds = gff.temp1[which(gff.temp1$type == "CDS"),]
          gffo = gff.temp1[which(gff.temp1$type != "CDS"),]
          if(nrow(gffcds) >= 1){
            if(isInGff(data = data, gffs = gffcds)){
              type = AminoAcidChange(data = data, gffg = gffcds)
              type = paste0(geneList[i], "|", type)
            }
          }
          if(nrow(gffo) >= 1){
            if(isInGff(data = data, gffs = gffo)){
              for(j in 1:nrow(gffo)){
                if(data$Pos[1] >= gffo$start[j] & data$Pos[1] <= gffo$end[j]){
                  type = gffo$type[j]
                }
              }
              type = paste0(geneList[i], "|", type)
            }
          }
          if(nrow(gffcds) >= 1 & nrow(gffo) >= 1 & !isInGff(data = data, gffs = gffo) & !isInGff(data = data, gffs = gffcds)) {
            type = "intron"
            type = paste0(geneList[i], "|", type)
          }
          if(length(geneList) >= 2){
            type1 = paste0(type1, "/", type)
          }
          if(length(geneList) == 1){
            type1 = type
          }
        }
        type1= gsub("^/", "", type1)
        return(type1)
      }
      if(length(grep(geneName,gff.near$attributes))>0){
        type1=NULL
        return(type1)
      }
    }


    ## function: intergenic distance definition
    intergenicDistance = function(data = NULL, gff = gff,geneName = geneName, ...){
      if(is.null(data)) stop("'data' is required!")
      if(is.null(gff)) stop("'gff' is required!")
      gff.near = nearestGeneIdentification(data = data, gff = gff)
      if(length(grep(geneName,gff.near$attributes))!=0){
        gff.near = gff.near[grep(geneName,gff.near$attributes),]
        distance = min(abs(data$Pos[1] - min(gff.near$start)),abs(data$Pos[1] - max(gff.near$end)))
        distance = paste0(geneName, "|", "intergenic;", distance, "bp")
        return(distance)
      }else{
        distance = NULL
        return(distance)
      }
    }
    # i=1
    for(i in 1:length(SNPList)){
      snp = data.table::fread(SNPList[i], header = T,check.names=F,sep = "\t")
      snpName = gsub(".vcf|.hmp|.hapmap|.csv|.txt", "", strsplit(SNPList[i], "/")[[1]][length(strsplit(SNPList[i], "/")[[1]])],ignore.case = T)
      geneName = strsplit(unlist(snpName),"_", fixed = TRUE)[[1]][[1]]
      #geneList <- append(geneList,geneName)
      if(is.hapmap(snp)){
        snp = snp[,c(1:4)]
        snp = tidyr::separate(snp, col = "alleles", into = c("Ref", "Alt"), sep = "/")
        names(snp) = c("ID", "Ref", "Alt", "Chrom", "Pos")
        snp = snp[,c("ID", "Chrom", "Pos", "Ref", "Alt")]
        snp$Gene = NA
        snp$Type = NA
        snp$NucleicChange = NA
        snp$AminoAcidChange = NA
      }else{
        if(is.vcf(snp)){
          snp = snp[,c(1:5)]
          names(snp) = c("Chrom", "Pos", "ID", "Ref", "Alt")
          snp = snp[,c("ID", "Chrom", "Pos", "Ref", "Alt")]
          snp$Gene = NA
          snp$Type = NA
          snp$NucleicChange = NA
          snp$AminoAcidChange = NA
        }else{
          if(!is.vcf(snp) & !is.hapmap(snp)){
            collist = c("Chrom|chrom|chr|Chr|chromosome|ID|Marker|Marker_ID|marker|marker_ID|Position|pos|Pos|position|Reference|ref|Ref|reference|Alternative|Alt|alt|alternative")
            if(length(grep(collist, names(snp))) == 5){
              names(snp) = gsub("Marker|Marker_ID|marker|marker_ID", "ID", names(snp))
              names(snp) = gsub("chrom|chr|Chr|chromosome", "Chrom", names(snp))
              names(snp) = gsub("Position|pos|position", "Pos", names(snp))
              names(snp) = gsub("Reference|ref|reference", "Ref", names(snp))
              names(snp) = gsub("Alternative|alt|alternative", "Alt", names(snp))
              snp = snp[,c("ID", "Chrom", "Pos", "Ref", "Alt")]
            }else{
              stop(paste0("Please check the column name of ", SNPList[i], ". It should contain at least 5 column to define 'maker_ID', 'Chrom', 'Pos', 'Ref', and 'Alt'."))
            }
          }
        }
      }
      chrList = unique(snp$Chrom)
      SNPDefi = data.frame(matrix(nrow = 0, ncol = 9),check.names=F)
      names(SNPDefi) = c("ID", "Chrom", "Pos", "Ref", "Alt", "Gene", "Type", "NucleicChange", "AminoAcidChange")
      # k = 493;j = 1
      for(j in 1:length(chrList)){
        if(length(chrList)==1){
          snp1 = snp
        }else snp1 = snp[which(snp$Chrom == as.character(chrList[j]))]
        for(k in 1:nrow(snp1)){
          if(isInGene(data =  snp1[k,], gff = gff, geneScope = geneScope)){
            type = SNPDefinition(data = snp1[k,], gff = gff, geneScope = geneScope,geneName=geneName)
            if(!is.null(type)){
              type = strsplit(type, "/")[[1]]
              snp.temp = data.frame(matrix(nrow = length(type), ncol = ncol(snp)),check.names=F)
              names(snp.temp) = names(snp)
              snp.temp[,c(1:5)] = snp1[k,c(1:5)]
              for(l in 1:length(type)){
                types = strsplit(type[l], "\\|")[[1]]
                typess = strsplit(types[2], ";")[[1]]
                snp.temp$Gene[l] = types[1]
                snp.temp$Type[l] = typess[1]
                if(length(typess) >= 2){
                  snp.temp$NucleicChange[l] = typess[2]
                  snp.temp$AminoAcidChange[l] = typess[3]
                }
              }
              SNPDefi = rbind(SNPDefi, snp.temp)
            }else next
          }
          if(!isInGene(data =  snp1[k,], gff = gff, geneScope = geneScope)){
            type = intergenicDistance(data = snp1[k,], gff = gff, geneName = geneName)
            if(!is.null(type)){
              type = strsplit(type, "\\|")[[1]]
              snp.temp = data.frame(matrix(nrow = 1, ncol = ncol(snp)),check.names=F)
              names(snp.temp) = names(snp)
              snp.temp[,c(1:5)] = snp1[k,c(1:5)]
              snp.temp$Gene[1] = type[1]
              snp.temp$Type[1] = type[2]
              SNPDefi = rbind(SNPDefi, snp.temp)
            }else next
          }
        }
      }
      # SNPDefi1 <- SNPDefi
      SNPDefi$ID = NULL;names(SNPDefi)[c(1:6)] = c("Chrom","Pos","REF","ALT","Gene_ID","Effect_type")
      SNPDefi <- SNPDefi[, c("Gene_ID","Chrom","Pos","REF","ALT","Effect_type","NucleicChange","AminoAcidChange")]
      #SNPDefi <- SNPDefi1
      if(length(unique(SNPDefi$Gene_ID))>1){
        SNPDefiGeno = SNPDefi[which(SNPDefi$Gene_ID==geneName),]
        SNPDefiGeno1 = SNPDefi[-which(SNPDefi$Gene_ID==geneName),]
        SNPDefi = SNPDefiGeno1[which(SNPDefiGeno1$Gene_ID==sort(unique(SNPDefiGeno1$Gene_ID))[1]),]
        SNPDefi <- dplyr::distinct(SNPDefi)
      }
      SNPDefi$Gene_ID = geneName
      s1 <- grep("nonsynonymous", SNPDefi$Effect_type,ignore.case = T)
      s2 <- grep("synonymous", SNPDefi$Effect_type,ignore.case = T)
      s3 <- grep("splicing", SNPDefi$Effect_type,ignore.case = T)
      s4 <- grep("stopgain", SNPDefi$Effect_type,ignore.case = T)
      s5 <- grep("stoploss", SNPDefi$Effect_type,ignore.case = T)
      s6 <- grep("intergenic", SNPDefi$Effect_type,ignore.case = T)
      s7 <- grep("intronic", SNPDefi$Effect_type,ignore.case = T)
      s8 <- grep("UTR", SNPDefi$Effect_type,ignore.case = T)
      s9 <- grep("upstream", SNPDefi$Effect_type,ignore.case = T)
      s10 <- grep("downstream", SNPDefi$Effect_type,ignore.case = T)
      gene.related <- unique(c(s1,s2,s3,s4,s5))
      intergenic <- s6
      intron.utr <- c(s7, s8)
      upstream <- s9
      downstream <- s10
      SNP.anno.gene.related <- SNPDefi[gene.related,]
      SNP.anno.gene.related$Description <- paste0("c.",SNP.anno.gene.related$NucleicChange,";","p.",SNP.anno.gene.related$AminoAcidChange)
      SNP.anno.gene.related$NucleicChange = NULL;SNP.anno.gene.related$AminoAcidChange = NULL
      SNPDefi$NucleicChange <- NULL;names(SNPDefi)[names(SNPDefi)=="AminoAcidChange"] <- "Description"
      SNP.anno.intergenic <- SNPDefi[intergenic,]
      SNP.anno.intergenic$Description <- stringr::str_split_fixed(SNP.anno.intergenic$Effect_type,";", 2)[,2]
      SNP.anno.intergenic$Effect_type <- stringr::str_split_fixed(SNP.anno.intergenic$Effect_type,";", 2)[,1]
      SNP.anno.intron.utr <- SNPDefi[intron.utr,]
      SNP.anno.upstream <- SNPDefi[upstream,]
      SNP.anno.downstream <- SNPDefi[downstream,]
      SNP.anno.gene.related <- SNP.anno.gene.related[order(SNP.anno.gene.related$Pos),]
      SNP.anno.intergenic <- SNP.anno.intergenic[order(SNP.anno.intergenic$Pos),]
      SNP.anno.intron.utr <- SNP.anno.intron.utr[order(SNP.anno.intron.utr$Pos),]
      SNP.anno.upstream <- SNP.anno.upstream[order(SNP.anno.upstream$Pos),]
      SNP.anno.downstream <- SNP.anno.downstream[order(SNP.anno.downstream$Pos),]
      SNP.anno.list <- list(SNP.anno.gene.related,SNP.anno.intergenic,SNP.anno.intron.utr,SNP.anno.upstream,SNP.anno.downstream)
      SNP.anno.in.total <- data.table::rbindlist(SNP.anno.list)
      SNP.anno.in.total <- SNP.anno.in.total[order(SNP.anno.in.total$Pos),]
      data.table::fwrite(SNP.anno.gene.related, file = paste0(dirPath, "/", "SNP.Annotation/", geneName, ".SNP.annotation.gene.related.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
      data.table::fwrite(SNP.anno.intergenic, file = paste0(dirPath, "/", "SNP.Annotation/", geneName, ".SNP.annotation.intergenic.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
      data.table::fwrite(SNP.anno.intron.utr, file = paste0(dirPath, "/", "SNP.Annotation/", geneName, ".SNP.annotation.intron.utr.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
      data.table::fwrite(SNP.anno.upstream, file = paste0(dirPath, "/", "SNP.Annotation/", geneName, ".SNP.annotation.upstream.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
      data.table::fwrite(SNP.anno.downstream, file = paste0(dirPath, "/", "SNP.Annotation/", geneName, ".SNP.annotation.downstream.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
      data.table::fwrite(SNP.anno.in.total, file = paste0(dirPath, "/", "SNP.Annotation/", geneName, ".SNP.annotation.in.total.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
    }
  }
  if(nthreads == 1){
    snpannotation(SNPList, genome, gff, geneScope)
  }
  if(nthreads > 1){
    parLapply(cl.core, SNPList, snpannotation, genome, gff, geneScope)
    stopCluster(cl.core)
  }
}
