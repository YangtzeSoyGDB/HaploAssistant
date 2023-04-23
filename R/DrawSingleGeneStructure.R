#' This is some description of this function.
#' @title to draw gene structure of each single gene.
#'
#' @description By using this package, you could use the function of DrawSingleGeneStructure to draw gene structure of each single gene.
#'
#' @details see above
#'
#' @param gffFile: generic feature format file of certain genome, and this file defines the start and end position of specific gene, which should contain at least 9 columns such as "seqid", "source", "type", "start", "end", "score", "strand", "phase" and "attributes";
#' @param haploResults: haplotype analysis results, could be either a folder stores gene list files or a single file of gene list;
#' @param size: height of gene body, default value is 5;
#' @param lcolor: color of the middle line, default value is "red".
#' 
#' @return files extracted and stored in the folder created by default.
#' @export DrawSingleGeneStructure
#' @examples DrawSingleGeneStructure(gffFile = "./Wm82.a2.v1.gene.gff", haploResults = "./Haplotype.Analysis/",size = 5, lcolor= "red")

DrawSingleGeneStructure = function(gffFile = NULL, haploResults = NULL, size = 5, lcolor= "red", pdfWidth = 16, pdfHeight = 4, ...){
  library(gggenes)
  library(ggplot2)
  library(data.table)
  library(ggrepel)
  library(stringr)
  dirPath <- getwd()
  typelist = c("downstream|upstream|intron.utr|gene.related|in.total|intergenic")
  if(is.null(gffFile)) stop("'gffFile' is required!")
  if(is.null(haploResults)) stop("'haploResults' is required!")
  if(is.null(size)) size = 5
  if(is.null(lcolor)) lcolor = "red"

  if(!dir.exists(paste0(dirPath, "/", "Gene.Structure"))){
    dir.create(paste0(dirPath, "/", "Gene.Structure"))
  }

  if(dir.exists(haploResults)){
    annoList = paste0(haploResults, list.files(haploResults))
  }else{
    annoList = haploResults
  }
  genelist = data.frame(str_split(as.list(annoList), "/"))[]
  genelist = as.character(genelist[nrow(genelist),])

  gff <- data.frame(fread(file =gffFile, header = T, sep = "\t"))
  names(gff) = c("chr","source","type","start","end","score","strand","phase","attributes")

  # i = 6;j = 1; k = 1
  for(i in 1:length(genelist)){
    gene.gff = gff[grep(genelist[i], gff$attributes),]
    annoList.temp = annoList[grep(genelist[i], annoList)]
    if(length(annoList.temp) == 0){
      next
    }
    if(dir.exists(annoList.temp)){
      annoList.temp = paste0(annoList.temp, "/", list.files(annoList.temp))
    }
    if(length(grep("\\.log", annoList.temp))>0){
      annoList.temp = annoList.temp[-grep("\\.log", annoList.temp)]
    }
    annoList.temp = annoList.temp[grep("SNP.annotation.csv", annoList.temp)]
    if(nrow(gene.gff) > 0) {
      if(!dir.exists(paste0(dirPath, "/", "Gene.Structure/", genelist[i]))){
        dir.create(paste0(dirPath, "/", "Gene.Structure/", genelist[i]))
      }
      gene.gff.total <- unique(gene.gff[,9])
      gene.gff.total <- gene.gff.total[-which(gene.gff.total == genelist[i])]
      if(length(annoList.temp) > 0){
        for(j in 1:length(annoList.temp)){
          SNP.anno <- data.frame(fread(file = annoList.temp[j], header = T, sep = "\t"))
          SNP.anno = SNP.anno[!is.na(SNP.anno[,3]),]
          type = str_extract_all(annoList.temp[j], typelist)[[1]]
          fileName = gsub(".haplotype.SNP.annotation.csv", "", strsplit(annoList.temp[j], "/")[[1]][length(strsplit(annoList.temp[j], "/")[[1]])])
          # a = 3
          SNP.anno$nucleic_acid_change = NA
          SNP.anno$amino_acid_change = NA
          for(a in 1:nrow(SNP.anno)){
            if(!is.na(SNP.anno$Description[a])){
              SNP.des <- strsplit(SNP.anno$Description[a], ":|;")[[1]]
              if(length(grep("c\\.", SNP.des)) >= 1){
                nucleic <- strsplit(SNP.des[grep("c.", SNP.des)],",")[[1]]
                SNP.anno$nucleic_acid_change[a] <-  nucleic[grep("c.", nucleic)]
              }
              if(length(grep("p\\.", SNP.des)) >= 1){
                amino <- strsplit(SNP.des[grep("p.", SNP.des)],",")[[1]]
                SNP.anno$amino_acid_change[a] <-  amino[grep("p.", amino)]
              }
            }else next
          }
          SNP.anno$nucleic_acid_change <- gsub("c.", "", SNP.anno$nucleic_acid_change)
          SNP.anno$amino_acid_change <- gsub("p.", "", SNP.anno$amino_acid_change)
          # b = 3
          for(b in 1:nrow(SNP.anno)){
            if(is.na(SNP.anno$nucleic_acid_change[b])){
              SNP.anno$nucleic_acid_change[b] = SNP.anno$Effect_type[b]
            }
            if(is.na(SNP.anno$amino_acid_change[b])){
              SNP.anno$amino_acid_change[b] = SNP.anno$Effect_type[b]
            }
          }
          for(k in 1:length(gene.gff.total)){
            gene.gff.n <- gene.gff[which(gene.gff[,9] == gene.gff.total[k]),]
            gene.gff.n <- gene.gff.n[which(gene.gff.n[,3] != "mRNA"),]
            gene.gff.n.cds <-  gene.gff.n[which(gene.gff.n[,3] == "CDS"),]
            gene.gff.n.utr <- gene.gff.n[which(gene.gff.n[,3] != "CDS"),]
            gene.gff.n.cds$type <- paste0(gene.gff.n.cds$type, seq(1:nrow(gene.gff.n.cds)))
            gene.gff.n <- rbind(gene.gff.n.cds, gene.gff.n.utr)
            gene.gff.n <- gene.gff.n[order(gene.gff.n$start),]
            gene.gff.n$y <- 1
            chr <- unique(gene.gff.n$chr)
            strand <- unique(gene.gff.n$strand)
            gene.symbol <- unique(gene.gff.n$attributes)
            if (strand == "-"){
              three_end = min(gene.gff.n$start)-100
            } else {
              three_end = max(gene.gff.n$end)+100
            }
            if (strand == "+"){
              five_start = min(min(gene.gff.n$start), min(SNP.anno$Pos))-100
            } else {
              five_start = max(max(gene.gff.n$end), max(SNP.anno$Pos))+100
            }
            POS = sort(unique(c(gene.gff.n$start, gene.gff.n$end, SNP.anno$pos)))
            End = max(POS); Start = min(POS)

            p = ggplot(gene.gff.n,aes(x = start,y = y)) +
              labs(title = paste0(gene.symbol, " (", strand, "); ", type)) +
              geom_segment(data = gene.gff.n,
                           mapping = aes(x = Start,y = 0, xend = End, yend = 0),
                           color = lcolor) +
              geom_segment(data = gene.gff.n,
                           mapping = aes(x = end[1],y = 0, xend = three_end, yend = 0),
                           arrow = arrow(length = unit(0.2,"cm")),color = "black") +
              geom_segment(data = gene.gff.n,
                           mapping = aes(x = start[1],y = 0, xend = five_start, yend = 0),
                           color = "black") +
              geom_segment(data = SNP.anno,
                           mapping = aes(x = pos, xend = pos, y = 0, yend = 0.02, color = "black"),
                           color = lcolor) +
              geom_segment(data = gene.gff.n,
                           mapping = aes(x = start,xend = end,y = 0, yend = 0, color = type),
                           size = size) +
              ylim(-0.05, 0.1) + xlab(chr) +
              geom_text_repel(data = SNP.anno,
                              mapping = aes(label = paste0(nucleic_acid_change, "(", amino_acid_change, ")"),
                                            x = pos, y = 0.03),
                              nudge_x=0,
                              nudge_y=0.01,
                              hjust = 0,
                              direction = "x",
                              angle=90,
                              force_pull = 0,
                              size=2.5,
                              segment.inflect=T,
                              max.overlaps=80,segment.linetype=1,
                              segment.color = "grey",
                              max.iter = 1e4, max.time = 1)+
              geom_text_repel(data = SNP.anno,
                              mapping = aes(label = pos,x = pos, y = -0.01),
                              nudge_x=0,
                              nudge_y=-0.01,
                              angle=90,
                              hjust = 1,
                              direction = "x",
                              size=2.5,
                              min.segmet.length = 0,
                              segment.inflect=T,
                              max.overlaps=80,segment.linetype=1,
                              segment.color = "grey",
                              max.iter = 1e4, max.time = 1)+
              scale_x_continuous(
                expand = expansion(mult = 0.1))+
              theme(axis.ticks.y.left = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.line.y = element_blank(),
                    axis.line.x.top = element_blank(),
                    panel.grid = element_blank(),
                    axis.text.y = element_blank(),
                    plot.background = element_blank(),
                    panel.background = element_blank(),
                    axis.title.y = element_blank(),
                    axis.line = element_line(colour = "black",size = 0.5))
            pdf(paste0(dirPath, "/", "Gene.Structure", "/", genelist[i], "/", gene.symbol,"_", type, "_gene.structure.pdf"), height = pdfHeight, width = pdfWidth)
            print(p)
            dev.off()
          }
        }
      }
    }
  }
}
