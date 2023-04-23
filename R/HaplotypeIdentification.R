#' This is some description of this function.
#' @title to identify haplotypes from fasta file, and perform basic statistics of each haplotype
#'
#' @description By using this package, you could use functions of haplotype.identification to identify haplotypes from fasta file, and perform basic statistics of each haplotype.
#'
#' @details see above
#'
#' @param fastaFile: input fasta files;
#' @param categoryFile: file stores category of each sample;
#' @param phenoFile: trait value of each sample;
#' @param plotMargin: adjust graph edge margins, default is 2;
#' @param relativeLoc: intron.utr; intergenic; upstream; downstream; gene.related;
#' @param nthreads: to define how many cores of this computer you would like to run the task. Make sure not to excede the maximum of this computer. Default value is 1.

#' @return files extracted and stored in the folder created by default, files of pdf, csv
#' @export HaplotypeIdentification
#' @examples HaplotypeIdentification(fastaFile = "./Geno2Fasta/",categoryFile = "./Trait.category.csv",phenoFile = "./Trait.csv")

HaplotypeIdentification = function(fastaFile = NULL, categoryFile = NULL, phenoFile = NULL, plotMargin = 2, pdfWidth = 5, pdfHeight = 5, nthreads = 1, ...){
  library(ggplot2)
  library(ggpubr)
  library(ggseqlogo)
  library(dplyr)
  library(Biostrings)
  library(data.table)
  library(parallel)
  library(stringr)
  library(tidyr)
  
  if(nthreads > 1){
    cores = detectCores()
    if(nthreads <= cores){
      cl.core = makeCluster(getOption("cl.cores", as.numeric(nthreads)))
    }
    if(nthreads > cores) stop(paste0("'nthreads' to be used is out of the maximum core number of this computer: ", cores, ". Please refine the 'nthreads' parameter."))
  }
  
  # t = 1
  haplotype.identification = function(fastaFile, categoryFile, phenoFile, pdfWidth, pdfHeight, plotMargin){
    ### N replaced by reference
    NReplacedByRef = function(data = NULL, from = 3){
      if(is.null(data)) stop("'data' is required.")
      for(i in 2:nrow(data)){
        for(j in from:ncol(data)){
          if(data[i,j] == "N" | data[i,j] == "-"){
            data[i,j] = data[grep("REF", data[,1]), j]
          }
        }
      }
      return(data)
    }
    
    informativeSiteIdentification = function(dataSep, from = 3){
      informative.site = c()
      # i = 3
      for(i in from:ncol(dataSep)){
        df = unique(dataSep[,i])
        w = grep("W", df)
        r = grep("R", df)
        s = grep("S", df)
        m = grep("M", df)
        y = grep("Y", df)
        k = grep("K", df)
        na = grep("-", df)
        sx = c(w,r,s,m,y,k,na)
        if(length(sx) == 0){df = df}else{df = df[-sx]}
        if(length(df) > 1){informative.site = c(informative.site, (i-from+1)) }
      }
      return(informative.site)
    }
    
    dirPath = getwd()
    if(is.null(fastaFile)) stop("'fastaFile' is required!")
    
    if(!dir.exists(paste0(dirPath, "/", "Haplotype.Analysis"))){
      dir.create(paste0(dirPath, "/", "Haplotype.Analysis"))
    }
    if(dir.exists(fastaFile)){
      fastaList = paste0(fastaFile, list.files(fastaFile))
    }else{
      fastaList = fastaFile
    }
    annoList = fastaList[grep(".csv", fastaList)]
    fastaList = fastaList[grep(".fasta", fastaList)]
    
    if(!is.null(categoryFile)){
      category = data.table::fread(categoryFile, header = T)
      names(category) = gsub("type", "Type", names(category))
      categoryList = colnames(category)[2:ncol(category)]
      categoryList = categoryList[-grep("longitude|latitude", categoryList,ignore.case = T)]
      names(category)[1] = "Sample_ID"
      categoryN = length(categoryList)
    }else{stop("'categoryFile' is requied!")}
    
    if(!is.null(phenoFile)){
      pheno = data.table::fread(phenoFile, header = T)
    }else{
      pheno <- data.frame(category$Sample_ID)
      pheno$trait <- rnorm(nrow(pheno))
    }
    names(pheno)[1] = "Sample_ID"
    phenoN = ncol(pheno)-1
    phenoList = colnames(pheno)[2:ncol(pheno)]
    
    if(!is.null(categoryFile)){
      categoryPheno = merge(x = category, y = pheno, by = "Sample_ID", all.x = TRUE)
    }
    
    geneList = c()
    for(i in 1:length(fastaList)){
      gene = strsplit(strsplit(fastaList[i], "/")[[1]][length(strsplit(fastaList[i], "/")[[1]])], "_")[[1]][1]
      geneList = c(geneList, gene)
    }
    geneList = unique(geneList)
    
    # t = 2
    for(t in 1:length(geneList)){
      if(!dir.exists(paste0(dirPath, "/Haplotype.Analysis/", as.character(geneList[t])))){
        dir.create(paste0(dirPath, "/Haplotype.Analysis/", as.character(geneList[t])))
      }
      sink(paste0(dirPath, "/Haplotype.Analysis/", geneList[t], "/", geneList[t], ".haplotype.identification.log"), type = "output", append = TRUE)
      writeLines(paste0("\n####### haplotype identification of ", geneList[t], " ######\n"))
      writeLines("\tCountry\tTrait\tCategory\tHaploSeq\tSeqN")
      sink(type = "output")
      
      fastaList.temp = fastaList[grep(geneList[t], fastaList)]
      annoList.temp = annoList[grep(geneList[t], annoList)]
      
      # a = 1
      for(a in 1:length(fastaList.temp)){
        if(length(grep("\\.fasta", fastaList.temp[a])) != 0){
          relativeLoc = gsub(".fasta", "", strsplit(fastaList.temp[a], "/")[[1]][length(strsplit(fastaList.temp[a], "/")[[1]])])
          relativeLoc = stringr::str_extract_all(relativeLoc, "intergenic|downstream|upstream|gene.related|intron.utr|in.total")[[1]]
          for(d in 1:length(annoList.temp)){
            if(length(grep(relativeLoc, annoList.temp[d])) == 1){
              anno = data.frame(data.table::fread(annoList.temp[d], header = T))
              anno = anno[,-c(8:ncol(anno))]
            }
          }
          fas = data.frame(Biostrings::readDNAStringSet(fastaList.temp[a], format = "fasta"))
          data = data.frame(row.names(fas), fas[,1])
          colnames(data) = c("Sample_ID", "seq")
          if(nrow(data) > 0){
            ##### informative site selection
            seqLength = length(strsplit(as.character(data$seq[1]), "")[[1]])
            dataSep = data.frame(matrix(nrow = nrow(data), ncol = ncol(data) + seqLength))
            colnames(dataSep) = c(colnames(data), paste0("n",seq(1:seqLength)))
            dataSep[,1:2] = data
            for(i in 1:nrow(dataSep)){
              dataSep[i,c(3:ncol(dataSep))] = strsplit(as.character(dataSep$seq[i]), "")[[1]]
            }
            dataSep = NReplacedByRef(data = dataSep, from = 3)
            informative.site = informativeSiteIdentification(dataSep = dataSep, from = 3)
            if(length(informative.site) >= 1){
              data.selected = dataSep[,c(1,2,(informative.site+2))]
              anno.selected = anno[informative.site,]
              for(i in 1:nrow(data.selected)){
                temp1 = data.selected[i,3]
                for(j in 4:ncol(data.selected)){
                  temp1 = paste0(temp1, data.selected[i,j])
                }
                data.selected[i,2] = temp1
              }
              data.selected = data.selected[,-c(3:ncol(data.selected))]
              data.total = merge(x = data.selected, y =  categoryPheno, by = "Sample_ID", all = TRUE)
              ref.data = data.total[which(data.total[,1] == "REF"),]
              ref.seq = strsplit(ref.data$seq[1], "")[[1]]
              
              # i = 1
              #
              for(i in 1:nrow(data.total)){
                temp1 = strsplit(data.total$seq[i], "")[[1]]
                temp1[grep("-", temp1)] = ref.seq[grep("-", temp1)]
                temp2 = temp1[1]
                for(j in 2:length(temp1)){
                  temp2 = paste0(temp2, temp1[j])
                }
                data.total$seq[i] = temp2
              }
              
              # i = 1; j =1
              for(i in 1:categoryN){
                catcol = grep(categoryList[i], names(data.total))
                longcol = grep("longitude|Longitude", names(data.total))
                latcol = grep("latitude|Latitude", names(data.total))
                for(j in 1:phenoN){
                  phenocol = grep(phenoList[j], names(data.total))
                  temp1 = data.total[,c(1,2,catcol, longcol, latcol, phenocol)]
                  temp1 = temp1[which(temp1[,ncol(temp1)] != ""),]
                  temp1 = temp1[which(temp1[,ncol(temp1)] != "NA"),]
                  temp1 = temp1[which(temp1[,ncol(temp1)] != "na"),]
                  temp1 = temp1[which(temp1[,ncol(temp1)] != "-999"),]
                  temp1col = ncol(temp1)
                  ##### informative site selection
                  seqLength = length(strsplit(as.character(temp1$seq[1]), "")[[1]])
                  dataSep = data.frame(matrix(nrow = nrow(temp1), ncol = temp1col + seqLength))
                  colnames(dataSep) = c(names(temp1), paste0("n",seq(1:seqLength)))
                  dataSep[,1:temp1col] = temp1
                  for(c in 1:nrow(dataSep)){
                    if(length(grep("NA", dataSep$seq[c])) == 0){
                      dataSep[c,c((temp1col+1):ncol(dataSep))] = strsplit(as.character(dataSep$seq[c]), "")[[1]]
                    }
                  }
                  informative.site = informativeSiteIdentification(dataSep = dataSep, from = 7)
                  anno.selected.1 = anno.selected[informative.site,]
                  # ref.pattern
                  ref.seq = ref.seq[informative.site]
                  ref.pattern = ref.seq[1]
                  for(g in 2:length(ref.seq)){
                    ref.pattern = paste0(ref.pattern, ref.seq[g])
                  }
                  data.selected = dataSep[,c(1:temp1col,(informative.site+temp1col))]
                  for(e in 1:nrow(data.selected)){
                    temp4 = data.selected[e,temp1col+1]
                    for(f in (temp1col+2):ncol(data.selected)){
                      temp4 = paste0(temp4, data.selected[e,f])
                    }
                    data.selected[e,2] = temp4
                  }
                  data.selected = data.selected[,-c((temp1col+1):ncol(data.selected))]
                  temp1 = data.selected
                  temp1$haplotype = 0
                  haploFreq = data.frame(table(temp1$seq))
                  haploFreq = haploFreq[order(haploFreq$Freq, decreasing = T),]
                  # haplotype identification
                  haplo.list = unique(unlist(haploFreq$Var1))
                  haplo.list.1 = c()
                  for(h in 1:length(haplo.list)){
                    temp2 =  strsplit(as.character(haplo.list[h]), "")[[1]]
                    w = grep("W", temp2)
                    r = grep("R", temp2)
                    s = grep("S", temp2)
                    m = grep("M", temp2)
                    y = grep("Y", temp2)
                    k = grep("K", temp2)
                    heterozygous = c(w,r,s,m,y,k)
                    if(length(heterozygous) == 0) {
                      haplo.list.1 = c(haplo.list.1, haplo.list[h])
                    }
                  }
                  if(!is.null(haplo.list.1)){
                    haplo.list.1 = unique(haplo.list.1)
                    haplo.list = haplo.list.1
                    for(q in 1:nrow(temp1)){
                      seq.q = strsplit(as.character(temp1$seq[q]), "")[[1]]
                      w = grep("W", seq.q)
                      r = grep("R", seq.q)
                      s = grep("S", seq.q)
                      m = grep("M", seq.q)
                      y = grep("Y", seq.q)
                      k = grep("K", seq.q)
                      heterozygous = c(w,r,s,m,y,k)
                      if(length(heterozygous) >= 1){
                        temp1$haplotype[q] = c("deleted")
                      }
                    }
                    temp1 = temp1[which(temp1$haplotype != "deleted"),]
                    ######
                    ##### informative site selection
                    seqLength = length(strsplit(as.character(temp1$seq[1]), "")[[1]])
                    dataSep = data.frame(matrix(nrow = nrow(temp1), ncol = ncol(temp1) + seqLength))
                    colnames(dataSep) = c(colnames(temp1), paste0("n",seq(1:seqLength)))
                    dataSep[,1:ncol(temp1)] = temp1
                    for(c in 1:nrow(dataSep)){
                      dataSep[c,c((ncol(temp1)+1):ncol(dataSep))] = strsplit(as.character(dataSep$seq[c]), "")[[1]]
                    }
                    informative.site = informativeSiteIdentification(dataSep = dataSep, from = 8)
                    anno.selected.2 = anno.selected.1[informative.site,]
                    # ref.pattern
                    ref.seq = ref.seq[informative.site]
                    ref.pattern = ref.seq[1]
                    for(g in 2:length(ref.seq)){
                      ref.pattern = paste0(ref.pattern, ref.seq[g])
                    }
                    data.selected = dataSep[,c(1:ncol(temp1),informative.site+ncol(temp1))]
                    for(e in 1:nrow(data.selected)){
                      temp4 = data.selected[e,ncol(temp1)+1]
                      for(f in (ncol(temp1)+2):ncol(data.selected)){
                        temp4 = paste0(temp4, data.selected[e,f])
                      }
                      data.selected[e,2] = temp4
                    }
                    data.selected = data.selected[,-c((ncol(temp1)+1):ncol(data.selected))]
                    temp1 = data.selected
                    ######
                    pdf(paste0(dirPath, "/Haplotype.Analysis/", geneList[t], "/", geneList[t], "_", phenoList[j], "_", relativeLoc, ".seq.logo.pdf"), width = ifelse(nchar(temp1$seq[1]) <= 10, pdfWidth, nchar(temp1$seq[1])*pdfWidth/10), height = pdfHeight)
                    p1 = ggseqlogo::ggseqlogo(as.character(temp1$seq))
                    p1 = p1 + ggplot2::theme(plot.margin=ggplot2::unit(rep(plotMargin, 4),'lines'))
                    print(p1)
                    dev.off()
                    ###### haplotype identification
                    haploFreq = data.frame(table(temp1$seq))
                    haploFreq = haploFreq[order(haploFreq$Freq, decreasing = T),]
                    haplo.list = unique(haploFreq$Var1)
                    hap.freq = data.frame(matrix(nrow = length(haplo.list), ncol = 10))
                    colnames(hap.freq) = c("Haplotype", "Country", "Pattern", "Frequency", "Min", "1st.Qu", "Median", "Mean", "3rd.Qu", "Max")
                    hap.freq$Haplotype = paste0("Hap", seq(1:nrow(hap.freq)))
                    hap.freq$Pattern = haplo.list
                    for(l in 1:nrow(hap.freq)){
                      if(hap.freq$Pattern[l] == ref.pattern){
                        hap.freq$Class[l] = c("reference")
                      }else{
                        hap.freq$Class[l] = c("specific")
                      }
                    }
                    # b = 1; d =1
                    for(b in 1:nrow(hap.freq)){
                      for(d in 1:nrow(temp1)){
                        if(temp1$seq[d] == hap.freq$Pattern[b]){
                          temp1$haplotype[d] = hap.freq$Haplotype[b]
                        }
                      }
                      temp2 = temp1[which(temp1$haplotype == hap.freq$Haplotype[b]),]
                      hap.freq$Frequency[b] = nrow(temp2)
                      summary = data.frame(as.character(summary(temp2[,grep(phenoList[j], names(temp2))])))
                      summary.t = t(summary)
                      hap.freq[b,c(5:10)] = summary.t[1,]
                    }
                    data.table::fwrite(temp1, file = paste0(dirPath, "/Haplotype.Analysis/", geneList[t], "/", geneList[t], "_", phenoList[j], "_", categoryList[i], "_", relativeLoc, ".haplotype"), row.names = F, col.names = T, quote = F, sep = "\t")
                    data.table::fwrite(hap.freq, file = paste0(dirPath, "/Haplotype.Analysis/", geneList[t], "/", geneList[t], "_", phenoList[j], "_", categoryList[i], "_", relativeLoc, ".haplotype.frequency.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
                    
                    ########################
                    ### haplotype categorization
                    class = unique(temp1[,3])
                    haplo.type = data.frame(matrix(nrow = nrow(hap.freq), ncol = length(class)+3))
                    colnames(haplo.type) = c("Haplotype","Pattern", "Class", class)
                    haplo.type$Haplotype = hap.freq$Haplotype
                    haplo.type$Pattern = hap.freq$Pattern
                    haplo.type$Class = hap.freq$Class
                    o = 1; p = 1
                    for(o in 1:nrow(haplo.type)){
                      temp3 = temp1[which(temp1$haplotype == haplo.type$Haplotype[o]), ]
                      temp3 = temp3[which(temp3[,4] != "-999"),]
                      for(p in 1:length(class)){
                        haplo.type[o,p+3] = length(grep(class[p], temp3[,3]))
                      }
                    }
                    data.table::fwrite(haplo.type, file = paste0(dirPath, "/Haplotype.Analysis/", geneList[t], "/", geneList[t], "_", phenoList[j], "_", categoryList[i], "_", relativeLoc, ".haplotype.category.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
                    data.table::fwrite(anno.selected.2, file = paste0(dirPath, "/Haplotype.Analysis/", geneList[t], "/", geneList[t], "_", phenoList[j], "_", categoryList[i], "_", relativeLoc, ".haplotype.SNP.annotation.csv"), row.names = F, col.names = T, quote = F, sep = "\t")
                    sink(paste0(dirPath, "/Haplotype.Analysis/", geneList[t], "/", geneList[t], ".haplotype.identification.log"), type = "output", append = TRUE)
                    writeLines(paste0("\t",relativeLoc, "\t", phenoList[j], "\t", categoryList[i], "\t", hap.freq$Pattern[grep("reference", hap.freq$Class)], "\t", nchar(as.character(hap.freq$Pattern[1]))))
                    sink(type = "output")
                  }
                  if(is.null(haplo.list.1)){
                    sink(paste0(dirPath, "/Haplotype.Analysis/", geneList[t], "/", geneList[t], ".haplotype.identification.log"), type = "output", append = TRUE)
                    writeLines(paste0("\t",relativeLoc, "\t", phenoList[j], "\t", categoryList[i], "\tNA\t", 0))
                    sink(type = "output")
                  }
                }
              }
            }
            
          }
        }
      }
      sink(paste0(dirPath, "/Haplotype.Analysis/", geneList[t], "/", geneList[t], ".haplotype.identification.log"), type = "output", append = TRUE)
      writeLines("\n ###### THE END ######\n")
      sink(type = "output")
    }
    
  }
  
  if(nthreads == 1){
    haplotype.identification(fastaFile, categoryFile, phenoFile, pdfHeight, pdfWidth, plotMargin)
  }
  if(nthreads > 1){
    parLapply(cl.core, fastaFile, haplotype.identification, categoryFile, phenoFile, pdfHeight, pdfWidth, plotMargin)
    stopCluster(cl.core)
  }
  
  
}