

#' This is some description of this function.
#' @title to compare phenotype difference between haplotypes
#'
#' @description By using this package, you could use function of HaplotypePlotting to compare phenotype difference between haplotypes
#'
#' @details see above
#'
#' @param haploResults: haplotype analysis results, could be either a folder stores gene list files or a single file of gene list;
#' @param haploMin: the minimum samples of a specific haplotype taking into account;
#' @param plotMargin: adjust graph edge margins, default is 2;
#' @param nthreads: to define how many cores of this computer you would like to run the task. Make sure not to excede the maximum of this computer. Default value is 1.

#' @return files of pdf
#' @export HaplotypePlotting
#' @examples HaplotypePlotting(haploResults = "./Haplotype.Analysis/", nthreads = 6)

HaplotypePlotting = function(haploResults = NULL, haploMin = 10, pdfWidth = 10, pdfHeight = 10, plotMargin = 2, nthreads = 1, ...){
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(ggpubr)
  library(multcomp)
  library(parallel)
  library(stringr)

  if(nthreads > 1){
    cores = detectCores()
    if(nthreads <= cores){
      cl.core = makeCluster(getOption("cl.cores", as.numeric(nthreads)))
    }
    if(nthreads > cores) stop(paste0("'nthreads' to be used is out of the maximum core number of this computer: ", cores, ". Please refine the 'nthreads' parameter."))
  }

  haploplotting = function(haploResults, haploMin = haploMin, pdfWidth = pdfWidth, pdfHeight = pdfHeight, plotMargin = plotMargin){
    dirPath = getwd()
    if(!dir.exists(paste0(dirPath, "/Haplo_pheno_test"))){
      dir.create(paste0(dirPath, "/Haplo_pheno_test"))
    }
    if(is.null(haploResults)) stop("'haploResults' is required")

    if(dir.exists(haploResults)){
      haploList = paste0(haploResults, list.files(haploResults))
    }else{
      haploList = haploResults
    }

    genelist = data.frame(stringr::str_split(as.list(haploList), "/"))[]
    genelist = as.character(genelist[nrow(genelist),])
    ## comparison function definition
    mycomparison = function(data){
      df = unique(data)
      if(length(df) >= 3){
        list.n = combn(length(df),2)
        compare.list = list(list.n[,])[[1]]
        my.comp = list()
        for(i in 1:ncol(compare.list)){
          temp1 <- c(df[compare.list[1,i]], df[compare.list[2,i]])
          my.comp <- c(my.comp, list(temp1))
        }
      }
      if(length(df) == 2){
        my.comp = list(c(df[1], df[2]))
      }
      if(length(df) == 1){
        my.comp = NULL
      }
      return(my.comp)
    }
    # i=2;j=5
    for(i in 1:length(genelist)){
      if(!dir.exists(paste0(dirPath, "/Haplo_pheno_test/", genelist[i]))){
        dir.create(paste0(dirPath, "/Haplo_pheno_test/", genelist[i]))
      }
      haploList.temp = haploList[grep(genelist[i], haploList)]
      if(dir.exists(haploList.temp)){
        haploList.temp = paste0(haploList.temp, "/", list.files(haploList.temp))
        haploList.temp = haploList.temp[-grep(".csv|.pdf|.log", haploList.temp)]
      }
      for(j in 1:length(haploList.temp)){
        fileName = strsplit(haploList.temp[j], "/")[[1]][length(strsplit(haploList.temp[j], "/")[[1]])]
        fileName = gsub("haplotype","",fileName,ignore.case = T)
        haplo = data.table::fread(haploList.temp[j], header = T)
        data.table::setkeyv(haplo, names(haplo))
        hap.freq = data.frame(table(haplo$haplotype))
        names(hap.freq) = c("haplotype", "freq")
        if(length(which(hap.freq$freq > haploMin))==0) print(paste0(fileName,"'s freq are less than haploMin."))
        if(length(which(hap.freq$freq > haploMin))==0) next
        hap.freq = hap.freq[which(hap.freq$freq > haploMin),]
        haplos = haplo[which(haplo$haplotype %in% hap.freq$haplotype),]
        haplos$haplotype = gsub("Hap","",haplos$haplotype)
        haplos$haplotype = as.numeric(haplos$haplotype)
        haplos = haplos[order(c(haplos$haplotype)),]
        haplos$haplotype = paste0("Hap",haplos$haplotype)
        haplo.list = sort(unique(haplos$haplotype))
        my_comparisons = mycomparison(data = haplo.list)
        if(!is.null(my_comparisons)){
          pdf(paste0(dirPath, "/Haplo_pheno_test/", genelist[i], "/", fileName, "haplo_pheno_test.pdf"), width = pdfWidth, height = pdfHeight)
          p1 = ggpubr::ggviolin(haplos, x="haplotype", y= names(haplos)[ncol(haplos)-1], fill = "haplotype", add = "boxplot", add.params = list(fill="white"))+
            ggpubr::stat_compare_means(comparisons = my_comparisons)
          p1 = p1 + ggplot2::theme(plot.margin=ggplot2::unit(rep(plotMargin,4),'lines'))
          #
          print(p1)
          dev.off()
        }
        if(is.null(my_comparisons)){
          pdf(paste0(dirPath, "/Haplo_pheno_test/", genelist[i], "/", fileName, "haplo_pheno_test.pdf"), width = pdfWidth, height = pdfHeight)
          p2 = ggpubr::ggviolin(haplos, x="haplotype", y= names(haplos)[ncol(haplos)-1], fill = "haplotype")
          p2 = p2 + ggplot2::theme(plot.margin=ggplot2::unit(rep(plotMargin,4),'lines'))
          print(p2)
          dev.off()
        }
      }
    }
  }

  if(nthreads == 1){
    haploplotting(haploResults, haploMin, pdfWidth, pdfHeight, plotMargin)
  }
  if(nthreads > 1){
    parLapply(cl.core, haploResults, haploplotting, haploMin, pdfWidth, pdfHeight, plotMargin)
    stopCluster(cl.core)
  }
}
