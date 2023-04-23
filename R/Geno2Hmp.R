#' This is some description of this function.
#' @title transforming the .vcf or the .geno into .hapmap format
#'
#' @description By using this package, you could transform the .vcf or the .geno into .hapmap format.
#'
#' @details see above
#'
#' @param genoFile: the files in vcf format or .geno format you want to transform;
#' @param nthreads: to define how many cores of this computer you would like to run the task. Make sure not to excede the maximum of this computer. Default value is 1.
#' 
#' @return files extracted and stored in the folder created by default.
#' @export Geno2Hmp
#' @examples Geno2Hmp(genoFile = "./vcfFile/",nthreads = 4)

Geno2Hmp <- function(genoFile=NULL,nthreads = 1){
  library(data.table)
  library(parallel)
  
  if(nthreads > 1){
    cores = detectCores()
    if(nthreads <= cores){
      cl.core = makeCluster(getOption("cl.cores", as.numeric(nthreads)))
    }
    if(nthreads > cores) stop(paste0("'nthreads' to be used is out of the maximum core number of this computer: ", cores, ". Please refine the 'nthreads' parameter."))
  }
  geno2hmp = function(genoFile=NULL){
    dirPath = getwd()
    if(!dir.exists(paste0(dirPath, "/geno2hmp/"))){
      dir.create(paste0(dirPath, "/geno2hmp/"))
    }
    
    # function: judge whether a data is vcf format
    is.vcf = function(data){
      if(length(grep("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT", names(data),ignore.case = T)) == 9){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }
    # function: judge whether a data is vcf format
    is.geno = function(data){
      if(length(grep("#CHROM|POS|ID|REF|ALT|ANNO", names(data),ignore.case = T)) == 5){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }
    
    if(dir.exists(genoFile)){
      genoList = paste0(genoFile, list.files(genoFile))
    }else{
      genoList = genoFile
    }
    
    is.geno = function(data){
      if(length(grep("#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT", names(data),ignore.case = T)) == 9){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }
    if(dir.exists(genoFile)){
      genoList = paste0(genoFile, list.files(genoFile))
    }else{
      genoList = genoFile
    }
    
    # i = 1
    for (i in 1:length(genoList)) {
      genoName = gsub(".geno|.csv|.txt|.hmp|.vcf", "", strsplit(genoList[i], "/")[[1]][length(strsplit(genoList[i], "/")[[1]])])
      geno = data.frame(data.table::fread(genoList[i], header = T, fill=TRUE),check.names = F)
      hmp = data.frame(matrix(nrow = nrow(geno.sample), ncol = 11))
      names(hmp) = c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panel","QCcode")
      
      if(is.vcf(data=geno)){
        names(geno)[1:9] <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
        if(length(which(nchar(geno$REF)>1)&(nchar(geno$ALT)>1))>0){
          geno <- geno[-(which(nchar(geno$REF)>1)&(nchar(geno$ALT)>1)),]
        }
        geno.sample = data.table::data.table(geno[,names(geno)[c(10:ncol(geno))]])
        if(length(grep("\\/",geno.sample[1,]))>0){
          geno.sample[geno.sample=="./."]<-"NN"
          geno.sample[geno.sample=="0/0"]<-"00"
          geno.sample[geno.sample=="1/0"]<-"10"
          geno.sample[geno.sample=="0/1"]<-"01"
          geno.sample[geno.sample=="1/1"]<-"11"
        }else if(length(grep("\\|",geno.sample[1,]))>0){
          geno.sample[geno.sample=="0|0"]<-"00"
          geno.sample[geno.sample=="1|0"]<-"10"
          geno.sample[geno.sample=="0|1"]<-"01"
          geno.sample[geno.sample=="1|1"]<-"11"
        }
        geno.sample = cbind(geno[,c(1:9)],geno.sample)
        rm(geno)
        ## i=1
        for(i in 1:nrow(geno.sample)){
          geno.sample[i,10:ncol(geno.sample)] <- gsub("0",geno.sample$REF[i],as.character(geno.sample[i,10:ncol(geno.sample)]),)
          geno.sample[i,10:ncol(geno.sample)] <- gsub("1",geno.sample$ALT[i],as.character(geno.sample[i,10:ncol(geno.sample)]),)
        }
        hmp$`rs#` = paste0(geno.sample$`#CHROM`, "_", geno.sample$POS)
        hmp$alleles = paste0(geno.sample$REF,"/",geno.sample$ALT); hmp$chrom =geno.sample$`#CHROM`; hmp$pos = geno.sample$POS; hmp$strand = "+"
        hmp$`assembly#` = NA; hmp$center = NA; hmp$protLSID = NA; hmp$assayLSID = NA; hmp$panel = NA; hmp$QCcode = NA
        hmp <- cbind(hmp,geno.sample[,c(10:ncol(geno.sample))])
      }else if(is.geno(data=vcf)){
        names(geno)[c(1:5)] = c("ANNO","#CHROM","POS","REF","ALT")
        hmp$`rs#` = paste0(geno$`#CHROM`, "_", geno$POS)
        hmp$alleles = paste0(geno$REF,"/",geno$ALT); hmp$chrom =geno$`#CHROM`; hmp$pos = geno$POS; hmp$strand = "+"
        hmp$`assembly#` = NA; hmp$center = NA; hmp$protLSID = NA; hmp$assayLSID = NA; hmp$panel = NA; hmp$QCcode = NA
        hmp <- cbind(hmp,geno[,c(6:ncol(geno))])              
      }
      
      data.table::fwrite(hmp, file = paste0(dirPath, "/geno2hmp/", genoName, ".hmp"), row.names = F, col.names = T, quote = F, sep = "\t")
    }
  }
  
  if(nthreads == 1){
    geno2hmp(genoFile)
  }
  if(nthreads > 1){
    parLapply(cl.core, genoFile, geno2hmp)
    stopCluster(cl.core)
  }
  
}
