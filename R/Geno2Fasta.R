#' This is some description of this function.
#' @title to transform geno data into fasta sequence
#'
#' @description By using this package, you could use the function of Geno2Fasta to transform geno data into fasta sequence.
#'
#' @details see above
#'
#' @param genositeFile: geno site classification folder;
#' @param nthreads: to define how many cores of this computer you would like to run the task. Make sure not to excede the maximum of this computer. Default value is 1.
#' 
#' @return files extracted and stored in the folder created by default.
#' @export Geno2Fasta
#' @examples Geno2Fasta(genositeFile = "./Geno.Site.Classification/", nthreads = 4)


Geno2Fasta <- function(genoSiteFile = NULL, nthreads = 1, ...){
  library(data.table)
  library(parallel)
  library(Biostrings)

  dirPath <- getwd()
  if(!dir.exists(paste0(dirPath, "/Geno2Fasta"))){
    dir.create(paste0(dirPath, "/Geno2Fasta"))
  }

  if(dir.exists(genoSiteFile)){
    genoList = paste0(genoSiteFile, list.files(genoSiteFile))
  }else{
    genoList = genoSiteFile
  }

  if(nthreads > 1){
    cores = detectCores()
    if(nthreads <= cores){
      cl.core = makeCluster(getOption("cl.cores", as.numeric(nthreads)))
    }
    if(nthreads > cores) stop(paste0("'nthreads' to be used is out of the maximum core number of this computer: ", cores, ". Please refine the 'nthreads' parameter."))
  }
  alleleType = "double"
  
  ## function
  vector2string = function(data = NULL,sep = "", ...){
    if(is.null(data)) stop("'data' is required for vector2string!")
    df = data[1]
    if(length(data) >= 2){
      for(i in 2:length(data)){
        df = paste(df, data[i], sep = sep)
      }
    }
    return(df)
  }

  ## function
  diploid2haploid = function(str){
    str = gsub("AG|GA", "R", str)
    str = gsub("CT|TC", "Y", str)
    str = gsub("AC|CA", "M", str)
    str = gsub("GT|TG", "K", str)
    str = gsub("GC|CG", "S", str)
    str = gsub("AT|TA", "W", str)
    str = gsub("AA", "A", str)
    str = gsub("CC", "C", str)
    str = gsub("TT", "T", str)
    str = gsub("GG", "G", str)
    str = gsub("--", "N", str)
    str = gsub("NN", "N", str)
    return(str)
  }

  # function
  matrix2fasta = function(mat){
    mat = list(mat)
    str = mat[[1]][1]
    for(i in 2:length(mat[[1]])){
      str = paste0(str, mat[[1]][i])
    }
    return(str)
  }

  # a = 15; i = 18
  geno2fas = function(genoList){
    for(a in 1:length(genoList)){
      geno = data.table::fread(genoList[a], header = T)
      data.table::setkeyv(geno, names(geno))
      col_name = names(geno)
      fileName = gsub(".txt|.csv|.geno|.hmp", "", strsplit(genoList[a], "/")[[1]][length(strsplit(genoList[a], "/")[[1]])])
      from = grep("QCcode", names(geno)) + 1
      to = ncol(geno)
      if(nrow(geno) >= 1){
        if(alleleType == "double"){
          geno = data.frame(geno); names(geno) = col_name
          genos = geno[,c(grep("REF", names(geno)),from:to)]
          genos = apply(genos, c(1,2), diploid2haploid)
          geno = cbind(geno[,c(1:(from-1))], genos)
          data.table::fwrite(geno, file = paste0(dirPath, "/", "Geno2Fasta/", fileName, ".single.allele", ".csv"), row.names = F, col.names = T, quote = F, sep = "\t")
          fas = Biostrings::DNAStringSet(apply(genos, 2, matrix2fasta))
          Biostrings::writeXStringSet(fas, file = paste0(dirPath, "/Geno2Fasta/", fileName, ".fasta"))
        }
      }
    }
  }

  ## function
  if(nthreads == 1){
    geno2fas(genoList)
  }
  if(nthreads > 1){
    parLapply(cl.core, genoList, geno2fas)
    stopCluster(cl.core)
  }
}
