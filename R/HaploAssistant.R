#' This is some descriptions of this function.
#' @title haplotype analysis for specific genes.
#'
#' @description By using this package, you could use functions of HaploAssistant to haplotype analysis for specific genes.
#'
#' @details see above
#'
#' @param categoryTransforming: determine whether to convert latitude and longitude data format;
#' @param GenoTransforming: determine whether to convert the file format from .vcf/.geno to .hapmap. Default is 'F'.
#' @param GenoExtract: determine whether to extract sample data from original Hapmap file. Default is "T".
#' @param GeoPlotting: determine whether to conduct geographic mapping. Default is "T".
#' @param trait: the trait you are working with;
#' @param SNPAnnoFile: SNP annotation files for users self sequencing;
#' @param categoryFile: file stores category of each sample;
#' @param gffFile: generic feature format file of certain genome. This file defines the start and end position of specific gene, it should contains at least "seqid", "source", "type", "start", "end", "score", "strand", and "attributes";
#' @param genomeFasta: genomic sequence in fasta format, which is required;
#' @param phenoFile: trait value of each sample;
#' @param geneModel: gene model you want to perform haplotype analysis, or a file with a column named "Gene_ID";
#' @param genoFile: a folder stores the hapmap file of each chromosome, or a single hapmap file of specific chromosome;
#' @param mapData: geographic data in 'shp' or 'json' format. If 'mapData' is NULL, map data would be downloaded from rnaturalearth (http://www.naturalearth.com).
#' @param haploMin: the minimum samples of a specific haplotype taking into account;
#' @param pdfWidth: numeric value to indicate the width of pdf file, default is 15;
#' @param pdfHeight: numeric value to indicate the height of pdf file, default is 10;
#' @param plotMargin: adjust graph edge margins;
#' @param nthreads: to define how many cores of this computer you would like to run the task. Make sure not to excede the maximum of this computer. Default value is 1;
#' @param geneScope: a region to define the upstream and downstream of a gene, genomic position out of this scope would be defined as intergenic. Default value is 2000;
#' @param size: height of gene body, default value is 5;
#' @param lcolor: color of the middle line, default value is "red";
#' @param type: could be either coastline or countries, default is countries;
#' @param region: the region to display on the map, could be "world", "China", or other countries, default is 'China';
#' @param mapResolution: value of 10, 50 and 110, or string of 'small', 'medium', or 'large', default is 50;
#' @param mapBorderSize: border size of geographic mapping, default is 0.1.
#' @param longLim: limitation of longitude, default is c(-180.00, 180.00);
#' @param latLim: limitation of latitude, default is c(-90.00, 90.00);
#' @param landColor: the color of land, default is lightyellow;
#' @param oceanColor: the color of ocean, default is white;
#' @param pieMax: numeric value to indicate the maximum size of pie plot on geographic map, default is 2;
#' @param expand: make the axes not overlap with each other, default is "T";
#' @param pieLinetype: line type inherited from ggplot2, default is 2;
#' @param pieLinesize: line size inherited from ggplot2, default is 0;
#' @param pieBorderColor: border color of pie plot, default is NA;
#' @param pieAlpha: numeric value to define transparency of pie plot, default is 0.8;
#' @param label_graticule: parameter of coord_sf(), which is a character vector indicating which graticule lines should be labeled where. It could be 'N', 'S', 'W', 'E', or combination of these characters. Default is 'EW'.
#' @param MappingDataType: data types of geographical mapping, default is "haplotype";
#' @param dataFiltration: logic value to indicatee whether to filtrate data by 'China', or other keywords, default is NULL.
#'
#' @return Files extracted and stored in the folder created by default. File types of hapmap, log, csv, fasta, haplotype, pdf.
#' @export HaploAssistant
#' @examples HaploAssistant(categoryTransforming = F, GenoTransforming = F, GenoExtract = T, GeoPlotting = T,
#'                          trait = "Trait", categoryFile = "./Trait.category.csv",
#'                          gffFile = "./Wm82.a2.v1.gene.gff", genomeFasta = "./Gmax_275_v2.0.fa",
#'                          phenoFile ="./Trait.csv", geneModel="Glyma.18G056600",
#'                          genoFile = "./Hapmap_Gm18.txt", mapData = "./China.json")


HaploAssistant = function(categoryTransforming=F,GenoTransforming = F,GenoExtract=T,GeoPlotting=T,
                          trait = NULL,SNPAnnoFile=NULL,categoryFile = NULL,gffFile = NULL,genomeFasta = NULL,phenoFile = NULL,geneModel=NULL,
                          genoFile = NULL,mapData = NULL,haploMin = 10, pdfWidth = 10, pdfHeight = NULL, plotMargin = NULL, nthreads = 1,
                          geneScope = 2000,size = 5, lcolor= "red",type = "countries",region = "China",mapResolution = 50,mapBorderSize = 0.1,
                          longLim = NULL,latLim = NULL,landColor = "lightyellow",oceanColor = "white",pieMax = 2,expand = T,dataFiltration = NULL,
                          pieLinetype = 2,pieLinesize = 0, pieBorderColor = NA,pieAlpha = 0.8,label_graticule = "EW",MappingDataType = "haplotype"){
  dirPath = getwd()
  if(GenoTransforming){
    Geno2Hmp(genoFile=genoFile,nthreads = nthreads)
    genoFile = paste0(dirPath,"/geno2hmp/")
  }
  if(categoryTransforming){
    LongLatTransforming(categoryFile=categoryFile)
    file.list = list.files(paste0(dirPath, "/LongLatTransformation/"))
    categoryFile = paste0(dirPath, "/LongLatTransformation/", file.list)[1]
  }
  if(GenoExtract){
    GenoExtractFromHapmap(trait = trait, categoryFile = categoryFile, genoFile = genoFile)
    sampleHapmapFile = paste0(dirPath,"/", trait, "_geno.extracted.from.hapmap/")
  }else sampleHapmapFile = genoFile
  RegionalSNPextractionFromHapmap(sampleHapmapFile = sampleHapmapFile, gffFile= gffFile, geneModel = geneModel)
  hapmapFile = paste0(dirPath,"/Hapmap.Regional.Extracted.Files/")
  if(!is.null(SNPAnnoFile)){
    GeneSpecificSNPannotation(SNPAnnoFile = SNPAnnoFile, geneModel =  geneModel, gffFile = gffFile)
  }else SNPAnnotation (genomeFasta = genomeFasta,gffFile = gffFile, hapmapFile = hapmapFile, geneScope = 2000, nthreads = 1)
  GenoSiteClassification(annoFile =  paste0(dirPath,"/SNP.Annotation/"), hapmapFile = hapmapFile)
  Geno2Fasta(genoSiteFile = paste0(dirPath,"/Geno.Site.Classification/"), nthreads = 1)
  HaplotypeIdentification (fastaFile = paste0(dirPath,"/Geno2Fasta/"), categoryFile = categoryFile, phenoFile = phenoFile, plotMargin = 2, pdfWidth = 5, pdfHeight = 5, nthreads = 1)
  haploResults = paste0(dirPath,"/Haplotype.Analysis/")
  if(!is.null(phenoFile)) HaplotypePlotting (haploResults = haploResults, haploMin = 10, pdfWidth = 10, pdfHeight = 10, plotMargin = 2, nthreads = 1)
  DrawSingleGeneStructure(gffFile = gffFile , haploResults = haploResults, size = 5, lcolor= "red", pdfWidth = 16, pdfHeight = 4)
  NexForPopART(haploResults = haploResults, categoryFile = categoryFile)
  if(GeoPlotting){
    GeoMapping(haploResults = haploResults ,mapData = mapData,type = "countries",
                region = "China",mapResolution = 50,longLim = NULL,latLim = NULL,
                landColor = "lightyellow",oceanColor = "white",pieMax = 2,pdfWidth = 15,
                pdfHeight = 10,expand = T,pieLinetype = 2,pieLinesize = 0,pieBorderColor = NA,
                pieAlpha = 0.8,label_graticule = "EW",plotMargin = 2,MappingDataType = "haplotype",mapBorderSize = 0.1,dataFiltration = NULL)
  }

  ### log file ###
  sink(paste0(dirPath, "/", trait,"_HaploAssistant.log"), type = "output", append = T)
  writeLines("###################################")
  writeLines("### Logs for HaploAssistant ###")
  writeLines("###################################\n")
  writeLines(paste0("Date: ", date()))
  writeLines("##########################")
  writeLines("### PARAMETERS SETTING ###")
  writeLines("##########################\n")
  writeLines(paste0("dirPath: ",dirPath))
  writeLines(paste0("trait: ", trait))
  writeLines(paste0("categoryFile: ", categoryFile))
  writeLines(paste0("gffFile: ", gffFile))
  writeLines(paste0("genomeFasta: ", genomeFasta))
  if(!is.null(phenoFile)){writeLines(paste0("phenoFile: ", phenoFile))
  }else writeLines("phenoFile: Simulation of standard normal distribution")
  writeLines(paste0("geneModel: ", geneModel))
  writeLines(paste0("genoFile: ", genoFile))
  writeLines(paste0("GenoTransforming: ", GenoTransforming))
  writeLines(paste0("GenoExtract: ", GenoExtract))
  writeLines(paste0("GeoPlotting: ", GeoPlotting))
  writeLines(paste0("categoryTransforming: ", categoryTransforming))
  writeLines(paste0("haploMin: ", haploMin))
  writeLines(paste0("size: ", size))
  if(GeoPlotting){
    writeLines(paste0("mapData: ", mapData))
    writeLines(paste0("type: ", type))
    writeLines(paste0("region: ", region))
    writeLines(paste0("mapResolution: ", mapResolution))
    writeLines(paste0("longLim: ", longLim))
    writeLines(paste0("latLim: ", latLim))
    writeLines(paste0("mapBorderSize: ", mapBorderSize))
    writeLines(paste0("landColor: ",landColor))
    writeLines(paste0("oceanColor: ", oceanColor))
    writeLines(paste0("pieMax: ", pieMax))
    writeLines(paste0("expand: ", expand))
    writeLines(paste0("pieLinetype: ", pieLinetype))
    writeLines(paste0("pieBorderColor: ", pieBorderColor))
    writeLines(paste0("pieLinesize: ", pieLinesize))
    writeLines(paste0("pieAlpha: ", pieAlpha))
    writeLines(paste0("label_graticule: ", label_graticule))
    writeLines(paste0("MappingDataType: ", MappingDataType))
    writeLines(paste0("dataFiltration : ", dataFiltration))
  }
  writeLines("\n# END")
  sink(type = "output")
}
