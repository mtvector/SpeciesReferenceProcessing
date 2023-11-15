#' This modification of the original ArchR function comes from Nelson Johansen
#' Create a gene annotation object for ArchR
#' 
#' This function will create a gene annotation object that can be used for creating ArrowFiles or an ArchRProject, etc.
#' 
#' @param genome A string that specifies the genome (ie "hg38", "hg19", "mm10", "mm9"). If `genome` is not supplied,
#' `TxDb` and `OrgDb` are required. If genome is supplied, `TxDb` and `OrgDb` will be ignored.
#' @param TxDb A `TxDb` object (transcript database) from Bioconductor which contains information for gene/transcript coordinates.
#' For example, from `txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene`.
#' @param mapping Dataframe that maps between gene_id and gene_symbol typically from biomaRt as show in the example.
#' @param OrgDb An `OrgDb` object (organism database) from Bioconductor which contains information for gene/transcript symbols from ids.
#' For example, from `orgdb <- org.Hs.eg.db`.
#' @param genes A `GRanges` object containing gene coordinates (start to end). Must have a symbols column matching the symbols column of `exons`.
#' @param exons A `GRanges` object containing gene exon coordinates. Must have a symbols column matching the symbols column of `genes`.
#' @param TSS A `GRanges` object containing standed transcription start site coordinates for computing TSS enrichment scores downstream.
#' @param annoStyle annotation style to map between gene names and various gene identifiers e.g. "ENTREZID", "ENSEMBL".
#'
#' @example 
#'
#' genomeAnnotation <- createGenomeAnnotation(genome = "BSgenome.Sscrofa.UCSC.susScr11")
#' txdb <- makeTxDbFromEnsembl(organism = "Sus Scrofa")
#'
#' ## Manual mapping between gene_id and symbol for "unsupported" species
#' library("biomaRt")                                                                                                                            
#' httr::set_config(httr::config(ssl_verifypeer = FALSE))                                                                                                                
#' ensembl <- useMart("ensembl", dataset="sscrofa_gene_ensembl") ##listDatasets(ensembl)                                                                     
#' mart.genes <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"), values=keys(txdb), mart=ensembl)   
#'
#' seqlevels(txdb) <- paste0(c(seq(1,18), "X", "Y"))
#' seqlevels(txdb) <- paste0("chr", seqlevels(txdb))
#'
#' ## Using a custom version of createGeneAnnotation that takes a data.frame to map gene_id to hgnc_symbol
#' geneAnnotation <- createGeneAnnotationUnsupportedSpecies(TxDb = txdb,
#'                                                          mapping = mart.genes,
#'                                                          OrgDb = org.Ss.eg.db)
#'
#' ## Remove genes without symbol
#' loci <- grep("NA", geneAnnotation$genes$symbol)
#' gid <- geneAnnotation$genes$gene_id[-loci]
#' df <- select(txdb, keys = gid, columns="TXNAME", keytype="GENEID")
#'
#' genes <- geneAnnotation$genes[-loci]
#' exons <- geneAnnotation$exons[-grep("NA", geneAnnotation$exons$symbol)]
#' tss <- geneAnnotation$TSS[which(geneAnnotation$TSS$tx_name %in% df$TXNAME)]
#'
#' ## Finalize annotation
#' geneAnnotationSubset <- createGeneAnnotation(genes = genes, 
#'                                              exons = exons, 
#'                                              TSS = tss)
#' save(genomeAnnotation, geneAnnotation, geneAnnotationSubset, file = "susScr11.RData")
#'
#' @export
createGeneAnnotationUnsupportedSpecies <- function(
  genome = NULL,
  TxDb = NULL,
  mapping = NULL,
  OrgDb = NULL,
  genes = NULL,
  exons = NULL,
  TSS = NULL,
  annoStyle = NULL
  ){

  # .validInput(input = genome, name = "genome", valid = c("character", "null"))
  # .validInput(input = TxDb, name = "TxDb", valid = c("txdb", "character", "null"))
  # .validInput(input = OrgDb, name = "OrgDb", valid = c("orgdb", "character", "null"))
  # .validInput(input = genes, name = "genes", valid = c("GRanges", "null"))
  # .validInput(input = exons, name = "exons", valid = c("GRanges", "null"))
  # .validInput(input = TSS, name = "TSS", valid = c("GRanges", "null"))
  # .validInput(input = annoStyle, name = "annoStyle", valid = c("character", "null"))

  if(is.null(genes) | is.null(exons) | is.null(TSS)){

    inGenes <- genes
    inExons <- exons
    inTSS <- TSS

    # .requirePackage("GenomicFeatures", source = "bioc")
    library(GenomicFeatures)

    # if(is.null(genome)) {
    #   if (is.null(TxDb) | is.null(OrgDb)) {
    #       stop("If no provided genome then you need TxDb and OrgDb!")
    #   }
    # }

    #if(!is.null(genome)){
    #  TxDb <- .getTxDb(genome)
    #  OrgDb <- .getOrgDb(genome)
    #}

    ###########################
    message("Getting Genes..")
    genes <- GenomicFeatures::genes(TxDb)

    if(is.null(annoStyle)){
      isEntrez <- mcols(genes)$symbol <- tryCatch({
        suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
        TRUE
      }, error = function(x){
        FALSE
      })

      isEnsembl <- mcols(genes)$symbol <- tryCatch({
        suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
        TRUE
      }, error = function(x){
        FALSE
      })

      if(isEntrez){
        annoStyle <- "ENTREZID"
      }else if(isEnsembl){
        annoStyle <- "ENSEMBL"
      }else{
        if(!is.null(mapping)){
          annoStyle <- "MANUAL"
        }else{
          stop("Could not identify keytype for annotation format!")
        }
      }
    }
    annoStyle <- toupper(annoStyle)

    message("Determined Annotation Style = ", annoStyle)

    ###########################
    if(annoStyle == "MANUAL"){
      mcols(genes)$symbol = mapping$hgnc_symbol[match(mcols(genes)$gene_id, mapping$ensembl_gene_id)]
      mcols(genes)$symbol[mcols(genes)$symbol == ""] = NA
      mcols(genes)$symbol[is.na(mcols(genes)$symbol)] <- paste0("NA_", mcols(genes)$gene_id)[is.na(mcols(genes)$symbol)]
    }else{
      mcols(genes)$symbol <- suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = annoStyle, multiVals = "first"))
      mcols(genes)$symbol[is.na(mcols(genes)$symbol)] <- paste0("NA_", mcols(genes)$gene_id)[is.na(mcols(genes)$symbol)]
    }
    names(genes) <- NULL
    genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)

    ###########################
    message("Getting Exons..")
    exons <- unlist(GenomicFeatures::exonsBy(TxDb, by = "tx"))
    exons$tx_id <- names(exons)
    mcols(exons)$gene_id <- suppressMessages(AnnotationDbi::select(TxDb, keys = paste0(mcols(exons)$tx_id), column = "GENEID", keytype = "TXID")[, "GENEID"])
    exons <- exons[!is.na(mcols(exons)$gene_id), ]
    if(annoStyle == "MANUAL"){
      mcols(exons)$symbol = mapping$hgnc_symbol[match(mcols(exons)$gene_id, mapping$ensembl_gene_id)]
      mcols(exons)$symbol[mcols(exons)$symbol == ""] = NA
      mcols(exons)$symbol[is.na(mcols(exons)$symbol)] <- paste0("NA_", mcols(exons)$gene_id)[is.na(mcols(exons)$symbol)]
    }else{
      mcols(exons)$symbol <- suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(exons)$gene_id, column = "SYMBOL", keytype = annoStyle, multiVals = "first"))
      mcols(exons)$symbol[is.na(mcols(exons)$symbol)] <- paste0("NA_", mcols(exons)$gene_id)[is.na(mcols(exons)$symbol)]
    }
    names(exons) <- NULL
    mcols(exons)$exon_id <- NULL
    mcols(exons)$exon_name <- NULL
    mcols(exons)$exon_rank <- NULL
    mcols(exons)$tx_id <- NULL
    exons <- sort(sortSeqlevels(exons), ignore.strand = TRUE)

    ###########################
    message("Getting TSS..")
    TSS <- unique(resize(GenomicFeatures::transcripts(TxDb), width = 1, fix = "start"))

    # if(!is.null(inGenes)){
    #   genes <- .validGRanges(inGenes)
    # }

    # if(!is.null(inExons)){
    #   exons <- .validGRanges(inExons)
    # }

    # if(!is.null(inTSS)){
    #   TSS <- .validGRanges(inTSS)
    # }

  }else{

    genes <- .validGRanges(genes)
    exons <- .validGRanges(exons)
    TSS <- unique(.validGRanges(TSS))

  }

  SimpleList(genes = genes, exons = exons, TSS = TSS)

}

.getTxDb <- function(genome = NULL, filter = TRUE, install = TRUE){
  
  .validInput(input = genome, name = "genome", valid = "character")
  .validInput(input = filter, name = "filter", valid = "boolean")
  .validInput(input = install, name = "install", valid = "boolean")
  
  if(toupper(genome) == "HG19"){
    if(suppressWarnings(!require(TxDb.Hsapiens.UCSC.hg19.knownGene))){
      if(install){
        message("Package does not exist, now trying bioconductor..")
        BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", update=FALSE)
      }else{
        stop("TxDb.Hsapiens.UCSC.hg19.knownGene is not installed!")
      }
    }
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  }else if(toupper(genome) == "HG38"){
    if(suppressWarnings(!require(TxDb.Hsapiens.UCSC.hg38.knownGene))){
      if(install){
        message("Package does not exist, now trying bioconductor..")
        BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", update=FALSE)
      }else{
        stop("TxDb.Hsapiens.UCSC.hg38.knownGene is not installed!")
      }
    }
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  }else if(toupper(genome) == "MM9"){
    if(suppressWarnings(!require(TxDb.Mmusculus.UCSC.mm9.knownGene))){
      if(install){
        message("Package does not exist, now trying bioconductor..")
        BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene", update=FALSE)
      }else{
        stop("TxDb.Mmusculus.UCSC.mm9.knownGene is not installed!")
      }
    }
    library(TxDb.Mmusculus.UCSC.mm9.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
  }else if(toupper(genome) == "MM10"){
    if(suppressWarnings(!require(TxDb.Mmusculus.UCSC.mm10.knownGene))){
      if(install){
        message("Package does not exist, now trying bioconductor..")
        BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", update=FALSE)
      }else{
        stop("TxDb.Mmusculus.UCSC.mm10.knownGene is not installed!")
      }
    }
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  }else if(toupper(genome) == "SACCER3"){
    if(suppressWarnings(!require(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene))){
      if(install){
        message("Package does not exist, now trying bioconductor..")
        BiocManager::install("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", update=FALSE)
      }else{
        stop("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene is not installed!")
      }
    }
    library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
    txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
  }else if(toupper(genome) == "RHEMAC8"){
    if(suppressWarnings(!require(TxDb.Mmulatta.UCSC.rheMac8.refGene))){
      if(install){
        message("Package does not exist, now trying bioconductor..")
        BiocManager::install("TxDb.Mmulatta.UCSC.rheMac8.refGene", update=FALSE)
      }else{
        stop("TxDb.Mmulatta.UCSC.rheMac8.refGene is not installed!")
      }
    }
    library(TxDb.Mmulatta.UCSC.rheMac8.refGene)
    txdb <- TxDb.Mmulatta.UCSC.rheMac8.refGene
  }else{
    stop("Genome not recognized!")
  }
  
  if(filter){
    txdb <- filterChrGR(txdb)
  }
  
  return(txdb)
  
}

.validInput <- function(input = NULL, name = NULL, valid = NULL){
  
  valid <- unique(valid)
  
  if(is.character(valid)){
    valid <- tolower(valid)
  }else{
    stop("Validator must be a character!")
  }
  
  if(!is.character(name)){
    stop("name must be a character!")
  }
  
  if("null" %in% tolower(valid)){
    valid <- c("null", valid[which(tolower(valid) != "null")])
  }
  
  av <- FALSE
  
  for(i in seq_along(valid)){
    
    vi <- valid[i]
    
    if(vi == "integer" | vi == "wholenumber"){
      
      if(all(is.numeric(input))){
        #https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
        cv <- min(abs(c(input%%1, input%%1-1)), na.rm = TRUE) < .Machine$double.eps^0.5
      }else{
        cv <- FALSE
      }
      
    }else if(vi == "null"){
      
      cv <- is.null(input)
      
    }else if(vi == "bool" | vi == "boolean" | vi == "logical"){
      
      cv <- is.logical(input)
      
    }else if(vi == "numeric"){
      
      cv <- is.numeric(input)
      
    }else if(vi == "vector"){
      
      cv <- is.vector(input)
      
    }else if(vi == "matrix"){
      
      cv <- is.matrix(input)
      
    }else if(vi == "sparsematrix"){
      
      cv <- is(input, "dgCMatrix")
      
    }else if(vi == "character"){
      
      cv <- is.character(input)
      
    }else if(vi == "factor"){
      
      cv <- is.factor(input)
      
    }else if(vi == "rlecharacter"){
      
      cv1 <- is(input, "Rle")
      if(cv1){
        cv <- is(input@values, "factor") || is(input@values, "character")
      }else{
        cv <- FALSE
      }
      
    }else if(vi == "palette"){
      
      cv <- all(.isColor(input))
      
    }else if(vi == "timestamp"){
      
      cv <- is(input, "POSIXct")
      
    }else if(vi == "dataframe" | vi == "data.frame" | vi == "df"){
      
      cv1 <- is.data.frame(input)
      cv2 <- is(input, "DataFrame")
      cv <- any(cv1, cv2)
      
    }else if(vi == "fileexists"){
      
      cv <- all(file.exists(input))
      
    }else if(vi == "direxists"){
      
      cv <- all(dir.exists(input))
      
    }else if(vi == "granges" | vi == "gr"){
      
      cv <- is(input, "GRanges")
      
    }else if(vi == "grangeslist" | vi == "grlist"){
      
      cv <- .isGRList(input)
      
    }else if(vi == "list" | vi == "simplelist"){
      
      cv1 <- is.list(input)
      cv2 <- is(input, "SimpleList")
      cv <- any(cv1, cv2)
      
    }else if(vi == "bsgenome"){
      
      cv1 <- is(input, "BSgenome")
      cv2 <- tryCatch({
        library(input)
        eval(parse(text=input))
      }, error = function(e){
        FALSE
      })
      cv <- any(cv1, cv2)
      
    }else if(vi == "se" | vi == "summarizedexperiment"){
      
      cv <- is(input, "SummarizedExperiment")
      
    }else if(vi == "seurat" | vi == "seuratobject"){
      
      cv <- is(input, "Seurat")
      
    }else if(vi == "txdb"){
      
      cv <- is(input, "TxDb")
      
    }else if(vi == "orgdb"){
      
      cv <- is(input, "OrgDb")
      
    }else if(vi == "bsgenome"){
      
      cv <- is(input, "BSgenome")
      
    }else if(vi == "parallelparam"){
      
      cv <- is(input, "BatchtoolsParam")
      
    }else if(vi == "archrproj" | vi == "archrproject"){
      
      cv <- is(input, "ArchRProject")
      ###validObject(input) check this doesnt break anything if we
      ###add it. Useful to make sure all ArrowFiles exist! QQQ
      
    }else{
      
      stop("Validator is not currently supported by ArchR!")
      
    }
    
    if(cv){
      av <- TRUE
      break
    }   
    
  }
  
  if(av){
    
    return(invisible(TRUE))
    
  }else{
    
    stop("Input value for '", name,"' is not a ", paste(valid, collapse="," ), ", (",name," = ",class(input),") please supply valid input!")
    
  }
  
}

#https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
.isWholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}

#https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
.isColor <- function(x = NULL){
  unlist(lapply(x, function(y) tryCatch(is.matrix(col2rgb(y)), error = function(e) FALSE)))
}

.isDiscrete <- function(x = NULL){
  is.factor(x) || is.character(x) || is.logical(x)
}

.isGRList <- function(x){
  isList <- grepl("list", class(x), ignore.case=TRUE)
  if(!isList){
    FALSE
  }else{
    allGR <- all(unlist(lapply(x, function(x) is(x, "GRanges") )))
    if(allGR){
      TRUE
    }else{
      FALSE
    }
  }
}

#' Get/Validate BSgenome
#' 
#' This function will attempt to get or validate an input as a BSgenome.
#' 
#' @param genome This option must be one of the following: (i) the name of a valid ArchR-supported genome ("hg38", "hg19", or "mm10"),
#' (ii) the name of a `BSgenome` package (for ex. "BSgenome.Hsapiens.UCSC.hg19"), or (iii) a `BSgenome` object.
#' @param masked A boolean describing whether or not to access the masked version of the selected genome. See `BSgenome::getBSgenome()`.
#' @export
validBSgenome <- function(genome = NULL, masked = FALSE){
  
  .validInput(input = genome, name = "genome", valid = c("character", "bsgenome"))
  .validInput(input = masked, name = "masked", valid = c("boolean"))
  
  stopifnot(!is.null(genome))
  if(inherits(genome, "BSgenome")){
    return(genome)
  }else if(is.character(genome)){
    genome <- tryCatch({
      .requirePackage(genome)
      bsg <- eval(parse(text = genome))
      if(inherits(bsg, "BSgenome")){
        return(bsg)
      }else{
        stop("genome is not a BSgenome valid class!")
      }
    }, error = function(x){
      BSgenome::getBSgenome(genome, masked = masked)
    })  
    return(genome)
  }else{
    stop("Cannot validate BSgenome options are a valid BSgenome or character for getBSgenome")
  }  
}

.validTxDb <- function(TxDb = NULL){
  stopifnot(!is.null(TxDb))
  if(inherits(TxDb, "TxDb")){
    return(TxDb)
  }else if(is.character(TxDb)){
    return(getTxDb(TxDb)) #change
  }else{
    stop("Cannot validate TxDb options are a valid TxDb or character for getTxDb")
  }
}

.validOrgDb <- function(OrgDb = NULL){
  stopifnot(!is.null(OrgDb))
  if(inherits(OrgDb, "OrgDb")){
    return(OrgDb)
  }else if(is.character(OrgDb)){
    return(getOrgDb(OrgDb)) #change
  }else{
    stop("Cannot validate OrgDb options are a valid OrgDb or character for getOrgDb")
  }
}

.validGRanges <- function(gr = NULL){
  stopifnot(!is.null(gr))
  if(inherits(gr, "GRanges")){
    return(gr)
  }else{
    stop("Error cannot validate genomic range!")
  }
}

.validGeneAnnotation <- function(geneAnnotation = NULL){
  
  if(!inherits(geneAnnotation, "SimpleList")){
    if(inherits(geneAnnotation, "list")){
      geneAnnotation <- as(geneAnnotation, "SimpleList")
    }else{
      stop("geneAnnotation must be a list/SimpleList of 3 GRanges for : Genes GRanges, Exons GRanges and TSS GRanges!")
    }
  }
  if(identical(sort(tolower(names(geneAnnotation))), c("exons", "genes", "tss"))){
    
    gA <- SimpleList()
    gA$genes <- .validGRanges(geneAnnotation[[grep("genes", names(geneAnnotation), ignore.case = TRUE)]])
    gA$exons <- .validGRanges(geneAnnotation[[grep("exons", names(geneAnnotation), ignore.case = TRUE)]])
    gA$TSS <- .validGRanges(geneAnnotation[[grep("TSS", names(geneAnnotation), ignore.case = TRUE)]])
    
  }else{
    stop("geneAnnotation must be a list/SimpleList of 3 GRanges for : Genes GRanges, Exons GRanges and TSS GRanges!")
  }
  
  gA
  
}

.validGenomeAnnotation <- function(genomeAnnotation = NULL){
  
  if(!inherits(genomeAnnotation, "SimpleList")){
    if(inherits(genomeAnnotation, "list")){
      genomeAnnotation <- as(genomeAnnotation, "SimpleList")
    }else{
      stop("genomeAnnotation must be a list/SimpleList of 3 GRanges for : blacklist GRanges, chromSizes GRanges and genome BSgenome package string (ie hg38 or BSgenome.Hsapiens.UCSC.hg38)!")
    }
  }
  
  if(identical(sort(tolower(names(genomeAnnotation))), c("blacklist", "chromsizes", "genome"))){
    
    gA <- SimpleList()
    gA$blacklist <- .validGRanges(genomeAnnotation[[grep("blacklist", names(genomeAnnotation), ignore.case = TRUE)]])
    if(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]]=="nullGenome"){
      gA$genome <- "nullGenome"
    }else{
      bsg <- validBSgenome(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]])
      gA$genome <- bsg@pkgname
    }
    gA$chromSizes <- .validGRanges(genomeAnnotation[[grep("chromsizes", names(genomeAnnotation), ignore.case = TRUE)]])
    
  }else{
    
    stop("genomeAnnotation must be a list/SimpleList of 3 GRanges for : blacklist GRanges, chromSizes GRanges and genome BSgenome package string (ie hg38 or BSgenome.Hsapiens.UCSC.hg38)!")
    
  }
  
  gA
  
}

.validGeneAnnoByGenomeAnno <- function(geneAnnotation, genomeAnnotation){
  
  allSeqs <- unique(paste0(seqnames(genomeAnnotation$chromSizes)))
  
  geneSeqs <- unique(paste0(seqnames(geneAnnotation$genes)))
  if(!all(geneSeqs %in% allSeqs)){
    geneNotIn <- geneSeqs[which(geneSeqs %ni% allSeqs)]
    message("Found Gene Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(geneNotIn, collapse=","))
    geneAnnotation$genes <- .subsetSeqnamesGR(geneAnnotation$genes, names = allSeqs)
  }
  
  exonSeqs <- unique(paste0(seqnames(geneAnnotation$exons)))
  if(!all(exonSeqs %in% allSeqs)){
    exonNotIn <- exonSeqs[which(exonSeqs %ni% allSeqs)]
    message("Found Exon Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(exonNotIn, collapse=","))
    geneAnnotation$exons <- .subsetSeqnamesGR(geneAnnotation$exons, names = allSeqs)
  }
  
  TSSSeqs <- unique(paste0(seqnames(geneAnnotation$TSS)))
  if(!all(TSSSeqs %in% allSeqs)){
    TSSNotIn <- TSSSeqs[which(TSSSeqs %ni% allSeqs)]
    message("Found TSS Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(TSSNotIn, collapse=","))
    geneAnnotation$TSS <- .subsetSeqnamesGR(geneAnnotation$TSS, names = allSeqs)
  }
  
  geneAnnotation
  
}


.validArchRProject <- function(ArchRProj = NULL){
  if(!inherits(ArchRProj, "ArchRProject")){
    stop("Not a valid ArchRProject as input!")
  }else{
    ArchRProj
  }
}