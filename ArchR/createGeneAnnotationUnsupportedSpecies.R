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

    if(!is.null(genome)){
      TxDb <- .getTxDb(genome)
      OrgDb <- .getOrgDb(genome)
    }

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
