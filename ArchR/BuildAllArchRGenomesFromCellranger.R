#This script processes a CSV file containing genome annotation directory 
#information and generates ArchR reference files for ATAC-seq analysis.

#Use singularity /allen/programs/celltypes/workgroups/rnaseqanalysis/bicore/singularity/10x_multiome_qc_4.2.0.sif
#You need to have UCSC tools directory in path or give full path for faToTwoBit
#write("TMP = '/home/matthew.schmitz/Matthew/tmp/'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))
#.libPaths(new = '/home/matthew.schmitz/R')
.libPaths(new = '/home/matthew.schmitz/Matthew/utils/R')
#Set the working directory to where you want the genome and org.db package to go
intermediate_path <- '~/Matthew/archr_genomes/'
setwd(intermediate_path)
#CSV of the genomes spreadsheet
tab=read.csv2('~/Reference_Genome_tracking.csv',header = T,sep=',')
#tab=tab[file.exists(tab$CR6.ARC.2.0.reference.in.BICore.folder),]
#tab=tab[3,]

#English name
en_val <- 'English.Name'
#Scientific name
sp_val <- 'Species'
#Directory to cellranger ref top dir
dir_val <- 'CR6.ARC.2.0.reference.in.BICore.folder'
#ref version name
rn_val <- 'Top.NCBI.assembly'
#URL where reference was downloaded from
url_val <- 'RefSeq.anno.for.CR6'
#assembly date
ad_val <- 'Assembly.Date'
#species ID
spid_val <- 'NCBI.TaxID'
#provider
pr_val <- 'NCBI'
#subpath to the genome fasta file beneath dir_val
fa_subpath <- 'fasta/genome.fa'
#subpath to the gtf/gff file beneath dir_val (gff is untested)
gtf_subpath <- "genes/genes.gtf.gz"

if (!requireNamespace("BSgenome")){
  BiocManager::install("BSgenome",lib="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/R", force=T)
}
if (!requireNamespace("AnnotationForge")){
  BiocManager::install("AnnotationForge",lib="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/R", force=T)
}
library(BSgenome)
library(Biostrings)
library(rtracklayer)
library(AnnotationForge)
library(ArchR)
library(GenomicFeatures)
library(biomaRt)
library(httr)

source('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/code/SpeciesReferenceProcessing/ArchR/createGeneAnnotationUnsupportedSpecies.R')

smart_join <- function(...) {
  # Concatenate all arguments with a single slash
  joined_path <- paste(..., sep = "/")
  
  # Remove any double slashes
  cleaned_path <- gsub("//", "/", joined_path)
  
  return(cleaned_path)
}

parse_species_string <- function(input_string,sep="\\.") {
  parts <- unlist(strsplit(input_string, sep))
  initials <- substr(parts[-length(parts)], 1, 1)
  output <- paste0(tolower(paste(initials, collapse = "")), tolower(parts[length(parts)]))
  return(output)
}

create_seed_file <- function(english_name,species_name,assembly_date,
                             ref_name,ref_url,seqfile_name,
                             seqs_srcdir,ref_path,provider=pr_val) {
  
  package_name <- gsub(" ", ".", species_name)
  #Make the string to write the seed file
  seed_content <- paste0(
    "Package: ", package_name, "\n",
    "Title: ", species_name," (",english_name, ") genome\n",
    "Description: ", english_name, " genome built from 2bit file\n",
    "Version: 0.1.0\n",
    "organism: ", species_name, "\n",
    "common_name: ", english_name, "\n",
    "provider: ",provider,"\n",
    "genome: ", ref_name, "\n",
    "release_date: ", assembly_date, "\n",
    "source_url: ", ref_url, "\n",
    "seqfile_name: ", seqfile_name, "\n",
    "seqs_srcdir: ", seqs_srcdir, "\n",
    "organism_biocview: ", gsub(" ", ".", species_name), "\n",
    "BSgenomeObjname: ", gsub(" ", ".", species_name), "\n",
    "circ_seqs: character(0)"
  )
  
  # Write seed content to seed.txt file
  out_path=paste0(ref_path, gsub(" ", ".", species_name), "_seed.txt")
  write.table(seed_content, out_path,sep='\t',col.names=FALSE,quote=F,row.names = F)
  return(c(out_path,package_name))

}

for (rn in rownames(tab)) {
  cur=tab[rn,]
  archref_files <- c()#list.files(path = cur[[dir_val]], pattern = "_ArchRef\\.RData")
  
  # Check if any such files exist
  if (length(archref_files) > 0) {
    message(paste0("File(s) containing '_ArchRef.RData' found in ",cur[[dir_val]]))
  } else {
    
  print(cur[[sp_val]])
  print(cur[[dir_val]])
  fa_in=smart_join(cur[[dir_val]],fa_subpath)
  twob_out=gsub('\\.fa','.2bit',x=fa_in)
  command=paste("faToTwoBit",fa_in,twob_out)
  system(paste0("bash -c 'source ~/.bashrc;",command,"'"), intern = TRUE)
  output_elements=create_seed_file(english_name=cur[[en_val]],
                   species_name=cur[[sp_val]],
                   provider=pr_val,
                   assembly_date=cur[[ad_val]],
                   ref_name=cur[[rn_val]],
                   ref_url=cur[[url_val]],
                   seqfile_name=basename(twob_out),
                   seqs_srcdir=dirname(twob_out),
                   ref_path=cur[[dir_val]]
                   )
  print('making bsg')
  genome_package_path=smart_join('./',output_elements[2])
  
  if (dir.exists(genome_package_path) & output_elements[2]!="") {
    print(paste("already exists so recurively deleting genome_package_path", genome_package_path))
    unlink(genome_package_path, recursive = TRUE)
  }
  
  BSgenome::forgeBSgenomeDataPkg(output_elements[1])
  print('load bsg')
  print(output_elements[2])
  tools:::.build_packages(output_elements[2])
  print(genome_package_path)
  genome_package_env <- devtools::load_all(genome_package_path)
  genome_description_file <- smart_join(genome_package_path, 'DESCRIPTION')
  genome_package_info <- read.dcf(genome_description_file)
  bsg <- genome_package_env$env[[output_elements[2]]]
  grange <- GenomicRanges::GRanges(seqnames = bsg@seqinfo@seqnames,ranges = bsg@seqinfo@seqlengths,seqinfo=bsg@seqinfo)
  genomeAnnotation <- createGenomeAnnotation(genome = bsg,chromSizes = grange)
  
  print('gtf')
  txdb <- GenomicFeatures::makeTxDbFromGFF(file = smart_join(cur[[dir_val]],gtf_subpath), 
                          format="gtf",organism = cur[[sp_val]],taxonomyId = cur[[spid_val]])
  genes <- GenomicFeatures::genes(txdb)
  genus_name <- strsplit(cur[[sp_val]],split=" ")[[1]][1]
  species_name <- strsplit(cur[[sp_val]],split=" ")[[1]][2]
  print('make org')
  orgDBname <- paste0("org.",toupper(substr(genus_name,1,1)),species_name,".eg.db")
  org_path=smart_join('./',orgDBname,'/')
  
  if (dir.exists(org_path) & orgDBname!="") {
    print(paste("already exists, consider deleting org_path to rebuild package:", org_path))
    #print(unlink(org_path, recursive = TRUE))
  }else{
  #cache must be rebuilt each time unless you know which require
  #ARCHIVE NCBI version or not (requires unigene)
  #https://github.com/Bioconductor/AnnotationForge/issues/13
  AnnotationForge::makeOrgPackageFromNCBI(
    version = "0.1",
    author = "EvoGen <EvoGenes@evogen.edu>",
    maintainer = "EvoGen <EvoGenes@evogen.edu>",
    outputDir = ".",
    tax_id = cur[[spid_val]],
    genus = genus_name,
    species = species_name,
    rebuildCache = TRUE
  )

  #   tryCatch(
  #     {
  #       # First block of code to try
  #       AnnotationForge::makeOrgPackageFromNCBI(
  #         version = "0.1",
  #         author = "EvoGen <EvoGenes@evogen.edu>",
  #         maintainer = "EvoGen <EvoGenes@evogen.edu>",
  #         outputDir = ".",
  #         tax_id = cur[[spid_val]],
  #         genus = genus_name,
  #         species = species_name,
  #         rebuildCache = FALSE
  #       )
  #     },
  #     error = function(e) {
  #       # Code to run if an error occurs in the first block
  #       message("Error encountered: ", e$message, ". Trying with rebuildCache = TRUE.")
  #       AnnotationForge::makeOrgPackageFromNCBI(
  #         version = "0.1",
  #         author = "EvoGen <EvoGenes@evogen.edu>",
  #         maintainer = "EvoGen <EvoGenes@evogen.edu>",
  #         outputDir = ".",
  #         tax_id = cur[[spid_val]],
  #         genus = genus_name,
  #         species = species_name,
  #         rebuildCache = TRUE
  #       )
  #     }
  #   )
   }
  
  org_package_env <- devtools::load_all(org_path)
  org_description_file <- smart_join(org_path, 'DESCRIPTION')
  org_package_info <- read.dcf(org_description_file)
  org_db <- org_package_env$env[[orgDBname]]
  
  
  print(org_db)
  ## Using a custom version of createGeneAnnotation that takes a data.frame to map gene_id to hgnc_symbol
  httr::set_config(httr::config(ssl_verifypeer = FALSE))
  ensembl <- useMart("ensembl", dataset=paste0(parse_species_string(output_elements[2]),"_gene_ensembl")) ##listDatasets(ensembl)
  if('hgnc_symbol' %in% listAttributes(ensembl)$name){
    mart.genes <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"), values=keys(txdb), mart=ensembl)
  }
  else{
    mart.genes1 <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"), values=keys(txdb), mart=ensembl)
    mart.genes2 <- getBM(attributes=c("ensembl_gene_id","hsapiens_homolog_associated_gene_name"), values=keys(txdb), mart=ensembl)
    mart.genes <- merge(mart.genes1,mart.genes2)
    colnames(mart.genes) <- gsub("hsapiens_homolog_associated_gene_name","hgnc_symbol", colnames(mart.genes))
  }
  #Not sure when this would be helpful but leaving it just in case
  #seqlevels(txdb) <- paste0(c(seq(1,18), "X", "Y"))
  #seqlevels(txdb) <- paste0("chr", seqlevels(txdb))
  ## Using a custom version of createGeneAnnotation that takes a data.frame to map gene_id to hgnc_symbol
  geneAnnotation <- createGeneAnnotationUnsupportedSpecies(genome = NULL,TxDb = txdb,#genome_package_info[1,'Package']
                                                           OrgDb = org_db,
                                                           mapping=mart.genes,
                                                           annoStyle = "MANUAL")  
  ## Remove genes without symbol
  geneAnnotation$genes$symbol <- str_replace(geneAnnotation$genes$symbol,'^NA_','')
  geneAnnotation$exons$symbol <- str_replace(geneAnnotation$exons$symbol,'^NA_','')
  loci <- grepl("^NA", geneAnnotation$genes$symbol)
  gid <- geneAnnotation$genes$gene_id[!loci]
  df <- select(txdb, keys = gid, columns="TXNAME", keytype="GENEID")
  
  genes <- geneAnnotation$genes[!loci]
  exons <- geneAnnotation$exons[!grepl("^NA", geneAnnotation$exons$symbol)]
  tss <- geneAnnotation$TSS[which(geneAnnotation$TSS$tx_name %in% df$TXNAME)]
  
  ## Finalize annotation
  geneAnnotationSubset <- createGeneAnnotation(genes = genes, 
                                               exons = exons, 
                                               TSS = tss)
  genomeAnnotation$genome <- 'nullGenome'
  save(genomeAnnotation, geneAnnotation, geneAnnotationSubset, 
       file = smart_join(cur[[dir_val]],paste0(output_elements[2],"_ArchRef.RData")))
  }
}
