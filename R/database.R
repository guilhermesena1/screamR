#  Copyright (C) 2019 University of Southern California and
#                Guilherme de Sena Brandine and Andrew D. Smith
#
#  Authors: Guilherme de Sena Brandine and Andrew D. Smith
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#     
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#      
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# This Rscript conormalizes the database and generates plots
#' Database conormalization
#' 
#' This function takes as an input the centroids from the database and the
#' metadata relative to each centroid and co-normalizes them based on phenotype. 
#' @param db the database mtx, with columns as database centroids and rows as
#' genes.
#' @param metadata the metadata read from the SCREAM database, with columns as
#' SRP, SRR, phenotype and additional observations 
#' @return The qsmooth co-normalized matrix
SCREAM.CoNormalize <- function(db, 
                               metadata){
  meta <- meta[!duplicated(meta$SRR),]
  rownames(meta) <- meta$SRR

  # Extract SRR names from db colnames
  # TODO: Replace with actual cluster
  srr <- str_extract(colnames(db), "SRR[0-9]*")

  # Check if all srrs are in the metadata:
  not.in.database <- which(!(srr %in% rownames(meta)))
  if(length(not.in.database) > 0){
    stop(paste0("please add the following SRRs to the metadata: ",
         paste(srr[not.in.database], collapse = " ")))
  }

  pheno <- meta[srr,]$Phenotypes
  LogProcess(paste0("Qsmooth on ", length(unique(pheno)), " phenotypes..."))
  db <- qsmooth(as.matrix(db), factor(pheno))

  writeMM(as(db,'sparseMatrix'), OUTPUT.CONORMALIZED.MATRIX)

  LogProcess("Scaling...")
  db <- GeneScale(db)

  LogProcess("SVD...")
  db.svd <- ReduceDimension(db)
  LogProcess(paste0("Number of db pcs: ", ncol(db.svd)))

  LogProcess("Running UMAP...")
  um <- umap(db.svd)
  SCREAM.PlotDB(um$layout, pheno, file = "SCREAM_UMAP.pdf")

  LogProcess("Running tSNE")
  ts <- tsne(db.svd)
  SCREAM.PlotDB(ts, file = "SCREAM_TSNE.pdf")

  LogProcess("Analysis finished successfully!")
}
#' Function to insert into the database
#' 
#' This function takes as an input an mtx file, analyzes it and inserts into
#' the database
#' with all centroids concatenated and 
#' @param IN.MTX.RAW.FILE the raw count data from an mtx file
#' @param IN.GENES.TSV.FILE the gene names of each row in the mtx file
#' @param DATABASE.PATH the directory path to which the analysis should be saved
#' in
#' @param srp the Sequencing Read Project (SRP) for the dataset
#' @param srr the Sequencing Read Run (SRR) for the dataset
#' @param sample.name An optional name for the dataset
#' @return NULL
#' @export
InsertIntoDatabase <- function(IN.MTX.RAW.FILE, 
                               IN.GENES.TSV.FILE,
                               DATABASE.PATH,
                               srp = NULL,
                               srr = NULL,
                               sample.name = NULL){

  ############# MAKE DIRECTORY FOR DATASET ###############
  if(is.null(srp))
    srp <- "SRP000000"
  if(is.null(srr))
    srr <- "SRR0000000"
  if(is.null(sample.name))
    sample.name <- ""

  OUT.DIR <- paste0(DATABASE.PATH, "/", paste0(srp, "_", srr))
  if(sample.name != "")
    OUT.DIR <- paste0(OUT.DIR, "_", sample.name)

  OUT.MTX.RAW.FILE <- paste0(OUT.DIR, "/mtx_raw.mtx")
  OUT.GENES.TSV.FILE <- paste0(OUT.DIR, "/genes.tsv")
  OUT.MTX.CENTROIDS.FILE <- paste0(OUT.DIR, "/mtx_centroids.mtx")
  OUT.MTX.PROBS.FILE <- paste0(OUT.DIR, "/mtx_probs.mtx")
  OUT.MCLUST.RDS.FILE <- paste0(OUT.DIR, "/mclust.rds")
  OUT.ANALYSIS.STATS.FILE <- paste0(OUT.DIR, "/analysis.stats")
  OUT.CLUSTER.ASSIGN.FILE <- paste0(OUT.DIR, "/clusters.tsv")

  ############# RUN PIPELINE ################################
  # Reading and filtering matrix
  LogProcess(paste0("Reading ", IN.MTX.RAW.FILE, "..."))
  m.raw <- readMM(IN.MTX.RAW.FILE)

  ############# FIX FOR SOMETHING I DID WRONG: ###############
  # My reference transcriptome has sequences for EGFP and tdTomato. I want to
  # remove these from any count data
  ###########################################################
  LogProcess("Filtering out EGFP and tdTomato")
  artificial.genes <-  c("EGFP","tdTomato")
  genes <- readLines(IN.GENES.TSV.FILE)

  # In case it's a 10x with ENSG/ENSMUSG as first column
  # if it's not then it doesn't do anything
  genes <- gsub("^.*\t", "", genes)
  rownames(m.raw) <- genes
  m.raw <- m.raw[!(genes %in% artificial.genes),]

  # Filtering and printing info
  LogProcess(paste0("Matrix dimension before: ", nrow(m.raw), " x ", ncol(m.raw)))

  csums <- Matrix::colSums(m.raw)
  rsums <- Matrix::rowSums(m.raw)
  m.raw <- m.raw[rsums > 50, csums > 1000]

  LogProcess(paste0("Matrix dimension after: ", nrow(m.raw), " x ", ncol(m.raw)))

  # Normalization
  LogProcess(paste0("Log normalizing ", IN.MTX.RAW.FILE, "..."))
  m <- log(1 + ColumnScale(m.raw, norm.factor = 1e4))
  LogProcess(paste0("Scaling ", IN.MTX.RAW.FILE, "..."))
  m <- GeneScale(m)

  # Computes eigenvalues to infer how many dimensions to reduce to
  LogProcess("Reducing Dimension...")
  m <- ReduceDimensionTruncated(m, nv = min(100, min(nrow(m) - 1,ncol(m) - 1)))
  LogProcess(paste0("PCs to keep: ", ncol(m)))

  # Clustering
  library(mclust)
  init.cluster <- hc(m)
  res <- Mclust(m, G = 1:50, modelNames="VVV",
                initialization = list(hcPairs = init.cluster,
                                      use = "VARS",
                                      modelName = "VVV"))

  LogProcess(paste0("Number of clusters: ", res$G))
  LogProcess(paste0("Model for best clustering: ", res$modelName))
  LogProcess(paste0("BIC: ", as.numeric(res$bic)))
  LogProcess("Summarizing number of cells in each cluster:")
  print(summary(factor(res$classification)))

  #################################################################
  ################ ANALYSIS DONE, NOW INSERT  #####################
  #################################################################


  ################ CREATE DIR  ########################
  LogProcess("Creating directory if it does not exist")
  if(!dir.exists(OUT.DIR))
    dir.create(OUT.DIR, mode = "0777")

  ################ WRITE MTX RAW  #####################
  writeMM(m.raw, OUT.MTX.RAW.FILE)

  ################ WRITE GENES ########################
  writeLines(rownames(m.raw), OUT.GENES.TSV.FILE)

  ################ WRITE PROBABILITIES ################
  LogProcess("Calculating centroids in high dimension")
  res$z[res$z < 1e-4] <- 0.0
  res$z[res$z > .999] <- 1.0
  silent.return <- writeMM(as(res$z, "sparseMatrix"), OUT.MTX.PROBS.FILE)

  # Finds the centroids in high dimension by multiplying the cluster assignment
  # probabilities to the cell counts to get the average profile of the cluster
  # Zeros very small probability numbers to render centroids sparse. We don't lose
  # much by discarding low probabilities and save a lot of computational time in
  # matrix multiplications

  # The centroids we are saving is a 'pseudobulk' in which we add the counts for
  # all cells that belong to a cluster. However, counts may not necessarily be
  # discrete, as some cluster probability assignments are not 0 or 1. 
  m.centroids <- as.matrix(m.raw) %*% res$z
  m.centroids[m.centroids < 1e-2] <- 0

  ############## WRITE CENTROIDS ################
  m.centroids <- as(m.centroids, "sparseMatrix")
  silent.return <- writeMM(m.centroids, OUT.MTX.CENTROIDS.FILE)

  ############## WRITE MCLUST OBJ ################
  saveRDS(res, file = OUT.MCLUST.RDS.FILE)

  ############## WRITE CLUSTERING ASSIGNMENTS####
  writeLines(as.character(res$classification), OUT.CLUSTER.ASSIGN.FILE)

  ############## WRITE ANALYSIS STATISTICS################
  analysis.stats <- data.frame(
    cluster_number = res$G,
    cluster_bic = as.numeric(res$bic),
    cluster_model = res$modelName,
    cluster_entropy = -sum(res$parameters$pro*log(res$parameters$pro)) / log(res$G),
    cluster_separability = ClusterSeparation(m, res),
    intrinsic_dimension = ncol(m))

  write.table(t(analysis.stats), 
              file = OUT.ANALYSIS.STATS.FILE, 
              quote=F,sep="\t",
              col.names=F, row.names=T)
  ######################################################
  LogProcess("Successfully added new dataset into the database!")
}
