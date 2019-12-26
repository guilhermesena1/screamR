########################################################
## IO FUNCTIONS
########################################################

#' Reads mtx path as dgCMatrix
#' @export
readMM <- function(mtx.file) {
  if (!file.exists(mtx.file))
    stop(paste0("File not found: ", mtx.file))

  as(Matrix::readMM(mtx.fiel), 'dgCMatrix')
}

#' Reads 10X output matrix
#' @export
Read10X <- function(dir){
  dir.files <- list.files(dir)
  mtx.file <- grep(".mtx$",dir.files, value = T)
  if(length(mtx.file) > 1 )
    stop("Multiple mtx files in directory!")
  mtx.file <- as.character(mtx.file[1])

  mtx <- screamR::readMM(paste0(dir, "/",mtx.file))
  genes.file <- grep("^genes", dir.files, value=T)
  if(length(genes.file) > 1)
    stop("Multiple gene names files in directory!")

  # Read gene names if there is a genes file
  if(length(genes.file) == 1){
    genes <- readLines(paste0(dir,"/", genes.file))
    if(length(genes) != nrow(mtx))
      stop(paste0("genes file has ", length(genes), " genes, mtx has ", nrow(mtx),
                  "rows"))
    genes <- gsub("^.*_","",genes)
    genes <- gsub("^.*\t","",genes)
    rownames(mtx) <- genes

  } else {
    LogProcess("Warning: no genes.tsv file found, reading mtx w/o gene names");
  }

  colnames(mtx) <- paste0("cell_", 1:ncol(mtx))

  return(mtx)
}

########################################################
## MATRIX TRANSFORMATION FUNCTIONS
########################################################

#' Converts mtx file into a different set of gene symbols
#' Genes that do not originally exist in the matrix are treated as zeros across
#' all cells
#' @export
MtxFixGenes <- function(mtx, new.genes){
  mtx <- data.frame(as.matrix(mtx))
  mtx <- mtx[new.genes,]
  rownames(mtx) <- new.genes
  mtx[is.na(mtx)] <- 0
  mtx
}


