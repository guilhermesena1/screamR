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

#########################################################
#### CUSTOM FUNCTIONS USED THROUGHOUT THE ANALYSIS
#########################################################
#' reimplementing size factor estimation to avoid DESeq dependency
estimateSizeFactorsForMatrix <- function(counts, locfunc=stats::median) {
  loggeomeans <- apply(log(counts),1,mean)
  if (all(is.infinite(loggeomeans)))
    return (NULL)

  sf <-apply(counts, 2, 
    function(cnts) {
      exp(locfunc((log(cnts) - loggeomeans)
                  [is.finite(loggeomeans) & 
                   cnts > 0
                  ])
         )
    })

  sf/exp(mean(log(sf)))
}

#' Converts count data into RPX
#'
#' This function attempts to size factor normalize the raw counts. If no genes
#' are expressed across all cells, size factors are approximated by library
#' size. 
#' @param m.raw the raw data matrix
#' @param size.factors optional, if you know that size factors cannot be
#' estimated, size factor detection can be skipped
#' @return The log-normalized matrix
#' @export
#' @examples
#' m.raw <- matrix(rpois(1000, lambda = 10), ncol = 20)
#' m <- LogNormalize(m)
ColumnScale <- function(m.raw, size.factors = T, norm.factor = NULL){
  # Number of reads per cell, relative to the median
  sf <- NULL
  if(size.factors & is.null(norm.factor)){
    LogProcess("Estimating size factors...")
    sf <- estimateSizeFactorsForMatrix(m.raw)
    if(is.null(sf))
      LogProcess("Not enough nonzero genes to estimate size factors.",
                 "Using library size instead")
    else
      return(Matrix::t(Matrix::t(m.raw)/sf))
  }

  # size factors as library size
  sf <- ColSums(m.raw)
  if(is.null(norm.factor))
    sf <- sf / median(sf)
  else
    sf <- sf / norm.factor

  Matrix::t(Matrix::t(m.raw)/sf)
}

#' Gene-wise scaling of log-normalized data
#'
#' Subtracts the row mean and divides by the row standard deviation.
#' @param m the normalized matrix (eg: from LogNormalize)
#' @return The scaled matrix, where each row has mean = 0 and sd = 1
#' @export
#' @examples
#' m.raw <- matrix(rpois(1000, lambda = 10), ncol = 20)
#' m <- LogNormalize(m)
#' m.scale <- GeneScale(m)
GeneScale <- function(m){
  # Cannot divide by zero standard deviation
  sds <- apply(m,1,sd)
  m <- m[sds>0,]

  # Scale by row, ie, subtract each gene by mean and divide by gene sd
  m <- Matrix::t(scale(Matrix::t(m)))

  return(m)
}
#' Reduces the matrix to SVD space
#'
#' Performs svd and keeps only the dimensions whose eigenvector is larger than
#' the theoretical maximum of the same matrix with standard normal values
#' @param m.scale the scaled matrix (eg: from GeneScale)
#' @return The matrix product V * Sigma in the Singular Value Decomposition of
#' m.scale
#' @export
#' @examples
#' m.raw <- matrix(rpois(1000, lambda = 10), ncol = 20)
#' m <- LogNormalize(m)
#' m.scale <- GeneScale(m)
#' m.svd <- ReduceDimenson(m.scale)
ReduceDimension <- function (m.scale){
  M <- ncol(m.scale)
  N <- nrow(m.scale)
  m.svd <- gmodels::fast.svd(m.scale)

  eigenvalues <- m.svd$d/sqrt(N)

  # Eigenvalues larger than the maximum of the MP distribution
  p.values <- 1 - RMTstat::pmp(eigenvalues, ndf = N, pdim = M)
  keep <- which(p.values < 0.01)

  # Keep at least 2 dimensions
  if(length(keep) < 2)
    keep <- c(1,2)

  return(m.svd$v[,keep] %*% diag(m.svd$d[keep]))
}
#' Reduces the matrix to SVD space to a fixed number of dimensions
#'
#' Performs svd and keeps only the dimensions whose eigenvector is larger than
#' the theoretical maximum of the same matrix with standard normal values.
#' Unlike ReduceDimension, this function is a faster version when a 'guess' on
#' number of eigenvalues is given. If the number of significant eigenvalues is
#' smaller than the guess, then this function should be used with significant
#' performance and no error, otherwise a warning is printed to increase the
#' guess size. 
#'  
#' @param m.scale the scaled matrix (eg: from GeneScale)
#' @return The matrix product V * Sigma in the Singular Value Decomposition of
#' m.scale
#' @export
#' @examples
#' m.raw <- matrix(rpois(1000, lambda = 10), ncol = 20)
#' m <- LogNormalize(m)
#' m.scale <- GeneScale(m)
#' m.svd <- ReduceDimensonTruncated(m.scale, N = 100)
ReduceDimensionTruncated <- function (m.scale, nv = 100){
  M <- ncol(m.scale)
  N <- nrow(m.scale)
  m.svd <- irlba::irlba(m.scale, nv = nv)

  eigenvalues <- m.svd$d/sqrt(N)

  # Eigenvalues larger than the maximum of the MP distribution
  p.values <- 1 - RMTstat::pmp(eigenvalues, ndf = N, pdim = M)
  keep <- which(p.values < 0.01)

  # Prints a warning if all principal components are significant, which means
  # more PCs should be calculated
  if(length(keep) == nv){
    LogProcess(paste0("[WARNING] All PCs significant, consider increasing",
                      "the number of  PCs!"))

    mp.eig <- (1 + sqrt(ncol(m.scale)/nrow(m.scale)))^2
    LogProcess("Normalized eigenvalues: ", paste(round(eigenvalues,3),
                                           collapse=" ")) 
    LogProcess("MP eigenvalue:", mp.eig)
  }

  # Keep at least 2 dimensions
  if(length(keep) < 2)
    keep <- c(1,2)
  return(m.svd$v[,keep] %*% diag(m.svd$d[keep]))
}

##########################################################################
### CLUSTERING AND AUX CLUSTERING FUNCTIONS
##########################################################################
#' Clusters the single cells using gaussian mixture model for a fixed value of k
#' 
#' @param m.svd The reduced dimension dataset, rows as cells columns as PCs
#' @param k The number of clusters to run EM on 
#' @param init the Mclust::hc initialization function
#' @return an mclust object with the probabilistic clustering results
#' @return An mclust object with clustering results: Classification,
#' probabilities for each cell and mean/variance/probability parameters for each
#' cluster. 
#' @export
#' @examples
#' m.raw <- matrix(rpois(1000, lambda = 10), ncol = 20)
#' m <- LogNormalize(m)
#' m.scale <- GeneScale(m)
#' m.svd <- ReduceDimenson(m.scale)
#' clusters <- ClusterWithFixedK(m.svd, k = 3)
ClusterWithFixedK <- function(m.svd, k, init = NULL,verbose=T){
  if(is.null(init)){
    if(verbose)
      LogProcess("Computing initial guess based on hierarchical clustering...")
    init <- hc(m.svd)
  }

  ans <- mclust::Mclust(m.svd, 
                 G = k,
                 initialization = list(hcPairs = init,
                                        use = "VARS",
                                        modelName="VVV"),
                                        verbose=verbose)

  if(is.null(ans))
    return (NA)

  return (ans)
}
#' Clusters the single cells using gaussian mixture model
#' 
#' Runs Mclust the way it is supposed to be used in the scRNA-Seq context
#' @param m.svd The reduced dimension dataset, rows as cells columns as PCs
#' @param N (temporary) the maximum number of clusters to try
#' @return an mclust object with the probabilistic clustering results
#' @return An mclust object with clustering results: Classification,
#' probabilities for each cell and mean/variance/probability parameters for each
#' cluster. 
#' @export
#' @examples
#' m.raw <- matrix(rpois(1000, lambda = 10), ncol = 20)
#' m <- LogNormalize(m)
#' m.scale <- GeneScale(m)
#' m.svd <- ReduceDimenson(m.scale)
#' clusters <- ClusterCells(m.svd)
ClusterCells <- function(m.svd, N = 30,verbose=T){
  # Number of clusters cannot be higher than number of cells
  if(N > nrow(m.svd)){
    LogProcess("WARNING: Scaling down max number of clusters to number of cells")
    N <- nrow(m.svd)
  }

  # This instructs which cells to merge based on hc
  LogProcess("Computing initial guess based on hierarchical clustering...")
  init.cluster <- hc(m.svd)
 
  # Initialize list lengths
  cluster.attempts <- as.list(rep(NA, N))
  separations <- rep(NA, N)

  for(i in 1:N){
    if(verbose)
      LogProcess(paste0("i = ",i))

    # Individual clustering attempt
    cluster.attempts[[i]] <- ClusterWithFixedK(m.svd,
                                               k = i,
                                               init = init.cluster,
                                               verbose)

    # Separation measure for this attempt
    separations[i] <-  ClusterSeparation(m.svd, cluster.attempts[[i]])
    if(verbose)
      LogProcess(paste0("separation for i = ",i,": ", separations[i]))
  }

  separations <- smooth.spline(1:N, separations, spar = 0.5)$y
  first.derivative <- separations[2:N] - separations[1:(N-1)]
  second.derivative <- first.derivative[2:(N-1)] - first.derivative[1:(N-2)]
  curvature <- second.derivative/((1 + first.derivative[1:(N-2)]^2)^1.5)
  best.g <- which.min(separations)

  return(list(
              mclust = cluster.attempts[[best.g]],
              separations = separations,
              curvature = curvature))
}
#' Calcualtes a measure of separability
#' 
#' Function to calculate cluster separation, ie, the ratio between the variance
#' within clusters and the variance between clusters. A smaller value means the
#' cluster is more separated. Note that separation > 1 means that the variance
#' within clusters is higher than the variance between, and this is an
#' indication that the dataset was overclustered. 
#' @param m.svd the reduced dimension form the ReduceDimension function
#' @param clusters the mclust object with probabilistic cluster assignments
#' @return A value between 0 (clusters fully separated) and infinity
#' @export
#' @examples
#' m.raw <- matrix(rpois(1000, lambda = 10), ncol = 20)
#' m <- LogNormalize(m)
#' m.scale <- GeneScale(m)
#' m.svd <- ReduceDimenson(m.scale)
#' clusters <- ClusterCells(m.svd)
#' ClusterSeparation(m.svd,clusters)
ClusterSeparation <- function(m.svd, clusters){
  if(clusters$G == 1)
    return (1)

  if(any(is.na(clusters$classification)))
    return(1)

  n <- nrow(m.svd)
  # pairwise distance matrix
  D <- as.matrix(dist(m.svd))^2

  # Cell cluster probability
  P <- clusters$z

  # separation matrix
  S <- t(P) %*% D %*% P

  # Avg distance between cells in the same cluster
  distance.within <- mean(diag(S))

  # Avg distance between cells in different clusters
  distance.between <- mean(S[upper.tri(S)])

  return (distance.within / (distance.within + distance.between))
}

##########################################################################
### DIFFERENTIAL EXPRESSION AND MARKER FINDING FUNCTIONS
##########################################################################
#' Student's t-test to find differentially expressed genes on each cluster
#' 
#' This function uses the scaled values calculated by GeneScale to run student's
#' t-test and find genes that mark certain clusters
#' @param m.scale the N x M scaled normalized count matrix given by GeneScale
#' @param clusters the clustering result from ClusterCells with k clusters
#' @param which.cluster the cluster to be tested
#' @return An N x (2k - 2) table, showing the log fold change of each gene in
#' the tested cluster vs all other cluster, with respective p-values
#' @export
#' @examples
#' m.raw <- matrix(rpois(1000, lambda = 10), ncol = 20)
#' m <- LogNormalize(m)
#' m.scale <- GeneScale(m)
#' m.svd <- ReduceDimenson(m.scale)
#' clusters <- ClusterCells(m.svd)
#' # Finds DE genes for cluster 1 vs 2
#' tbl <- DiffExpr(m.scale, clusters$classification, 1,2)
DiffExpr <- function(m.scale, clusters,cl1,cl2){
  if(sum(cl1 %in% clusters) == 0)
    stop(paste0("Cluster ", cl1," is not a cluster",
                " in the given list"))
  if(sum(cl2 %in% clusters) == 0)
    stop(paste0("Cluster ", cl2," is not a cluster",
                " in the given list"))

  # Repeat for every gene
  t(apply(m.scale, 1, function(x){
    # Expressions in cluster A
    case <- x[clusters == cl1]

    # Expressions in any other cluster
    rest <- x[clusters == cl2]

    # Log fold-change
    logfc <- mean(case) - mean(rest)

    #p-value
    p <- 1

    # If logfc <= 0, then it's rejected no need to run t test
    if(logfc > 0.0)
      if(sd(case) > 0 | sd(rest) > 0)
        p <- t.test(case,rest,alternative = "greater")$p.value

    ans <- c(logfc,p)
    names(ans) <- c(paste0("logfc_",cl1,"_vs",cl2), paste0("p_",cl1,"_vs_",cl2))
    ans
  }))
}
#' Student's t-test to find differentially expressed genes on each cluster
#' 
#' This function uses the scaled values calculated by GeneScale to run student's
#' t-test and find genes that mark certain clusters
#' @param m.scale the N x M scaled normalized count matrix given by GeneScale
#' @param clusters the clustering result from ClusterCells with k clusters
#' @param which.cluster the cluster to be tested
#' @return An N x (2k - 2) table, showing the log fold change of each gene in
#' the tested cluster vs all other cluster, with respective p-values
#' @export
#' @examples
#' m.raw <- matrix(rpois(1000, lambda = 10), ncol = 20)
#' m <- LogNormalize(m)
#' m.scale <- GeneScale(m)
#' m.svd <- ReduceDimenson(m.scale)
#' clusters <- ClusterCells(m.svd)
#' # Finds DE genes for cluster 1
#' tbl <- DiffExprAll(m.scale, clusters$classification, 1)
DiffExprAll <- function(m.scale, clusters, which.cluster){
  if(sum(which.cluster %in% clusters) == 0)
    stop(paste0("Cluster ", which.cluster," is not a cluster",
                " in the Mclust result"))

  # All clusters that are not the test cluster
  others <- unique(clusters)
  others <- others[others != which.cluster]

  # Repeat for every gene
  t(apply(m.scale, 1, function(x){

    # Expressions in cluster A
    case <- x[clusters == which.cluster]
    ans <- c()
    for(i in others){
      rest <- x[clusters == i]
      logfc <- mean(case) - mean(rest)
      p <- 1

      # If logfc <= 0, then it's rejected no need to run t test
      if(logfc > 0.0)
        if(sd(case) > 0 | sd(rest) > 0)
          p <- t.test(case,rest,alternative = "greater")$p.value

      ans <- c(ans,logfc,p)
    }

    names(ans) <- rep(c("logfc_","p_"), length(others))
    append <- c()
    for(i in others)
      append <- c(append,c(i,i))
    names(ans) <- paste0(names(ans), append)
    ans
  }))
}
#' Finds marker genes from differential expression results
#' 
#' Given a test from DiffExpr, returns the genes which are significantly larger
#' in a cluster vs all other clusters, these are considered to be marker genes. 
#' @param diffexpr the result from the DiffExpr function
#' @param sig.level the maximum p-value to be considered significantly DE. Genes
#' where p < sig.level in every cluster will be returned. 
#' @return A list of genes given by the m.scale row names
#' @export
#' @examples
#' m.raw <- matrix(rpois(1000, lambda = 10), ncol = 20)
#' m <- LogNormalize(m)
#' m.scale <- GeneScale(m)
#' m.svd <- ReduceDimenson(m.scale)
#' clusters <- ClusterCells(m.svd)
#' # Finds DE genes for cluster 1
#' tbl <- DiffExpr(m.scale, clusters$classification, 1)
#' # Finds the markers of cluster 1 with p < 0.05 on all tests
#' markers <- GetMarkers(tbl, sig.level = 0.05)
GetMarkers <- function(diffexpr, sig.level = 0.05){
  diffexpr <- diffexpr[, grepl("^p_",colnames(diffexpr))]
  keep <- apply(diffexpr, 1,function(x){
          !any(x > sig.level) 
  })
  return(rownames(diffexpr)[keep])
} 

##########################################################################
### PAIRWISE DISTANCE FUNCTIONS
##########################################################################
#' Averages the reduced dimension dataset by a given feature
#' 
#' @param m the N x M count matrix
#' @param group the vector of length M with k factors
#' @return the N x k matrix with cells grouped by the group feature (eg: sum or
#' mean)
#' @export
GroupByFeature <- function(m, group, func = sum){
  if(length(group) != ncol(m))
    stop("[ERROR] Number of group factors different than number of matrix rows")

  sapply(unique(group), function(x){
         apply(m[,group == x], 1,func)
  })
}

#' Ratio between distance of matching clusters and distance between clusters
#' within the same dataset
#' 
#' @param dist.matrix the pairwise distance between clusters
#' @param pheno.dataset.tbl a table with two columns: phenotype and dataset of
#' origin
#' @return the k x d averages of cells from each group
#' @export
WithinBetweenRatio <- function(dist.matrix, pheno.dataset.tbl){
  same.cluster.dist <- 0
  within.dataset.dist <- 0
  ns <- 0
  nw <- 0
  for(i in 1:nrow(dist.matrix)){
    for(j in 1:nrow(dist.matrix)){
      if(i != j){
        pa <- pheno.dataset.tbl[i,1]
        pb <- pheno.dataset.tbl[j,1]
        da <- pheno.dataset.tbl[i,2]
        db <- pheno.dataset.tbl[j,2]

        if((pa == pb) & (da != db)){
          ns <- ns + 1
          same.cluster.dist <- same.cluster.dist + dist.matrix[i,j]
        }
        nw <- nw + 1
        within.dataset.dist <- within.dataset.dist + dist.matrix[i,j]
      }
    }
  }
  same.cluster.dist <- same.cluster.dist/ns
  within.dataset.dist <- within.dataset.dist/nw

  same.cluster.dist/within.dataset.dist
}


##########################################################################
### PHENOTYPE INFERENCE
##########################################################################
#' Hypergeometric test if marker genes significantly resemble a given phenotype
#' 
#' Outputs the probability that a set of K marker genes have k genes in common
#' with a set of n phenotype genes in a universe of N total genes
#' @param all.genes a set of N genes representing the ones tested for DE
#' @param pheno.genes a set of n genes, contained in all genes, that represent
#' some phenotype
#' @param marker.genes a set of K genes, contained in all.genes 
#' @return A list of genes given by the m.scale row names
#' @export
#' @examples
#' m.raw <- matrix(rpois(1000, lambda = 10), ncol = 20)
#' m <- LogNormalize(m)
#' m.scale <- GeneScale(m)
#' m.svd <- ReduceDimenson(m.scale)
#' clusters <- ClusterCells(m.svd)
#' # Finds DE genes for cluster 1
#' tbl <- DiffExpr(m.scale, clusters$classification, 1)
#' # Finds the markers of cluster 1 with p < 0.05 on all tests
#' markers <- GetMarkers(tbl, sig.level = 0.05)
HyperTest <- function(all.genes, pheno.genes, marker.genes){
  if(any(!(pheno.genes %in% all.genes))){
     bad.genes <- pheno.genes[which(!(pheno.genes %in% all.genes))]
     frac.bad <- 100 * length(bad.genes) / length(pheno.genes)

  #   LogProcess("[WARNING] Some of the  phenotype genes are not in the",
  #              "annotation. We will be working with the good genes only.")
  #   LogProcess("Fraction of bad genes: ", frac.bad)

     if(frac.bad == 100)
       return (NULL)
     pheno.genes <- intersect(pheno.genes, all.genes)
  }

  N <- length(all.genes)
  n <- length(pheno.genes)
  K <- length(marker.genes)
  hits <- intersect(marker.genes,pheno.genes)
  k <- length(hits) 

  # The minus one is because phyper is p <= x, so if there are zero the p-value
  # should be 1
  return (list (num.hits = k,
                hits = paste(hits, collapse=","),
                p = 1 - phyper(q = k - 1, m = n, n = N - n, k = K)))
}

#' Reads file downloaed from http://bio-bigdata.hrbmu.edu.cn/CellMarker
#' into a list of phenotypes, whose names are phenotypes and values are vectors
#' marker genes
#' @param the filename of a CellMarker table, eg
#' /path/to/db/mm10/metadata/markers.tsv
#' @return a list with phenotypes as names and genes as values
#' @export
ReadCellMarkersIntoList <- function(markers.file){
  df <- read.table(markers.file, sep="\t", quote= "\"", header=1, row.names=NULL) 

  ans <- list()
  for(i in 1:nrow(df)){
    row <- df[i,]
    pheno <- paste(c(as.character(row$tissueType),
                   as.character(row$UberonOntologyID),
                   as.character(row$cancerType),
                   as.character(row$cellType),
                   as.character(row$cellName),
                   as.character(row$cellOntologyId)),
                   collapse = ".")

    ans[[pheno]] <- do.call(c, strsplit(as.character(row$geneSymbol),
                             split = ", "))
  }

  return (ans)
}

#' Reads file downloaded from https://panglaodb.se/markers.html
#' into a list of phenotypes, whose names are phenotypes and values are vectors
#' marker genes
#' @param the filename of a Panglao table, eg
#' /path/to/db/mm10/metadata/panglao.tsv
#' @return a list with phenotypes as names and genes as values
#' @export
ReadPanglaoIntoList <- function(markers.file){
  df <- read.table(markers.file, sep="\t", quote="\"",
                   header=1, row.names=NULL) 

  ans <- list()
  for(i in 1:nrow(df)){
    row <- df[i,]
    pheno <- paste(c(as.character(row$organ),
                     as.character(row$pheno),
                     as.character(row$germlayer)),
                   collapse = ".")

    # TODO: Add aliases
    genes <- as.character(row$symbol)

    if(!is.null(ans[[pheno]]))
      ans[[pheno]] <- c(ans[[pheno]], genes)
    else
      ans[[pheno]] <- genes
  }

  return (ans)
}


#' Given a set of markers, guesses the phenotype based on smallest p value 
#' of a set of hypergeometric tests from a phenotype list
#' @export
GuessPhenotype <- function(markers, all.genes, pheno.list, verbose=T){
  p.list <- NULL
  best.p <- 1
  pheno <- NULL
  for(i in names(pheno.list)){
    test <- HyperTest(all.genes, pheno.list[[i]], markers)
    if(!is.null(test)){
      p.list <- rbind(p.list, data.frame(num.hits = test$num.hits,
                                         p.value = test$p,
                                         hits = test$hits,
                                         row.names = i))

      if(test$p < best.p){
        best.p <- test$p
        pheno <- i
      }
    } else {
      if(verbose)
        LogProcess(i,"is a bad phenotype")
    }
  }

  if(best.p > 0.01)
    pheno <- "UNKNOWN"

  return(list (pheno = pheno, best.p = best.p, full.table = p.list))
}

##########################################################################
### PATHWAY ANALYSIS
##########################################################################

#' @export
ReadPathwaysIntoList <- function(file.path){
  tbl <- read.table(file.path, 
                    sep="\t",
                    header=F, 
                    row.names=NULL,
                    stringsAsFactors=F)
  ans <- list()
  for(i in 1:nrow(tbl)){
    ans[[tbl[i,1]]] <- strsplit(tbl[i,2],",")[[1]]
  }

  ans
}

#' Makes an expression matrix of pathways 
#' @param count.mtx the n x m normalized matrix (eg from GeneScale)
#' @param pathway.list the list file with k pathways as names and gene lists as
#' @param verbose print more run info
#' values, eg: from ReadPathwaysIntoList function
#' @return a k x m matrix with pathway scores for each cell
#' @export
MakePathwayMatrix <- function(count.mtx, pathway.list,verbose=F){
  ans <- NULL
  for(i in names(pathway.list)){
    if(verbose)
      LogProcess("Calculating eigengene for pathway",i,"...")

    p.genes <- pathway.list[[i]]
    p.genes <- p.genes[p.genes %in% rownames(count.mtx)]
    tmp <- count.mtx[p.genes,]
    if(length(p.genes)>1){
      eigengenes <- svd(tmp, nv = 1)$v
      ans <- cbind(ans, eigengenes)
    } else {
      ans <- cbind(ans,tmp)
    }
  }
  ans <- t(ans)
  LogProcess("Dim: ", paste(dim(ans), collapse = " x "))
  rownames(ans) <- names(pathway.list)
  colnames(ans) <- colnames(count.mtx)
  ans
}
