#  Copyright (C) 2019 University of Southern California and
#                Guilherme de Sena Brandine and Andrew D. Smith
#
#  Authors: Guilherme de Sena Brandine and Andrew D. Smith
#	 This program is free software: you can redistribute it and/or modify
#	 it under the terms of the GNU General Public License as published by
#	 the Free Software Foundation, either version 3 of the License, or
#	 (at your option) any later version.
#
#	 This program is distributed in the hope that it will be useful,
#	 but WITHOUT ANY WARRANTY; without even the implied warranty of
#	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	 GNU General Public License for more details.
#
#	 You should have received a copy of the GNU General Public License
#	 along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Methods to load, normalize and summarize data onto a single cell database
suppressMessages(library(Matrix))
suppressMessages(library(scales)) # for alpha function

##################################################################
######################### DB CONSTANTS ###########################
##################################################################

PATH.TO.SCREAM.DB <- "/auto/cmb-06/as/desenabr/scream/db" 

# Color palette to plot things nicely
SCREAM.COLORS <- c('#377eb8', '#ff7f00', '#4daf4a',
                   '#f781bf', '#a65628', '#384eff',
                   '#999999', '#e41a1c', '#dede00',
                   '#ffe119', '#000000', '#ffbea3',
                   '#911eb4', '#46f0f0', '#f032e6',
                   '#d2f53c', '#008080', '#e6beff',
                   '#aa6e28', '#800000', '#aaffc3',
                   '#808000', '#ffd8b1', '#000080',
                   '#808080', '#fabebe', '#a3f4ff')


##################################################################
######################### DB FUNCTIONS ###########################
##################################################################

#' Plots the dataset in reduced dimension
#'
#' Plots the reduced dimension database
#' res = any 2D dimension reduction (svd, tsne, umap)
#' samples = the discrete values to color points by
#' @param res the N x 2 reduced dimension, where N is the number of cells
#' @param samples optional, the samples to color by
#' @param points if false, points will become the first letter of the sample
#' name. Used when there are too many samples to color
#' @param file optional, a file path to save the plot as a pdf file, if
#' provided
#' @param main optional, the title of the plot
#' @param edges optional, if provided, connects the sample centroids by edges,
#' often used when representing transitions between clusters within the plot
#' @param plot.colors optional, the color palette to use to color points. 
#' @export
#' @examples
#' library(umap)
#' um <- umap(m.svd)
#' PlotReduction(um$layout, samples = cluster$classification)
PlotReduction <- function(res, 
                          samples = NULL, 
                          file = NULL, 
                          symbols = NULL,
                          points = T,
                          centroids = T,
                          width = 16,
                          height = 9,
                          main = "", 
                          edges = NULL,
                          edge.colors = NULL,
                          plot.colors = SCREAM.COLORS) {

  # If samples not provided, plot them the same color
  if(is.null(samples))
    samples <- factor(rep(1,nrow(res)))

  # Puts the label of each SRP in the most common location
  if(!is.factor(samples))
    samples <- factor(samples)

  if(centroids){
  cents <- NULL
    for(i in levels(samples)){
      if(sum(samples == i) > 1){
        # Finds the place with maximum density in 2d
        kernel <- MASS::kde2d(res[samples==i,1], res[samples==i,2])
        max.val <- which(kernel$z == max(kernel$z), arr.ind=T)

        cents <- rbind (cents, data.frame(x = kernel$x[max.val[1]], 
                                          y = kernel$y[max.val[2]]))
      }
      else {
        cents <- rbind (cents,res[samples == i,])
      }
    }
    colnames(cents) <- c("x","y")
    rownames(cents) <- levels(samples)
  }

  if(!is.null(file))
    pdf(file, width = width, height = height)

  # plot itself with no points (just xy limits)
  palette(plot.colors)
  plot(NA,
       xlim = c(min(res[,1]),max(res[,1])),
       ylim = c(min(res[,2]),max(res[,2])),
       xlab = "Dimension 1", ylab = "Dimension 2",
       main = main,
       xaxt = 'n', yaxt='n')

  # Points
  if(is.null(symbols))
    symbols.plot <- rep(19, nrow(res))
  else
    symbols.plot <- as.integer(factor(symbols))
  if(!points)
    text(res,
			 label = substring(as.character(samples),1,1), cex = .5, 
			 col = as.integer(samples))
  else
    points(res, cex = .5,
           pch = symbols.plot,
           col = as.integer(samples))

  # Tree edges
  if(!is.null(edges))
    for(i in 1:nrow(edges)){
      col <- "black"
      if(!is.null(edge.colors))
        col <- edge.colors[i]

      xs <- cents[edges[i,1],]
      xe <- cents[edges[i,2],]
      lines(c(xs[1], xe[1]), c(xs[2], xe[2]), col = col)
    }

  # Sample modes
  if(centroids){
    points(cents, col = "black", pch = 19, cex = 2.3)
    points(cents, col = "white", pch = 19, cex = 1.8)
    text(cents, col = "black", cex = 0.8, label = levels(samples))
  }

  # Legend of sample colors
  legend("topright", 
         legend = levels(samples), 
         col = 1:length(unique(samples)),
         pch = 19,
         cex = 1)

   # Legend of sample sym,bols
  if(!is.null(symbols))
    legend("bottomright", 
         legend = levels(factor(symbols)), 
         pch = 1:length(unique(symbols)),
         col = "black",
         cex = 1)
  

  if(!is.null(file))
    dev.off()

}
#' Plots a dataset gene in reduced dimension
#'
#' Plots a gene in the reduced dimension database
#' res = any 2D dimension reduction (svd, tsne, umap)
#' samples = the discrete values to color points by
#' @param res the N x 2 reduced dimension, where N is the number of cells
#' @param m.scale the scaled matrix from GeneScale, the values of which will be
#' used to color the points by gene expression
#' @param gene the gene to plot
#' @param file optional, a path to save the plot as a pdf file, if provided
#' @export
#' @examples
#' library(umap)
#' um <- umap(m.svd)
#' SCREAM.PlotGene(um$layout, m.scale, "Cd24")
PlotGene <- function(res,
                     m.scale,gene = NULL, 
                     symbols = NULL,
                     file = NULL, 
                     xlab = NULL,
                     ylab = NULL,
                     colors = c("lightgray","blue"),
                     width = 16, 
                     height =9) {
  # makes sure we plot a gene within the matrix
  if (!(gene %in% rownames(m.scale))) 
      stop(paste0("Gene ", gene, " is not in the matrix row names"))
  values <- m.scale[gene, ]

  # makes continuous colors for expression values
  crp <- colorRampPalette(colors)
  col.gradient <- crp(10)[as.numeric(cut(values, breaks = 10))]
  if (!is.null(file)) 
      pdf(file, width = width, height = height)
  xlab <- if (is.null(xlab)) 
      "Dimension 1"
  else 
       xlab
  ylab <- if (is.null(ylab)) 
      "Dimension 2"
  else 
       ylab
  symbols <- if(is.null(symbols))
      19
  else
      as.integer(symbols)
  plot(res, cex = 0.5, pch = symbols, 
      col = col.gradient)
  if (!is.null(gene)) 
      title(main = gene)
  if (!is.null(file)) 
      dev.off()

}


