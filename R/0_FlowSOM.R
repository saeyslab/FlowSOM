#' Run the FlowSOM algorithm
#'
#' Method to run general FlowSOM workflow. 
#' Will scale the data and uses consensus meta-clustering by default.
#'
#' @param input         a flowFrame, a flowSet, a matrix with column names or
#'                      an array of paths to files or directories
#' @param pattern       if input is an array of file- or directorynames, select 
#'                      only files containing pattern
#' @param compensate    logical, does the data need to be compensated
#' @param spillover     spillover matrix to compensate with
#'                      If NULL and compensate = TRUE, we will look for $SPILL 
#'                      description in FCS file.
#' @param transform     logical, does the data need to be transformed with the
#'                      transformation given in \code{transformFunction}.
#' @param toTransform   column names or indices that need to be transformed.
#'                      Will be ignored if \code{transformList} is given.
#'                      If \code{NULL} and transform = \code{TRUE}, column names
#'                      of \code{$SPILL} description in FCS file will be used.
#' @param transformFunction Defaults to logicleTransform()
#' @param transformList transformList to apply on the samples.
#' @param scale         logical, does the data needs to be rescaled. 
#'                      Default = FALSE
#' @param scaled.center see \code{\link{scale}}
#' @param scaled.scale  see \code{\link{scale}}
#' @param silent        if \code{TRUE}, no progress updates will be printed
#' @param colsToUse     Markers, channels or indices to use for building the SOM. 
#'                      Default (NULL) is all the columns used to build the 
#'                      FlowSOM object.
#' @param importance    array with numeric values. Parameters will be scaled 
#'                      according to importance
#' @param nClus         Exact number of clusters for meta-clustering. 
#'                      Ignored if maxMeta is specified.
#'                      Default = 10.
#' @param maxMeta       Maximum number of clusters to try out for 
#'                      meta-clustering. If \code{NULL} (default), only one 
#'                      option will be computed (\code{nClus}).
#' @param seed          Set a seed for reproducible results
#' @param ...           options to pass on to the SOM function 
#'                      (xdim, ydim, rlen, mst, alpha, radius, init, distf)
#'
#' @return A \code{list} with two items: the first is the flowSOM object 
#'         containing all information (see the vignette for more detailed 
#'         information about this object), the second is the metaclustering of 
#'         the nodes of the grid. This is a wrapper function for 
#'         \code{\link{ReadInput}}, \code{\link{BuildSOM}}, 
#'         \code{\link{BuildMST}} and \code{\link{MetaClustering}}. 
#'         Executing them separately may provide more options.
#'
#' @seealso \code{\link{scale}}, 
#'          \code{\link{ReadInput}}, 
#'          \code{\link{BuildSOM}},
#'          \code{\link{BuildMST}}, 
#'          \code{\link{MetaClustering}}
#' @examples
#' # Read from file
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- FlowSOM(fileName, compensate = TRUE, transform = TRUE,
#'                       scale = TRUE, colsToUse = c(9, 12, 14:18), nClus = 10)
#' # Or read from flowFrame object
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
#' ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(flowCore::keyword(ff)[["SPILL"]]),
#'                                 flowCore::logicleTransform()))
#' flowSOM.res <- FlowSOM(ff, 
#'                        scale = TRUE, 
#'                        colsToUse = c(9, 12, 14:18), 
#'                        nClus = 10)
#' 
#' # Plot results
#' PlotStars(flowSOM.res,
#'           backgroundValues = flowSOM.res$metaclustering)
#' 
#' # Get metaclustering per cell
#' flowSOM.clustering <- GetMetaclusters(flowSOM.res)
#' 
#' @importFrom BiocGenerics colnames
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @importFrom flowCore read.FCS compensate transform logicleTransform exprs 
#'             transformList write.FCS 'exprs<-' keyword fr_append_cols
#' @importFrom igraph graph.adjacency minimum.spanning.tree layout.kamada.kawai
#'             plot.igraph add.vertex.shape get.edges shortest.paths E V 'V<-'
#'             igraph.shape.noclip
#' @importFrom rlang .data
#' @importFrom stats prcomp
#' @importFrom utils capture.output packageVersion
#' @importFrom XML xmlToList xmlParse
#' @importFrom dplyr group_by summarise_all select
#' @importFrom stats median
#' 
#' @export
FlowSOM <- function(input, pattern = ".fcs", 
                    compensate = FALSE, spillover = NULL, 
                    transform = FALSE, toTransform = NULL, 
                    transformFunction = flowCore::logicleTransform(), 
                    transformList = NULL, scale = FALSE, 
                    scaled.center = TRUE, scaled.scale = TRUE, silent = TRUE, 
                    colsToUse = NULL, nClus = 10, maxMeta = NULL, importance = NULL, 
                    seed = NULL, ...){
  # Method to run general FlowSOM workflow. 
  # Will scale the data and uses consensus meta-clustering by default.
  #
  # Args:
  #    input: dirName, fileName, array of fileNames, flowFrame or 
  #           array of flowFrames
  #    colsToUse: column names or indices to use for building the SOM
  #    maxMeta: maximum number of clusters for meta-clustering
  #
  # Returns:
  #    list with the FlowSOM object and an array with final clusterlabels
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  t <- system.time(fsom <- ReadInput(input, pattern = pattern, 
                                     compensate = compensate, 
                                     spillover = spillover, 
                                     transform = transform, 
                                     toTransform = toTransform, 
                                     transformFunction = transformFunction, 
                                     transformList = transformList,
                                     scale = scale,
                                     scaled.center = scaled.center, 
                                     scaled.scale = scaled.scale, 
                                     silent = silent))
  if(!silent) message(t[3], "\n")
  t <- system.time(fsom <- BuildSOM(fsom, colsToUse, silent = silent, 
                                    importance = importance, ...))
  if(!silent) message(t[3], "\n")
  t <- system.time(fsom <- BuildMST(fsom, silent = silent))
  if(!silent) message(t[3], "\n")
  if(is.null(maxMeta)){
    t <- system.time(cl <- as.factor(
      metaClustering_consensus(fsom$map$codes, nClus, seed = seed)))
  } else {
    t <- system.time(cl <- as.factor(MetaClustering(fsom$map$codes,
                                                    "metaClustering_consensus", 
                                                    maxMeta,
                                                    seed = seed)))
  }
  fsom$map$nMetaclusters <- length(levels(cl))
  fsom$metaclustering <- cl
  fsom <- UpdateDerivedValues(fsom)
  fsom$info$parameters <- match.call()
  fsom$info$date <- as.character(Sys.time())
  fsom$info$version <- as.character(utils::packageVersion("FlowSOM"))
  if(!silent) message(t[3], "\n")
  return(fsom)
}

#' Print FlowSOM object
#' 
#' @param x FlowSOM object to print information about
#' @param ...  Further arguments, not used
#' @examples 
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- FlowSOM(fileName, compensate = TRUE, transform = TRUE,
#'                       scale = TRUE, colsToUse = c(9, 12, 14:18), nClus = 10)
#' print(flowSOM.res)
#' 
#' @export
print.FlowSOM <- function(x, ...){
  if(!is.null(x$map)){
    cat("FlowSOM model trained on", nrow(x$data), "cells and", 
        length(x$map$colsUsed), "markers, \n using a",
        paste0(x$map$xdim,"x",x$map$ydim), paste0("grid (",NClusters(x)), "clusters) and",
        NMetaclusters(x), "metaclusters.")
    
    cat("\n\nMarkers used: ", paste(x$prettyColnames[x$map$colsUsed], collapse =", "))
  } else {
    cat("FlowSOM model to train on", nrow(x$data), "cells.")
  }
  
  if(!is.null(x$metaclustering)){
    cat("\n\nMetacluster cell count:\n")
    counts <- GetCounts(x)
    print(counts)
  }
  
  if(!is.null(x$outliers)){
    n_outliers <- sum(x$outliers$perCluster$Number_of_outliers)
    n_mad <- round((x$outliers$perCluster[1,"Threshold"] - 
                      x$outliers$perCluster[1,"Median_distance"]) / 
                     x$outliers$perCluster[1,"Median_absolute_deviation"])
    cat("\n", n_outliers, paste0("cells (",
                                 round(100*n_outliers/nrow(x$data),2),
                                 "%)"),
        "are further than", n_mad,"MAD from their cluster center.")
  }
}


#' Aggregate multiple FCS files together
#' 
#' Aggregate multiple FCS files to analyze them simultaneously. 
#' A new FCS file is written, which contains about \code{cTotal} cells,
#' with \code{ceiling(cTotal/nFiles)} cells from each file. Two new columns
#' are added: a column indicating the original file by index, and a noisy 
#' version of this for better plotting opportunities (index plus or minus a 
#' value between 0 and 0.1).
#' 
#' @param fileNames   Character vector containing full paths to the FCS files
#'                    or a flowSet to aggregate
#' @param cTotal      Total number of cells to write to the output file
#' @param channels    Channels/markers to keep in the aggregate. Default NULL 
#'                    takes all channels of the first file.
#' @param writeOutput Whether to write the resulting flowFrame to a file. 
#'                    Default FALSE
#' @param outputFile  Full path to output file. Default "aggregate.fcs"
#' @param keepOrder If TRUE, the random subsample will be ordered in the same
#'                  way as they were originally ordered in the file. Default =
#'                  FALSE.
#' @param silent If FALSE, prints an update every time it starts processing a
#'               new file. Default = FALSE. 
#' @param sampleWithReplacement If TRUE and more cells per file are requested
#'                              than actually present, all cells will be included
#'                              plus additional resampling. Otherwise, at most 
#'                              all cells will be included once. Default = FALSE.
#' @param ...     Additional arguments to pass to read.FCS
#'                  
#' @return This function does not return anything, but will write a file with
#'         about \code{cTotal} cells to \code{outputFile}
#'
#' @seealso \code{\link{ceiling}}
#'
#' @examples
#' # Define filename
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' # This example will sample 2 times 500 cells.
#' ff_new <- AggregateFlowFrames(c(fileName, fileName), 1000)
#' 
#' @importFrom flowCore read.FCS fr_append_cols keyword colnames markernames 
#'             exprs write.FCS sampleNames
#' @importFrom stats rnorm
#' @importFrom methods is
#' 
#' @export
AggregateFlowFrames <- function(fileNames, 
                                cTotal,
                                channels = NULL,
                                writeOutput = FALSE, 
                                outputFile  = "aggregate.fcs",
                                keepOrder = FALSE, 
                                silent = FALSE,
                                sampleWithReplacement = FALSE,
                                ...){
  
  # Compute number of cells per file
  nFiles <- length(fileNames)
  cFile <- ceiling(cTotal/nFiles)
  
  flowFrame <- NULL
  fileMatrix <- NULL
  
  diffNumberChannels <- FALSE
  diffMarkers <- FALSE
  
  for(i in seq_len(nFiles)){
    if(is(fileNames, "flowSet")) {
      file_name <- sampleNames(fileNames)[i]
      f <- fileNames[[i]]
    } else {
      file_name <- fileNames[i]
      if(!silent) {message("Reading ", file_name)}
      f <- flowCore::read.FCS(fileNames[i], ...)
    }
    # Random sampling
    if(sampleWithReplacement & (nrow(f) < cFile)){
      ids <- c(seq_len(nrow(f)),
               sample(seq_len(nrow(f)), cFile-nrow(f), replace = TRUE))
    } else {
      ids <- sample(seq_len(nrow(f)), min(nrow(f), cFile))
    }
   
    
    if(keepOrder) ids <- sort(ids)
    
    
    
    file_ids <- rep(i, min(nrow(f), cFile))
    
    m <- cbind(file_ids,
               file_ids + stats::rnorm(length(file_ids), 0, 0.1),
               ids)
    
    f <- f[ids,] 
    if(!all(channels %in% colnames(f))){
      diffNumberChannels <- TRUE
      channelsToAdd <- channels[ ! channels %in% colnames(f)]
      extracols <- matrix(0, 
                          nrow = flowCore::nrow(f), 
                          ncol = length(channelsToAdd),
                          dimnames = list(NULL, channelsToAdd))
      f <- flowCore::fr_append_cols(f, extracols)
    }
    
    if(is.null(flowFrame)){
      fileMatrix <- m
      if(is.null(channels)){
        channels <- colnames(f)
        flowFrame <- f
      } else {
        channels <- GetChannels(f, channels)
        flowFrame <- f[, channels, drop = FALSE]
      }
      flowCore::keyword(flowFrame)[["$FIL"]] <- basename(outputFile)
      flowCore::keyword(flowFrame)[["FILENAME"]] <- basename(outputFile)
    } else {
      cols_f <- flowCore::colnames(f)
      cols_flowFrame <- flowCore::colnames(flowFrame)
      commonCols <- intersect(cols_f, cols_flowFrame)
      
      if (length(commonCols) == 0) stop("No common channels between files")
      if (!diffNumberChannels && 
          length(cols_flowFrame) != length(commonCols)){
        diffNumberChannels <- TRUE
      }
      
      if (!diffMarkers && 
          any(!flowCore::markernames(f)[commonCols] %in% 
              flowCore::markernames(flowFrame)[commonCols])){
        diffMarkers <- TRUE
      }
      
      flowCore::exprs(flowFrame) <- 
        rbind(flowCore::exprs(flowFrame)[, commonCols, drop = FALSE], 
              flowCore::exprs(f)[, commonCols, drop = FALSE])
      
      fileMatrix <- rbind(fileMatrix, m)
    }
  }
  
  colnames <- c("File", "File_scattered", "Original_ID")
  prev_agg <- length(grep("File[0-9]*$", colnames(flowFrame)))
  if(prev_agg > 0){
    colnames[c(1, 2)] <- paste0(colnames[c(1, 2)], prev_agg + 1)
  }
  prev_ids <- length(grep("Original_ID[0-9]*$", colnames(flowFrame)))
  if(prev_ids > 0){
    colnames[3] <- paste0(colnames[3], prev_ids + 1)
  }
  colnames(fileMatrix) <- colnames
  flowFrame <- flowCore::fr_append_cols(flowFrame, fileMatrix)
  
  if (diffNumberChannels){
    warning("Files do not contain the same number of channels/markers. ",
            "Zeros might have been imputed for missing values.")
  }
  
  if (diffMarkers){ 
    warning("Files do not contain the same markers")
  }
  
  if(writeOutput){
    flowCore::write.FCS(flowFrame, filename = outputFile)
  }
  
  return(flowFrame)
}



#' Process a FlowJo workspace file
#'
#' Reads a FlowJo workspace file using the flowWorkspace library 
#' and returns a list with a matrix containing gating results and a vector with 
#' a label for each cell from a set of specified gates
#'
#' @param files       The FCS files of interest
#' @param wspFile    The FlowJo wsp file to read
#' @param group       The FlowJo group to parse. Default "All Samples".
#' @param cellTypes  Cell types to use for final labeling the cells. Should
#'                    correspond with a subset of the gate names in FlowJo.
#' @param getData    If true, flowFrames are returned as well.
#' @param ...         Extra arguments to pass to CytoML::flowjo_to_gatingset
#'
#' @return This function returns a list, which for every file contains a list
#' in which the first element ("matrix") is a matrix containing filtering 
#' results for each specified gate and the second element ("manual") is a vector
#' which assigns one label to each cell. If only one file is given, only one
#' list is returned instead of a list of lists.
#'
#' @seealso \code{\link{PlotPies}}
#'
#' @examples
#'
#' # Identify the files
#' fcs_file <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' wspFile <- system.file("extdata", "gating.wsp", package = "FlowSOM")
#' 
#' # Specify the cell types of interest for assigning one label per cell
#' cellTypes <- c("B cells",
#'                 "gd T cells", "CD4 T cells", "CD8 T cells",
#'                 "NK cells", "NK T cells")
#'
#' # Parse the FlowJo workspace   
#' gatingResult <- GetFlowJoLabels(fcs_file, wspFile,
#'                                 cellTypes = cellTypes,
#'                                 getData = TRUE)
#'
#' # Check the number of cells assigned to each gate
#' colSums(gatingResult$matrix)
#' 
#' # Build a FlowSOM tree
#' flowSOM.res <- FlowSOM(gatingResult$flowFrame,
#'                        colsToUse = c(9, 12, 14:18),
#'                        nClus = 10,
#'                        seed = 1)
#'    
#'  # Plot pies indicating the percentage of cell types present in the nodes
#'  PlotPies(flowSOM.res,
#'           gatingResult$manual,
#'           backgroundValues = flowSOM.res$metaclustering)
#'
#' @export
GetFlowJoLabels <- function(files,
                            wspFile,
                            group = "All Samples",
                            cellTypes = NULL,
                            getData = FALSE,
                            ...) {
  if (requireNamespace("CytoML", quietly = TRUE) & 
      requireNamespace("flowWorkspace", quietly = TRUE)) {
    ws <- CytoML::open_flowjo_xml(wspFile, sample_names_from = "sampleNode")
    samples_in_ws <- CytoML::fj_ws_get_samples(ws)
    subset <- unlist(lapply(basename(files), function(f) which(samples_in_ws[,"name"] == f)))
    gates <- CytoML::flowjo_to_gatingset(ws, 
                                         name = group,
                                         subset = subset,
                                         ...)
    files_in_wsp <- flowWorkspace::sampleNames(gates)
    counts <- as.numeric(gsub(".*_([0-9]*)$", "\\1", files_in_wsp)) 
    files_in_wsp <- gsub("_[0-9]*$", "", files_in_wsp)
    result <- list()
    for(file in files){
      print(paste0("Processing ", file))
      file_id <- grep(paste0("^\\Q", basename(file), "\\E$"), 
                      files_in_wsp)
      if(length(file_id) == 0) {stop("File ", basename(file), 
                                     " not found. Files available: \n",
                                     paste0(files_in_wsp, "\n"))}
      gate_names <- flowWorkspace::gs_get_pop_paths(gates[[file_id]], 
                                                    path = "auto")
      
      gatingMatrix <- matrix(NA,
                             nrow = counts[file_id],
                             ncol = length(gate_names),
                             dimnames = list(NULL,
                                             gate_names))
      for(gate in gate_names){
        gatingMatrix[, gate] <- 
          flowWorkspace::gh_pop_get_indices(gates[[file_id]], gate)
      }
      
      if(is.null(cellTypes)){
        cellTypes_tmp <- flowWorkspace::gs_get_leaf_nodes(gates[[file_id]],
                                                          path = "auto")
      } else {
        cellTypes_tmp <- cellTypes
      }
      
      manual <- ManualVector(gatingMatrix, cellTypes_tmp)
      
      result[[file]] <- list("matrix" = gatingMatrix,
                             "manual" = manual)
      
      if (getData) {
        result[[file]]$flowFrame <- 
          flowWorkspace::gh_pop_get_data(gates[[file_id]])
      }
    }
    
    if (length(files) == 1){
      result <- result[[1]]
    }
    
    return(result)
  } else {
    message(paste0("The packages \"CytoML\" and \"flowWorkspace\" are necessary", 
                   " for the function GetFlowJoLabels, but are not available."))
  }
}

#' Summarize the gating matrix into one vector, only including the cell types of
#' interest
#'
#' Extract the compensated and transformed data and all gate labels.
#'
#' @param manualMatrix Matrix containing boolean values, indicating for every
#'                      gate (column) whether the cell (row) is part of it or not.
#' @param cellTypes Cell types to use in the summary vector. All others will be
#'                   ignored and cells which do not fall in one of these gates
#'                   will get the label "Unknown". Order is important!
#'
#' @return A factor with one label for every cell
#'
#' @export
ManualVector <- function(manualMatrix, cellTypes){
  
  if(is.list(manualMatrix)){ 
    manualMatrix <- do.call(rbind, manualMatrix) 
  }
  
  manual <- rep("Unlabeled", nrow(manualMatrix))
  for(cellType in cellTypes){
    manual[manualMatrix[, cellType]] <- cellType
  }
  manual <- factor(manual, levels=c("Unlabeled", cellTypes))
  return(manual)
}

#' GetChannels
#' 
#' Get channel names for an array of markers, given a flowFrame or a FlowSOM
#' object. As available in "name". \code{\link{grep}} is used to look for the 
#' markers. Other regex can be added.
#' 
#' @param object  The flowFrame or the FlowSOM object of interest 
#' @param markers Vector with markers or channels of interest. Also accepts the
#'                index of the marker found in the object.
#' @param exact   If TRUE (default), the grep pattern will be extended to
#'                start with ^\\\\Q and end with \\\\E$, so only exact matches 
#'                are possible.
#'                  
#' @return Corresponding channel names
#'
#' @seealso \code{\link{GetMarkers}}
#'
#' @examples
#' 
#'    # Read the flowFrame
#'    fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#'    ff <- flowCore::read.FCS(fileName)
#'    GetChannels(ff, c("FSC-A", "CD3", "FITC-A"))
#'    GetMarkers(ff, c("FSC-A", "CD3", "FITC-A"))
#'
#' @export
GetChannels <- function(object, markers, exact = TRUE) { 
  if (is(object, "flowFrame")) {
    object_channels <- unname(flowCore::parameters(object)@data[["name"]])
    object_markers <- unname(flowCore::parameters(object)@data[["desc"]])
  } else if (is(object, "FlowSOM")) {
    object_channels <- names(object$prettyColnames)
    object_markers <- unname(gsub(" <.*", "", object$prettyColnames))
  } else {
    stop("Object should be of class flowFrame or FlowSOM")
  }
  
  if (is.logical(markers)) markers <- which(markers)
  
  channelnames <- c()
  for (marker in markers){
    if(is.numeric(marker)) {
      iChannel <- marker
    } else {
      if(exact) marker <- paste0("^\\Q", marker, "\\E$")
      iChannel <- grep(marker, object_markers)
    }
    if (length(iChannel) != 0){
      for (i in iChannel){
        channel <- object_channels[i]
        names(channel) <- object_markers[i]
        channelnames <- c(channelnames, channel)
      }
    } else {
      iChannel <- grep(marker, object_channels)
      if (length(iChannel) != 0){
        channel <- object_channels[iChannel]
        names(channel) <- channel
        channelnames <- c(channelnames, channel)
      } else {
        stop(paste("Marker", marker, "could not be found"))
      }
    }
  }
  return(channelnames)
}

#' GetMarkers
#' 
#' Get marker names for an array of channels, given a flowFrame or a FlowSOM 
#' object. As available in "desc". If this is NA, defaults to channel name. 
#' \code{\link{grep}} is used to look for the markers. Other regex can be added.
#' 
#' @param object   The flowFrame or the FlowSOM object of interest 
#' @param channels Vector with markers or channels of interest. Also accepts the
#'                 index of the channel in the object.
#' @param exact   If TRUE (default), the grep pattern will be extended to
#'                start with ^\\\\Q and end with \\\\E$, so only exact matches 
#'                are possible.
#'                                  
#' @return Corresponding marker names
#'
#' @seealso \code{\link{GetChannels}}
#'
#' @examples
#' 
#'    # Read the flowFrame
#'    fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#'    ff <- flowCore::read.FCS(fileName)
#'    GetChannels(ff, c("FSC-A", "CD3", "FITC-A"))
#'    GetMarkers(ff, c("FSC-A", "CD3", "FITC-A"))
#'
#' @export
GetMarkers <- function(object, channels, exact = TRUE) { 
  if (is(object, "flowFrame")) {
    object_channels <- unname(flowCore::parameters(object)@data[["name"]])
    object_markers <- unname(flowCore::parameters(object)@data[["desc"]])
  } else if (is(object, "FlowSOM")) {
    object_channels <- names(object$prettyColnames)
    object_markers <- unname(gsub(" <.*", "", object$prettyColnames))
  } else {
    stop("Object should be of class flowFrame or FlowSOM")
  }
  
  if (is.logical(channels)) channels <- which(channels)
  
  markernames <- c()
  for (channel in channels){
    if (is.numeric(channel)) {
      iMarker <- channel
    } else {
      if (exact) channel <- paste0("^\\Q", channel, "\\E$")
      iMarker <- grep(channel, object_channels)
    }
    if (length(iMarker) != 0){
      for (i in iMarker){
        marker <- object_markers[i]
        if (is.na(marker)) marker <- object_channels[i]
        names(marker) <- object_channels[i]
        markernames <- c(markernames, marker)
      }
    } else {
      iMarker <- grep(channel, object_markers)
      if (length(iMarker) != 0){
        marker <- object_markers[iMarker]
        names(marker) <- marker
        markernames <- c(markernames, marker)
      } else {
        stop(paste("Channel", channel, "could not be found"))
      }
    }
  }
  return(markernames)
}

#' UpdateFlowSOM 
#' 
#' Update old FlowSOM object to a new one and checks if it is a flowSOM object
#' 
#' Determines whether or not the fsom input is of class "FlowSOM" and returns  
#' the FlowSOM object and metaclustering object inside fsom
#' 
#' @param fsom  FlowSOM object, as generated by \code{\link{BuildMST}} or
#'              \code{\link{FlowSOM}}
#'                          
#' @return A FlowSOM object
#' 
#' @seealso \code{\link{PlotFlowSOM}}
#' 
#' @importFrom dplyr group_by summarise_all select
#' @importFrom stats  median
#' 
#' @export
UpdateFlowSOM <- function(fsom){
  if (is(fsom,"list") && !is.null(fsom$FlowSOM)) {
    fsom$FlowSOM$metaclustering <- fsom$metaclustering
    fsom <- fsom$FlowSOM
  }
  if (!is(fsom,"FlowSOM")) {
    stop("fsom should be a FlowSOM object.")
  }
  fsom$prettyColnames <- gsub("\\((.*)\\)", "<\\1>", fsom$prettyColnames)
  if (is.null(fsom$map$pctgs)){
    pctgs <- rep(0, fsom$map$nNodes)
    names(pctgs) <- as.character(seq_len(fsom$map$nNodes))
    pctgs_tmp <- table(fsom$map$mapping[, 1]) / nrow(fsom$map$mapping)
    pctgs[names(pctgs_tmp)] <- pctgs_tmp
    fsom$map$pctgs <- pctgs
  } 
  if (is.null(fsom$map$nMetaclusters)){
    fsom$map$nMetaclusters <- length(levels(fsom$metaclustering))
  }
  if (is.null(fsom$map$metaclusterMFIs) && !is.null(fsom$metaclustering)){
    fsom$map$metaclusterMFIs <- 
      data.frame(fsom$data, 
                 mcl = fsom$metaclustering[fsom$map$mapping[, 1]],
                 check.names = FALSE) %>% 
      dplyr::group_by(.data$mcl, .drop = FALSE) %>% 
      dplyr::summarise_all(stats::median) %>% 
      dplyr::select(-.data$mcl) %>% 
      data.frame(row.names = levels(fsom$metaclustering),
                 check.names = FALSE)
  }
  return(fsom)
}
