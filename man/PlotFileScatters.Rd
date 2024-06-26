% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0_FlowSOM.R
\name{PlotFileScatters}
\alias{PlotFileScatters}
\title{PlotFileScatters}
\usage{
PlotFileScatters(
  input,
  fileID = "File",
  channels = NULL,
  yLim = NULL,
  yLabel = "marker",
  quantiles = NULL,
  names = NULL,
  groups = NULL,
  color = NULL,
  legend = FALSE,
  maxPoints = 50000,
  ncol = NULL,
  nrow = NULL,
  width = NULL,
  height = NULL,
  silent = FALSE,
  plotFile = "FileScatters.png",
  ...
)
}
\arguments{
\item{input}{Either a flowSet, a flowFrame with a file ID column (e.g. 
output from the \code{\link{AggregateFlowFrames}} includes
a "File" column) or a vector of paths pointing to FCS files}

\item{fileID}{Name of the file ID column when the input is a flowFrame, 
default to "File" (File ID column in the 
\code{\link{AggregateFlowFrames}} flowFrame output).}

\item{channels}{Vector of channels or markers that need to be plotted, 
if NULL (default), all channels from the input will be 
plotted}

\item{yLim}{Optional vector of a lower and upper limit of the y-axis}

\item{yLabel}{Determines the label of the y-axis. Can be "marker" and\\or
"channel" or abbrevations. Default = "marker".}

\item{quantiles}{If provided (default NULL), a numeric vector with values
between 0 and 1. These quantiles are indicated on the plot}

\item{names}{Optional parameter to provide filenames. If \code{NULL} 
(default), the filenames will be numbers. Duplicated 
filenames will be made unique.}

\item{groups}{Optional parameter to specify groups of files, should have
the same length as the \code{input}. Id \code{NULL} 
(default), all files will be plotted in the same color}

\item{color}{Optional parameter to provide colors. Should have the same
lengths as the number of groups (or 1 if \code{groups} is 
\code{NULL})}

\item{legend}{Logical parameter to specify whether the group levels 
should be displayed. Default is \code{FALSE}}

\item{maxPoints}{Total number of data points that will be plotted per 
channel, default is 50000}

\item{ncol}{Number of columns in the final plot, optional}

\item{nrow}{Number of rows in the final plot, optional}

\item{width}{Width of png file. By default NULL the width parameter is 
estimated based on the input.}

\item{height}{Height of png file. By default NULL the width parameter is 
estimated based on the input.}

\item{silent}{If FALSE, prints an update every time it starts processing 
a new file. Default = FALSE.}

\item{plotFile}{Path to png file, default is "FileScatters.png". If 
\code{NULL}, the output will be a list of ggplots}

\item{...}{Arguments for read.FCS (e.g. truncate_max_range)}
}
\value{
List of ggplot objects if \code{plot} is \code{FALSE}, 
        otherwise \code{filePlot} with plot is created.
}
\description{
Make a scatter plot per channel for all provided files
}
\examples{
# Preprocessing
fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
ff <- flowCore::read.FCS(fileName)
ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
ff <- flowCore::transform(ff,
         flowCore::transformList(colnames(flowCore::keyword(ff)[["SPILL"]]),
                                flowCore::logicleTransform()))

flowCore::write.FCS(ff[1:1000, ], file = "ff_tmp1.fcs")
flowCore::write.FCS(ff[1001:2000, ], file = "ff_tmp2.fcs")
flowCore::write.FCS(ff[2001:3000, ], file = "ff_tmp3.fcs")
 
# Make plot
PlotFileScatters(input = c("ff_tmp1.fcs", "ff_tmp2.fcs", "ff_tmp3.fcs"),
                 channels = c("Pacific Blue-A", 
                              "Alexa Fluor 700-A", 
                              "PE-Cy7-A"), 
                 maxPoints = 1000)

}
