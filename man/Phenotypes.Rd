\name{Phenotypes-class}
\Rdversion{1.1}
\docType{class}
\alias{Phenotypes-class}
\alias{Phenotypes}

\title{Class \code{"Phenotypes"}}
\description{
The return data from running \code{\link{flowType}}, containing counts 
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Phenotypes", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{CellFreqs}:}{Numeric vector containing counts of the number of cells belonging to each cell type. }
    \item{\code{MFIs}:}{Matrix of MFIs, with rows for cell types and columns for markers.}
    \item{\code{PhenoCodes}:}{Vector of character strings representing the codes of each cell type (phenotype).}
    \item{\code{PropMarkers}:}{Numeric vector specifying which markers were used for combinatorial gating.}
    \item{\code{MFIMarkers}:}{Numeric vector specifying for which markers MFIs were computed for each cell type.}
    \item{\code{MarkerNames}:}{A character vector of the names of all markers in the flowFrame given}
    \item{\code{Partitions}:}{The first level partitions that each cell in the flowFrame belong to in each channel.}
    \item{\code{MaxPopSize}:}{MaxMarkersPerPop}
    \item{\code{PartitionsPerMarker}:}{Vector of number of partitions used for each marker}
    \item{\code{Thresholds}:}{A list of vectors with the calculated thresholds for each marker (if a clustering algorithm was used) or the thresholds provided by the user.}
  }
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{flowType}}
}
\examples{
showClass("Phenotypes")
}
\keyword{classes}
