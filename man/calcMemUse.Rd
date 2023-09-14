\name{calcMemUse}
\alias{calcMemUse}
\title{Function \code{calcMemUse} in Package \pkg{flowType}}
\description{
Estimates the memory usage in bytes for running flowType with a given set of parameters.
}
\usage{
calcMemUse(NumPops, NumPropMarkers, NumMFIMarkers, NumCells, MaxMarkersPerPop, PartitionsPerChannel = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NumPops}{
  Number of cell types which will be returned. Can be computed using \code{\link{calcNumPops}}
%%     ~~Describe \code{NumPops} here~~
}
  \item{NumPropMarkers}{
  Numer of markers to use for combinatorial gating
%%     ~~Describe \code{NumPropMarkers} here~~
}
  \item{NumMFIMarkers}{
  Number of markers to determine the MFIs of for every cell type
%%     ~~Describe \code{NumMFIMarkers} here~~
}
  \item{NumCells}{
  Number of cells in the flowFrame passed to flowType
%%     ~~Describe \code{NumCells} here~~
}
  \item{MaxMarkersPerPop}{
  Maximum number of markers to use at once in combinatorial gating (ie all cell types over 1:MaxMarkersPerPop will be counted)
%%     ~~Describe \code{MaxMarkersPerPop} here~~
}
  \item{PartitionsPerChannel}{
  Number of partitions per channel. 
%%     ~~Describe \code{PartitionsPerChannel} here~~
}
}
\details{
If you use different numbers of partitions for different channels, try providing the highest number as \code{PartitionsPerChannel}, and expect an over-estimate.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Estimated memory use in bytes.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Kieran O'Neill
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{calcNumPops}}, \code{\link{flowType}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities}
