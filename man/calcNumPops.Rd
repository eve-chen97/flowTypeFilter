\name{calcNumPops}
\alias{calcNumPops}
\title{Function \code{calcNumPops} in Package \pkg{flowType}
}
\description{ Compute the number of populations that will be produced by running flowType with a given set of parameters. This is especially useful for estimating memory but is also used internally to determine the size of return objects to pass down to C++. This may also be useful for determining the cutoff for number of markers to use to make phenotypes (in terms of statistical power for later testing).}

\usage{
calcNumPops(PartitionsPerMarker, MaxMarkersPerPop)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PartitionsPerMarker}{
Integer vector specifying the number of partitions for each marker, in order.
}
  \item{MaxMarkersPerPop}{
Integer speciying the threshold chosen
}
}

\value{
Integer specifying the number of populations the given paramters would produce.
}

\author{
Kieran O'Neill
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{calcMemUse}}, \code{\link{flowType}}
}
\examples{
calcNumPops(c(2,2,3,2,2,4), 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }

