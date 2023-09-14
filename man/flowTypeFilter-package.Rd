\name{flowType-package}
\alias{flowType-package}
\docType{package}
\title{
  flowType: Phenotyping Flow Cytometry Assays
}

\description{
flowType uses a simple threshold, Kmeans, flowMeans, or flowClust to partition every channel to a positive and a negative cell population.
These partitions are then combined to generate a set of multi-dimensional phenotypes.
}
\details{
\tabular{ll}{
Package: \tab flowType\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2011-04-27\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
Depends: \tab methods\cr
}
For a given FCS file, the flowType function extracts a the phenotypes
and reports their cell frequencies (number of cells) and mean fluorescence intensity (MFI)s.
}
\author{
  Nima Aghaeepour, Kieran O'Neill, Adrin Jalali
}
\references{
  Please cite the following for the current version of flowType:
  
  O'Neill K, Jalali A, Aghaeepour N, Hoos H, Brinkman RR. Enhanced 
  flowType/RchyOptimyx: a BioConductor pipeline for discovery in
  high-dimensional cytometry data.  Bioinformatics. 
  2014 May 1;30(9):1329-30. doi: 10.1093/bioinformatics/btt770
  
  The original paper and description can be found at:
  
  Nima Aghaeepour, Pratip K. Chattopadhyay, Anuradha Ganesan, Kieran
  O'Neill, Habil Zare, Adrin Jalali, Holger H. Hoos, Mario Roederer,
  and Ryan R. Brinkman. Early Immunologic Correlates of HIV Protection
  can be Identified from Computational Analysis of Complex Multivariate
  T-cell Flow Cytometry Assays. Bioinformatics, 2011.
}


\keyword{ clustering }
\keyword{ classification }
\keyword{ HIV }
\keyword{ FlowCytData }

\examples{
#Load the library
library(flowType)
data(DLBCLExample)
MarkerNames <- c('Time', 'FSC-A','FSC-H','SSC-A','IgG','CD38','CD19','CD3','CD27','CD20', 'NA', 'NA')

#These markers will be analyzed
PropMarkers <- 3:5
MFIMarkers <- PropMarkers
MarkerNames <- c('FS', 'SS','CD3','CD5','CD19')

#Run flowType
Res <- flowType(DLBCLExample, PropMarkers, MFIMarkers, 'kmeans', MarkerNames);

MFIs=Res@MFIs;
Proportions=Res@CellFreqs;
Proportions <- Proportions / max(Proportions)
names(Proportions) <- unlist(lapply(Res@PhenoCodes, 
                      function(x){return(decodePhenotype(
                      x,Res@MarkerNames[PropMarkers],
                      Res@PartitionsPerMarker))}))

#Select the 30 largest phenotypes
index=order(Proportions,decreasing=TRUE)[1:30]
bp=barplot(Proportions[index], axes=FALSE, names.arg=FALSE)
text(bp+0.2, par("usr")[3]+0.02, srt = 90, adj = 0, labels = names(Proportions[index]), xpd = TRUE, cex=0.8)
axis(2);
axis(1, at=bp, labels=FALSE);
title(xlab='Phenotype Names', ylab='Cell Proportion')

#These phenotype can be analyzed using a predictive model (e.g., classification or regression)
}
