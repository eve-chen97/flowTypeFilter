\name{getLabels}
\alias{getLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  getLabels: Returns the labels of the cells in a given phenotype.
}
\description{
  Returns the labels of the cells in a given phenotype in a Phenotypes object.
}
\usage{
getLabels(Phenotypes, PhenotypeNumber)
}
\arguments{
  \item{Phenotypes}{
    An object of class \code{Phenotypes} as produced by the
    \code{flowType} function.
  }
  \item{PhenotypeNumber}{
    A numeric or character value representing the phenotypes number of
    name, respectively. 
  }
}

\value{
  \item{\code{Membership Labels}:}{A vector of length of the number of
  events. 1 and 2 represent the cells that are not-included and included
  in the phenotype respectively.}
}

\references{
  Nima Aghaeepour, Pratip K. Chattopadhyay, Anuradha Ganesan, Kieran
  O'Neill, Habil Zare, Adrin Jalali, Holger H. Hoos, Mario Roederer,
  and Ryan R. Brinkman. Early Immunologic Correlates of HIV Protection
  can be Identified from Computational Analysis of Complex Multivariate
  T-cell Flow Cytometry Assays. submitted to Bioinformatics, 2011.
}
\author{
  Nima Aghaeepour
}

\examples{
#See the vigentte
}

\keyword{ clustering }
\keyword{ classification }
\keyword{ HIV }
\keyword{ FlowCytData }
