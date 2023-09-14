\name{decodePhenotype}
\alias{decodePhenotype}
\alias{decodePhenotype-methods}
\alias{decodePhenotype,character,character,numeric-method}
\alias{decodePhenotype,character, character,numeric-method}
\title{Method \code{decodePhenotype} in Package \pkg{flowType}}
\description{ Method to decode phenotypes back to a human-readable string.

}
\section{Methods}{
\describe{

\item{\code{signature(pheno.code = "character", marker.names = "character", partitions.per.marker = "numeric")}}{
%%  ~~describe this method here~~
In \pkg{flowType}, phenotypes themselves are represented by codes (e.g. 012) rather than full strings (CD4+CD8-), in order to save memory when a very large number of phenotypes are considered. \code{decodePhenotype} serves to translate the codes back to a human-readable string.
}
}
}

\details{
FlowType's encoding is as follows:

0 -- marker not considered in phenotype (don't care about its value)
1 -- marker is negative (e.g. CD4-)
2 -- marker is positive (e.g. CD4+)
3 -- marker is positive, but brighter than 2 (CD4++)
4 -- marker is even brighter (CD4+++)
etc

Note that this encoding system does not allow for "dim" markers -- dim positivity is denoted by the first level of positivity.

Also note that the encoding is performed from the dimmest to the brightest partition, but the location of thresholds will dictate the interpretation of the code. (e.g. if you only set one threshold, but you place it between the positive and the bright population, then both positive and negative events will be considered negative.)
}

\seealso{
\code{\link{encodePhenotype}}, \code{\link{flowType}}
}

\examples{
decodePhenotype('1034',c('CD34','CD3','CD45','CD19'), 4)
decodePhenotype('20013',c('CD34','CD3','CD45','CD19', 'CD20'), 4)
}
\keyword{ utilities }

\author{
Kieran O'Neill
}