\name{encodePhenotype}
\alias{encodePhenotype}
\alias{encodePhenotype-methods}
\alias{encodePhenotype,character,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function \code{encodePhenotype} in Package \pkg{flowType}
}
\description{

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
In \pkg{flowType}, phenotypes themselves are represented by codes (e.g. 012) rather than full strings (CD4+CD8-), in order to save memory when a very large number of phenotypes are considered. \code{encodePhenotype} serves to translate a human-readable string down to flowType's internal coded representation. 
}
\usage{
encodePhenotype(pheno.string, marker.names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno.string}{
character vector containing containing the phenotype string to be encoded 
}
  \item{marker.names}{
vector of character vectors each specifying the name of a channel, in order
%%     ~~Describe \code{marker.names} here~~
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

\value{
Character vector containing the encoded phenotype, with one character per channel.
}

\author{
Kieran O'Neill
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{decodePhenotype}}, \code{\link{flowType}}
}
\examples{
encodePhenotype('CD34++CD3-CD45+++',c('CD34','CD3','CD45','CD19'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }

