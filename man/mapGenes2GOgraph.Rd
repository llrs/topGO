\name{mapGenes2GOgraph}
\alias{mapGenes2GOgraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  TODO: This function builds for each node a vector containing all the genes/probes
  that can be annotated to that node.
  It starts with the nodes on the lowest level, and then pushes their genes
  to the parents/ancestors
}

\usage{
mapGenes2GOgraph(dag, mostSpecificGOs, nodeLevel = buildLevels(dag, leafs2root = TRUE))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dag}{ ~~Describe \code{dag} here~~ }
  \item{mostSpecificGOs}{ ~~Describe \code{mostSpecificGOs} here~~ }
  \item{nodeLevel}{ ~~Describe \code{nodeLevel} here~~ }
}

\details{
  ~~ If necessary, more details than the description above ~~
}

\value{
  An object of class \code{\link[graph]{graphNEL-class}} is returned.
  The attribute of each node in the graph contains a mapping
  of the genes/probes.
}

\author{Adrian Alexa}

\seealso{
  \code{\link{topGOdata-class}},
  \code{\link{buildLevels}},
  \code{\link{buildGOgraph.topology}},
  \code{\link{annFUN}}
}


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}

\keyword{graphs}
