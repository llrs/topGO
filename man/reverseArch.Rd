\name{reverseArch}
\alias{reverseArch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  TODO: Simple function to invert the direction of edges in an directed graph.
  The returned graph is of class graphNEL.
  It can use either simple matrices or sparse matrices (SparseM library)
}

\usage{
reverseArch(dirGraph, useAlgo = "sparse", useWeights = TRUE)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{dirGraph}{The graph to be transformed }
  \item{useAlgo}{"sparse" or "normal"}
  \item{useWeights}{If weights should be used (if useAlgo = 'normal' that the weigths are used anyway)}
}

\details{
  ~~ If necessary, more details than the description above ~~
}

\value{
    An object of class \code{\link[graph]{graphNEL-class}} is returned.
}

\author{Adrian Alexa}

\seealso{
  \code{\link{buildLevels}},
  \code{\link{mapGenes2GOgraph}},
  \code{\link{inducedGraph}}
}


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}

\keyword{graphs}
