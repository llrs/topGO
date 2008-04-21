\name{buildGOgraph.topology}

\alias{buildGOgraph.topology}

\title{buils GO graph starting from the most specific GO terms}

\description{
  This function is building the GO graph starting from the most specific terms.
  The structure of the GO graph is build recursively.
}

\usage{
buildGOgraph.topology(knownNodes, whichOnto = "BP")
}


\arguments{
  \item{knownNodes}{character vector of GO terms}
  \item{whichOnto}{character string specifying one of the three GO
    ontologies: \code{"BP"}, \code{"MF"}, \code{"CC"}}
}


\value{
  An object of class \code{\link[graph]{graphNEL-class}} is returned.
  The graph is directed (the edges are from leaves to the root)
  and it contains all GO terms specific to \code{whichOnto} ontology.
}

\author{Adrian Alexa}

\seealso{
  \code{\link{topGOdata-class}},
  \code{\link{buildLevels}},
  \code{\link{mapGenes2GOgraph}},
  \code{\link{annFUN}}
}


\keyword{graphs}
