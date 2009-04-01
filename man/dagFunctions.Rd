\name{Determines the levels of a Directed Acyclic Graph (DAG)}

\alias{buildLevels}
\alias{getNoOfLevels}
\alias{getGraphRoot}
\alias{reverseArch}



\title{Utility functions for Directed Acyclic Graphs (DAG)}

\description{

  Determines the levels of a Directed Acyclic Graph (DAG)

  
  TODO: This function take the a directed graph and constructs a named
  vector which contain the level on which a node is. The root has level
  1.

  TODO: Find the root(roots) of the DAG

  TODO: Simple function to invert the direction of edges in an directed graph.
  The returned graph is of class graphNEL.
  It can use either simple matrices or sparse matrices (SparseM library)

}

\usage{
buildLevels(dag, root = NULL, leafs2root = TRUE)
getNoOfLevels(graphLevels)
getGraphRoot(dag, leafs2root = TRUE)
reverseArch(dirGraph, useAlgo = "sparse", useWeights = TRUE)
}


\arguments{
  \item{dag}{ ~~Describe \code{dag} here~~ }
  \item{root}{ ~~Describe \code{root} here~~ }
  \item{leafs2root}{The leafs2root parameter tell if the graph has edges
    directed from the leaves to the root, or vice-versa}
  \item{graphLevels}{ ~~Describe \code{graphLevels} here~~ }  
  
  \item{dirGraph}{The graph to be transformed }
  \item{useAlgo}{"sparse" or "normal"}
  \item{useWeights}{If weights should be used (if useAlgo = 'normal' that the weigths are used anyway)}
  
}


\details{
  .....
}

\value{
  it returns a list containing:
  \item{level2nodes}{Environment where the key is the level number
    with the value being the nodes on that level.}
  \item{nodes2level}{Environment where the key is the node label (the GO
    ID) and the value is the level on which that node lies.}
  \item{noOfLevels}{The number of levels}
  \item{noOfNodes}{The number of nodes}

  An object of class \code{\link[graph]{graphNEL-class}} is returned.
      
}

\author{Adrian Alexa}

\seealso{
  \code{\link{topGOdata-class}},
  \code{\link{inducedGraph}}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}

\keyword{graphs}
