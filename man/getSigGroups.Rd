\name{getSigGroups}

\alias{getSigGroups}

\alias{getSigGroups-methods}
\alias{getSigGroups,topGOdata,classicCount-method}
\alias{getSigGroups,topGOdata,classicScore-method}
\alias{getSigGroups,topGOdata,elimCount-method}
\alias{getSigGroups,topGOdata,elimScore-method}
\alias{getSigGroups,topGOdata,weightCount-method}


%- Also NEED an '\alias' for EACH other topic documented here.

\title{Algorithms for scoring GO terms}

\description{
  TODO: This function is use for dispatching each algorithm
}

\usage{
getSigGroups(object, test.stat, ...)
}

%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{test.stat}{ ~~Describe \code{test.stat} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}


\author{Adrian Alexa}

\seealso{
  \code{\link{topGOdata-class}},
  \code{\link{classicCount-class}},
  \code{\link{classicScore-class}}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(object, test.stat, ...) standardGeneric("getSigGroups")
}

\keyword{methods}
