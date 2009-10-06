\name{groupGOTerms}
\alias{groupGOTerms}

\alias{GOCCTerm}
\alias{GOBPTerm}
\alias{GOMFTerm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  TODO: Function that split GOTERM in different ontologies.
  Every new environment contain only the terms from one
  of the ontologies 'BP', 'CC', 'MF'
}

\usage{
groupGOTerms(where)
}

\arguments{
  \item{where}{The the environment where you wantto bind the results}
}

%- maybe also 'usage' for other objects documented here.
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
  \code{\link[GO.db]{GOTERM}}
}

\examples{
groupGOTerms()
}

\keyword{misc}
