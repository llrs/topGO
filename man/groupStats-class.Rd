\name{groupStats-class}
\docType{class}
\alias{groupStats-class}

\alias{allMembers<-}
\alias{allMembers}
\alias{members<-}
\alias{members}
\alias{Name<-}
\alias{Name}
\alias{numAllMembers}
\alias{numMembers}
\alias{runTest}
\alias{testStatistic}
\alias{testStatPar}
\alias{updateGroup}



\alias{allMembers<-,groupStats-method}
\alias{allMembers,groupStats-method}
\alias{initialize,groupStats-method}
\alias{members<-,groupStats-method}
\alias{members,groupStats-method}
\alias{Name<-,groupStats-method}
\alias{Name,groupStats-method}
\alias{numAllMembers,groupStats-method}
\alias{numMembers,groupStats-method}
\alias{runTest,groupStats-method}
\alias{testStatistic,groupStats-method}
\alias{testStatPar,groupStats-method}
\alias{updateGroup,groupStats,character,character-method}


\title{ Class "groupStats"}
\description{A virtual class containing basic group (GO term) data: gene
  names, genes scores, etc...
}

\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
	 \describe{
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
    \item{\code{members}:}{Object of class \code{"character"} ~~ }
    \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{allMembers<-}{\code{signature(object = "groupStats")}: ... }
    \item{allMembers}{\code{signature(object = "groupStats")}: ... }
    \item{initialize}{\code{signature(.Object = "groupStats")}: ... }
    \item{members<-}{\code{signature(object = "groupStats")}: ... }
    \item{members}{\code{signature(object = "groupStats")}: ... }
    \item{Name<-}{\code{signature(object = "groupStats")}: ... }
    \item{Name}{\code{signature(object = "groupStats")}: ... }
    \item{numAllMembers}{\code{signature(object = "groupStats")}: ... }
    \item{numMembers}{\code{signature(object = "groupStats")}: ... }
    \item{runTest}{\code{signature(object = "groupStats")}: ... }
    \item{testStatistic}{\code{signature(object = "groupStats")}: ... }
	 }
}
\author{Adrian Alexa}

\seealso{
  \code{\link{classicCount-class}},
  \code{\link{getSigGroups-methods}}
}

\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{classes}
