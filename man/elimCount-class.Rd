\name{elimCount-class}
\docType{class}

\alias{elimCount-class}
\alias{removeCount-class}

\alias{cutOff}
\alias{cutOff,elimCount-method}
\alias{cutOff<-}
\alias{cutOff<-,elimCount-method}

\alias{elim}
\alias{elim,elimCount-method}
\alias{elim,removeCount-method}
\alias{elim<-}
\alias{elim<-,elimCount-method}
\alias{elim<-,removeCount-method}

\alias{initialize,elimCount-method}
\alias{initialize,removeCount-method}

\alias{contTable,elimCount-method}
\alias{numSigAll,elimCount-method}
\alias{numSigMembers,elimCount-method}
\alias{numMembers,elimCount-method}
\alias{numAllMembers,elimCount-method}
\alias{sigAllMembers,elimCount-method}
\alias{sigMembers<-,elimCount-method}
\alias{sigMembers,elimCount-method}


\alias{numAllMembers,removeCount-method}
\alias{numMembers,removeCount-method}
\alias{numSigAll,removeCount-method}
\alias{numSigMembers,removeCount-method}
\alias{sigAllMembers,removeCount-method}
\alias{sigMembers,removeCount-method}

\alias{GOFisherTest,elimCount-method}


  
\title{Classes "elimCount" and "removeCount"}
\description{	 ~~ A concise (1-5 lines) description of what the class is.  ~~}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("elimCount", testStatistic, name, allMembers, groupMembers, sigMembers, elim, cutOff, ...)}.
	 ~~ describe objects here ~~ 
}
\section{Slots}{
	 \describe{
    \item{\code{elim}:}{Object of class \code{"integer"} ~~ }
    \item{\code{cutOff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{significant}:}{Object of class \code{"integer"} ~~ }
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
    \item{\code{members}:}{Object of class \code{"character"} ~~ }
    \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
    \item{\code{testStatPar}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{classicCount}"}, directly.
Class \code{"\linkS4class{groupStats}"}, by class "classicCount", distance 2.
}
\section{Methods}{
No methods defined with class "elimCount" in the signature.
}

\details{
  TODO: Some datails here.....
  
}

\author{Adrian Alexa}

\seealso{
  \code{\link{classicScore-class}},
  \code{\link{groupStats-class}},
  \code{\link{getSigGroups-methods}}
}

\examples{
##---- Should be DIRECTLY executable !! ----
}

\keyword{classes}
