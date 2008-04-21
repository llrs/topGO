\name{elimExpr-class}
\docType{class}

\alias{elimExpr-class}
\alias{removeExpr-class}


\alias{cutOff<-,elimExpr-method}
\alias{cutOff,elimExpr-method}
\alias{initialize,elimExpr-method}


\alias{allMembers,removeExpr-method}
\alias{elim<-,removeExpr-method}
\alias{elim,removeExpr-method}
\alias{initialize,removeExpr-method}
\alias{members,removeExpr-method}
\alias{numAllMembers,removeExpr-method}
\alias{numMembers,removeExpr-method}


\title{Class "elimExpr" ~~~ }
\description{	 ~~ A concise (1-5 lines) description of what the class is.  ~~}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("elimExpr", testStatistic, name, groupMembers, exprDat, pType, elim, cutOff, ...)}.
	 ~~ describe objects here ~~ 
}
\section{Slots}{
	 \describe{
    \item{\code{cutOff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{elim}:}{Object of class \code{"integer"} ~~ }
    \item{\code{eData}:}{Object of class \code{"environment"} ~~ }
    \item{\code{pType}:}{Object of class \code{"factor"} ~~ }
    \item{\code{name}:}{Object of class \code{"character"} ~~ }
    \item{\code{allMembers}:}{Object of class \code{"character"} ~~ }
    \item{\code{members}:}{Object of class \code{"character"} ~~ }
    \item{\code{testStatistic}:}{Object of class \code{"function"} ~~ }
    \item{\code{testStatPar}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{removeExpr}"}, directly.
Class \code{"\linkS4class{classicExpr}"}, by class "removeExpr", distance 2.
Class \code{"\linkS4class{groupStats}"}, by class "removeExpr", distance 3.
}
\section{Methods}{
  \describe{
    \item{cutOff<-}{\code{signature(object = "elimExpr")}: ... }
    \item{cutOff}{\code{signature(object = "elimExpr")}: ... }
    \item{getSigGroups}{\code{signature(object = "topGOdata", test.stat = "elimExpr")}: ... }
    \item{initialize}{\code{signature(.Object = "elimExpr")}: ... }
	 }
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
showClass("elimExpr")
}
\keyword{classes}
