\name{topGOresult-class}
\docType{class}

\alias{topGOresult-class}

\alias{description,topGOresult-method}
\alias{description<-,topGOresult,ANY-method}


\alias{testName,topGOresult-method}
\alias{testClass,topGOresult-method}
\alias{testName<-,topGOresult-method}
\alias{testClass<-,topGOresult-method}

\alias{score,topGOresult-method}
\alias{score<-,topGOresult-method}

\alias{initialize,topGOresult-method}

\alias{score}
\alias{testName}
\alias{testClass}
\alias{testName<-}
\alias{testClass<-}


\title{Class "topGOresult"}

\description{Class instance created by \code{\link{getSigGroups-methods}}}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("topGOresult", description, score, testName, testClass)}.
}

\section{Slots}{
  \describe{
    \item{\code{description}:}{Object of class \code{"character"} ~~ }
    \item{\code{score}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{testName}:}{Object of class \code{"character"} ~~ }
    \item{\code{testClass}:}{Object of class \code{"character"} ~~ }
  }
}

\section{Methods}{
  \describe{
    
    \item{\code{score}:}{ ~~describe this method here }
}}


\author{Adrian Alexa}

\seealso{
  \code{\link{classicScore-class}},
  \code{\link{groupStats-class}},
  \code{\link{getSigGroups-methods}}
}


\keyword{classes}
