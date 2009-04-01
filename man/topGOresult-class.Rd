\name{topGOresult-class}
\docType{class}

\alias{topGOresult-class}


\alias{initialize,topGOresult-method}

\alias{description,topGOresult-method}
\alias{description<-,topGOresult,ANY-method}

\alias{score}
\alias{score,topGOresult,character-method}
\alias{score,topGOresult,missing-method}
\alias{score<-,topGOresult-method}

\alias{testName}
\alias{testName<-}
\alias{testName,topGOresult-method}
\alias{testName<-,topGOresult-method}

\alias{algorithm}
\alias{algorithm<-}
\alias{algorithm<-,topGOresult-method}
\alias{algorithm,topGOresult-method}

\alias{show,topGOresult-method}
\alias{print,topGOresult-method}

\alias{geneData<-}
\alias{geneData}
\alias{geneData<-,topGOresult-method}
\alias{geneData,topGOresult-method}



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
