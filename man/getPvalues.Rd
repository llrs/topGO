\name{getPvalues}

\alias{getPvalues}


\title{Function to compute p-values of a t-test for a gene expression matrix.}
\description{

}

\usage{
   getPvalues(edata, classlabel, test = "t", alternative = c("greater", "two.sided", "less")[1],
   genesID = NULL, correction = c("none", "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD",
   "BH", "BY")[8]) 
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{edata}{Gene expression matrix.}
  \item{classlabel}{The phenotype of the data}
  \item{test}{Which test statistic to use}
  \item{alternative}{The alternative of the test statistic}
  \item{genesID}{if a subset of genes is provided}
  \item{correction}{Multiple testing correction procedure}
}

\details{
  ~~ If necessary, more details than the description above ~~
}

\value{
  An named vector of p-values is returned.
}

\author{Adrian Alexa}

\seealso{
}


\examples{

library(ALL)
data(ALL)

## discriminate B-cell from T-cell
classLabel <- as.integer(sapply(ALL$BT, function(x) return(substr(x, 1, 1) == 'T')))

## Differentially expressed genes
geneList <- getPvalues(exprs(ALL), classlabel = classLabel,
                       alternative = "greater", correction = "BY")

hist(geneList, 50)
}

\keyword{graphs}
