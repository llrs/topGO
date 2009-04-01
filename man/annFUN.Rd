\name{annFUN}

\alias{annFUN.db}
\alias{annFUN}
\alias{annFUN.GO2genes}
\alias{annFUN.gene2GO}
\alias{annFUN.file}
\alias{annFUN.org}

\alias{inverseList}
\alias{readMappings}


\title{functions to map gene IDs to GO terms}
\description{
  These functions are used to compile a list of GO terms and their
  mappings to gene identifiers.
}


\usage{
annFUN.db(whichOnto, feasibleGenes = NULL, affyLib)
annFUN.org(whichOnto, feasibleGenes = NULL, mapping, ID = "entrez") 
annFUN(whichOnto, feasibleGenes = NULL, affyLib)
annFUN.gene2GO(whichOnto, feasibleGenes = NULL, gene2GO)
annFUN.GO2genes(whichOnto, feasibleGenes = NULL, GO2genes)
annFUN.file(whichOnto, feasibleGenes = NULL, file, ...)

readMappings(file, sep = "\t", IDsep = ",")
inverseList(l)
}


\arguments{
  \item{whichOnto}{character string specifying one of the three GO
    ontologies: \code{"BP"}, \code{"MF"}, \code{"CC"}}

  \item{feasibleGenes}{character vector containing a subset of gene
    identifiers. Only these genes will be used to annotate GO
    terms. Default value is \code{NULL} which means all gene identifiers
    will be used.}

  \item{affyLib}{character string containing the name of the Affymetrix
    chip.}

  \item{gene2GO}{named list of character vectors. The list names are
    genes identifiers. For each gene the character vector contains the
    GO terms IDs it maps to. Only the most specific annotations are required.}
  
  \item{GO2genes}{named list of character vectors. The list names are
    GO terms IDs. For each GO the character vector contains the
    genes identifiers which are mapped to it. Only the most specific
    annotations are required.}
  
  \item{mapping}{......}
  \item{ID}{.....}
  \item{file}{.....}
  \item{...}{other parameters}

  \item{sep}{.....}
  \item{IDsep}{....}
  \item{l}{......}
}

\details{
  The function \code{annFUN.db} uses the mappings provided
  in the Bioconductor annotation data packages. For example, if the
  Affymetrix hgu133a chip it is used, then the user should set
  \code{affyLib = "hgu133a.db"}.

  The functions \code{annFUN.gene2GO} and \code{annFUN.GO2genes} are
  used when the user provide his own annotations.

  All these function restrict the GO terms to the ones belonging
  to the specified ontology. 
}


\value{
  A named(GO terms IDs) list of character vectors. 
}


\author{Adrian Alexa}

\seealso{
  \code{\link{topGOdata-class}}
}

\examples{

library(hgu133a.db)
set.seed(111)


## generate a gene list and the GO annotations
numGenes <- 50
selGenes <- sample(ls(hgu133aGO), numGenes)
gene2GO <- lapply(mget(selGenes, envir = hgu133aGO), names)
gene2GO[sapply(gene2GO, is.null)] <- NA

## the annotation for the first three genes
gene2GO[1:3]

## inverting the annotations
go2genes <- annFUN.gene2GO(whichOnto = "CC", gene2GO = gene2GO)



## generate a GO list with the genes annotations
numGO <- 30
selGO <- sample(ls(hgu133aGO2PROBE), numGO)
GO2gene <- lapply(mget(selGO, envir = hgu133aGO2PROBE), as.character)

GO2gene[1:3]

## select only the GO terms for a specific ontology
go2gene <- annFUN.GO2genes(whichOnto = "CC", GO2gene = GO2gene)

}

\keyword{misc}
