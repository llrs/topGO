#################### topGOdata methods ####################

##signature(.Object = "topGOdata",
##          ontology = "character",
##          allGenes = "ANY",
##          geneSelectionFun = "function",
##          description = "character"),
##          annotationFun = "function",

setMethod("initialize", "topGOdata",
          function(.Object,
                   ## which Ontology to be used
                   ontology,
                   ## a named numeric or factor, the names are the genes ID
                   allGenes, 
                   ## function to select the signif. genes
                   geneSelectionFun = NULL,
                   ## description of the class
                   description = character(0),
                   ## expression matrix         #!
                   expressionMatrix = NULL,     #!
                   ## phenotype information     #!
                   phenotype = NULL,            #!
                   ## annotation function
                   annotationFun,
                   ## additional parameters for the annotationFun
                   ...) {
            
            .Object@description <- description
            .Object@ontology <- ontology

            .Object@expressionMatrix <- expressionMatrix  #!
            .Object@phenotype        <- phenotype         #!
            
            ## some checking
            if(is.null(names(allGenes)))
              stop("allGenes must be a named vector")
            
            if(!is.factor(allGenes) && !is.numeric(allGenes))
              stop("allGenes should be a factor or a numeric vector")

            .Object@allGenes <- names(allGenes)
            
            if(is.factor(allGenes)) {
              if(length(levels(allGenes)) != 2)
                stop("allGenes must be a factor with 2 levels")
              .Object@allScores <- factor(as.character(allGenes))
              .Object@geneSelectionFun <- function(x) {
                return(as.logical(as.integer(levels(x)))[x])
              }
            }
            else {
              .Object@allScores <- as.numeric(allGenes)
              
              ## function to select which genes are significant
              if(is.null(geneSelectionFun))
                warning("No function to select the significant genes provided!")
              .Object@geneSelectionFun <- geneSelectionFun
            }
            
            ## loading some required libraries 
            require('GO') || stop('package GO is required')

            ## this function is returning a list of GO terms from the specified ontology
            ## whith each entry being a vector of genes
            cat("\nBuilding most specific GOs .....")
            mostSpecificGOs <- annotationFun(ontology, .Object@allGenes, ...)
            cat("\t(", length(mostSpecificGOs), "GO terms found. )\n")
            
            ## the the GO graph is build started from the most specific terms
            cat("\nBuild GO DAG topology ..........")
            g <- buildGOgraph.topology(mostSpecificGOs, ontology)
            cat("\t(",  numNodes(g), "GO terms and", numEdges(g), "relations. )\n")
                
            ## probably is good to store the leves but for the moment we don't 
            .nodeLevel <- buildLevels(g, leafs2root = TRUE)
            
            ## annotate the nodes in the GO graph with genes
            cat("\nAnnotating nodes ...............")
            g <- mapGenes2GOgraph(g, mostSpecificGOs, nodeLevel = .nodeLevel) ## leafs2root
            
            ## select the feasible genes
            gRoot <- getGraphRoot(g)
            feasibleGenes <- ls(nodeData(g, n = gRoot, attr = "genes")[[gRoot]])
            cat("\t(", length(feasibleGenes), "genes annotated to the GO terms. )\n")

            .Object@feasible <- .Object@allGenes %in% feasibleGenes
            .Object@graph <- g
            
            ## some clean-up
            ##detach(pos = which(search() == "package:GO"))
            
            .Object
          })



### temp  accessor methods ###
if(!isGeneric("expressionMatrix"))
  setGeneric("expressionMatrix", function(object) standardGeneric("expressionMatrix"))
                               
setMethod("expressionMatrix", "topGOdata", function(object) object@expressionMatrix)

if(!isGeneric("phenotype"))
  setGeneric("phenotype", function(object) standardGeneric("phenotype"))
                               
setMethod("phenotype", "topGOdata", function(object) object@phenotype)
######################################################################




#################### the accessor functions ####################

if(!isGeneric("description"))
  setGeneric("description", function(object) standardGeneric("description"))
                               
setMethod("description", "topGOdata", function(object) object@description)

if(!isGeneric("ontology"))
  setGeneric("ontology", function(object) standardGeneric("ontology"))
                               
setMethod("ontology", "topGOdata", function(object) object@ontology)

if(!isGeneric("allGenes"))
  setGeneric("allGenes", function(object) standardGeneric("allGenes"))
                               
setMethod("allGenes", "topGOdata", function(object) object@allGenes)

##if(!isGeneric("allScores"))
##  setGeneric("allScores", function(object) standardGeneric("allScores"))
                               
##setMethod("allScores", "topGOdata", function(object) object@allScores)

if(!isGeneric("feasible"))
  setGeneric("feasible", function(object) standardGeneric("feasible"))
                               
setMethod("feasible", "topGOdata", function(object) object@feasible)

if(!isGeneric("graph"))
  setGeneric("graph", function(object) standardGeneric("graph"))
                               
setMethod("graph", "topGOdata", function(object) object@graph)

if(!isGeneric("geneSelectionFun"))
  setGeneric("geneSelectionFun", function(object) standardGeneric("geneSelectionFun"))
                               
setMethod("geneSelectionFun", "topGOdata", function(object) object@geneSelectionFun)


#################### the replacement functions ####################

if(!isGeneric("description<-"))
  setGeneric("description<-", function(object, value) standardGeneric("description<-"))
                               
setMethod("description<-", "topGOdata", function(object, value) {object@description <- value; object})

if(!isGeneric("ontology<-"))
  setGeneric("ontology<-", function(object, value) standardGeneric("ontology<-"))
                               
setMethod("ontology<-", "topGOdata", function(object, value) {object@ontology <- value; object})

## not a good idea to modify the list of allGenes, see updateGenes()
#if(!isGeneric("allGenes<-"))
#  setGeneric("allGenes<-", function(object, value) standardGeneric("allGenes<-"))
#                               
#setMethod("allGenes<-", "topGOdata", function(object, value) {object@allGenes <- value; object})


if(!isGeneric("feasible<-"))
  setGeneric("feasible<-", function(object, value) standardGeneric("feasible<-"))
                               
setMethod("feasible<-", "topGOdata", function(object, value) {object@feasible <- value; object})


if(!isGeneric("geneSelectionFun<-"))
  setGeneric("geneSelectionFun<-", function(object, value) standardGeneric("geneSelectionFun<-"))
                               
setMethod("geneSelectionFun<-", "topGOdata",
          function(object, value) {object@geneSelectionFun <- value; object})


if(!isGeneric("graph<-"))
  setGeneric("graph<-", function(object, value) standardGeneric("graph<-"))
                               
setMethod("graph<-", "topGOdata", function(object, value) {object@graph <- value; object})



#################### methods to update/modify the objects ####################

## this function is used for updating the genes
## for the moment it allows to change the gene score and to restrict the gene names
## to the set of feasible genes
if(!isGeneric("updateGenes"))
  setGeneric("updateGenes", function(object, geneList, geneSelFun) standardGeneric("updateGenes"))

## for the case in which each gene has a score
setMethod("updateGenes",
          signature(object = "topGOdata", geneList = "numeric", geneSelFun = "function"),
          function(object, geneList, geneSelFun) {
            
            fGenes <- genes(object)
            
            if(is.null(names(geneList)))
              stop("geneList must be a named vector")
            
            if(!all(fGenes %in% names(geneList)))
              stop("Please provide a feasible geneList")
               
            object@allGenes <- names(geneList)

            object@allScores <- as.numeric(geneList)
            object@geneSelectionFun <- geneSelFun

            ## set up the feasible genes index vector
            object@feasible <- object@allGenes %in% fGenes
            
            ## TODO ........
            ## we need to update the graph 
            ## object@graph <- updateGraph(object@graph, object@feasible)

            object
          })


setMethod("updateGenes",
          signature(object = "topGOdata", geneList = "factor", geneSelFun = "missing"),
          function(object, geneList) {

            fGenes <- genes(object)
            
            if(is.null(names(geneList)))
              stop("geneList must be a named vector")
            
            if(!all(genes(object) %in% names(geneList)))
              stop("Please provide a feasible geneList")
                        
            object@allGenes <- names(geneList)
            
            if(length(levels(geneList)) != 2)
              stop("geneList must be a factor with 2 levels")
            object@allScores <- factor(as.character(geneList))
            object@geneSelectionFun <- function(x) {
              return(as.logical(as.integer(levels(x)))[x])
            }

            ## set up the feasible genes index vector
            object@feasible <- object@allGenes %in% fGenes
            
            ## TODO ........
            ## we need to update the graph 
            ## object@graph <- updateGraph(object@graph, object@feasible)

            object
          })



#################### methods to obtain genes lists ####################

## return the genes that are used in the analysis (the one that can be
## mapped to the specific Ontology)
if(!isGeneric("genes"))
  setGeneric("genes", function(object) standardGeneric("genes"))

setMethod("genes", "topGOdata", function(object) object@allGenes[object@feasible])


if(!isGeneric("numGenes"))
  setGeneric("numGenes", function(object) standardGeneric("numGenes"))

setMethod("numGenes", "topGOdata", function(object) sum(object@feasible))


## 
if(!isGeneric("geneScore"))
  setGeneric("geneScore", function(object, whichGenes, ...) standardGeneric("geneScore"))

setMethod("geneScore",
          signature(object = "topGOdata", whichGenes = "missing"),
          function(object, use.names = FALSE) {
            if(is.factor(object@allGenes))
              retList <- (as.numeric(levels(object@allScores))[object@allScores])[object@feasible]
            else
              retList <- as.numeric(object@allScores)[object@feasible]
            
            if(use.names)
              names(retList) <- genes(object)

            return(retList)
          })

setMethod("geneScore",
          signature(object = "topGOdata", whichGenes = "character"),
          function(object, whichGenes, use.names = TRUE) {
            index <- object@feasible & (object@allGenes %in% whichGenes)
            if(is.factor(object@allGenes))
              retList <- (as.numeric(levels(object@allScores))[object@allScores])[index]
            else
              retList <- as.numeric(object@allScores)[index]
            
            names(retList) <- object@allGenes[index]
            retList <- retList[whichGenes[whichGenes %in% genes(object)]]
            
            if(!use.names) 
              names(retList) <- NULL
            
            return(retList)
          })


if(!isGeneric("sigGenes"))
  setGeneric("sigGenes", function(object) standardGeneric("sigGenes"))

setMethod("sigGenes", "topGOdata",
          function(object) {
            
            ##if(is.null(object@geneSelectionFun))
            ##  return(NULL)
              
            ## select the significant genes and the feasible ones
            sGenesIndex <- object@geneSelectionFun(object@allScores) & object@feasible
            return(object@allGenes[sGenesIndex])
          })

if(!isGeneric("numSigGenes"))
  setGeneric("numSigGenes", function(object) standardGeneric("numSigGenes"))

setMethod("numSigGenes", "topGOdata",
          function(object) {
            return(sum(object@geneSelectionFun(object@allScores) & object@feasible))
          })


############## methods to access the information on GO terms ##############

if(!isGeneric("usedGO"))
  setGeneric("usedGO", function(object) standardGeneric("usedGO"))

setMethod("usedGO", "topGOdata", function(object) nodes(graph(object)))

if(!isGeneric("attrInTerm"))
  setGeneric("attrInTerm", function(object, attr, whichGO) standardGeneric("attrInTerm"))

setMethod("attrInTerm", 
          signature(object = "topGOdata", attr = "character", whichGO = "character"),
          function(object, attr, whichGO) {
            return(.getFromNode(graph(object), attr, whichGO))
          })

setMethod("attrInTerm", 
          signature(object = "topGOdata", attr = "character", whichGO = "missing"),
          function(object, attr, whichGO) {
            return(.getFromNode(graph(object), attr, nodes(graph(object))))
          })

## function that return for each GO term specified in the whichGO parameter the vector
## of annotated genes. If the whichGO is missing than the vector of annotated genes is
## returned for each GO term 
if(!isGeneric("genesInTerm"))
  setGeneric("genesInTerm", function(object, whichGO) standardGeneric("genesInTerm"))

setMethod("genesInTerm", 
          signature(object = "topGOdata", whichGO = "character"),
          function(object, whichGO) {
            return(.genesInNode(graph(object), whichGO))
          })


setMethod("genesInTerm", 
          signature(object = "topGOdata", whichGO = "missing"),
          function(object) {
            return(.genesInNode(graph(object), nodes(graph(object))))
          })

## similar like above but returns the scores in each GO
if(!isGeneric("scoresInTerm"))
  setGeneric("scoresInTerm", function(object, whichGO, ...) standardGeneric("scoresInTerm"))

setMethod("scoresInTerm", 
          signature(object = "topGOdata", whichGO = "character"),
          function(object, whichGO, use.names = FALSE) {
            l <- lapply(.genesInNode(graph(object), whichGO),
                        function(x) geneScore(object, x, use.names = use.names))
            return(l)
          })


setMethod("scoresInTerm", 
          signature(object = "topGOdata", whichGO = "missing"),
          function(object, use.names = FALSE) {
            return(scoreInNode(object, nodes(graph(object)), use.names = use.names))
          })



if(!isGeneric("countGenesInTerm"))
  setGeneric("countGenesInTerm", function(object, whichGO) standardGeneric("countGenesInTerm"))

setMethod("countGenesInTerm", 
          signature(object = "topGOdata", whichGO = "character"),
          function(object, whichGO) {
            return(.countsInNode(graph(object), whichGO))
          })

setMethod("countGenesInTerm", 
          signature(object = "topGOdata", whichGO = "missing"),
          function(object) {
            return(.countsInNode(graph(object), nodes(graph(object))))
          })
            
## function that return for each GO term specified in the whichGO parameter
## a dataframe withe the following informations:
##   annotated - the number of annotated genes 
##   significant  - the number of annotated significant genes
##   expected - how many genes are expected in a random case
if(!isGeneric("termStat"))
  setGeneric("termStat", function(object, whichGO) standardGeneric("termStat"))

setMethod("termStat", 
          signature(object = "topGOdata", whichGO = "character"),
          function(object, whichGO) {
            
            x <- .genesInNode(graph(object), whichGO)
            
            anno <- sapply(x, length)
            
            sGenes <- sigGenes(object)
            sig <- sapply(x, function(e) length(intersect(e, sGenes)))

            expect <- numSigGenes(object) / length(genes(object))
            expect <- round(anno * expect, 2)
            
            return(data.frame(Annotated = anno,
                              Significant = sig,
                              Expected = expect,
                              row.names = names(x)))
          })

setMethod("termStat", 
          signature(object = "topGOdata", whichGO = "missing"),
          function(object) {
            termStat(object, nodes(graph(object)))
          })



## write information in the graph nodes
if(!isGeneric("updateTerm<-"))
  setGeneric("updateTerm<-", function(object, attr, value)  standardGeneric("updateTerm<-"))

setMethod("updateTerm<-", 
          signature(object = "topGOdata", attr = "character", value = "ANY"),
          function(object, attr, value) {
            object@graph <- .writeToNodes(graph(object), attr, value)
            object
          })


##############################  printing   ##############################
## make use of both "print" and "show" generics


setMethod("print", "topGOdata", function(x, ...) .printTopGOdata(x))
setMethod("show", "topGOdata", function(object) .printTopGOdata(x = object))
          
.printTopGOdata <- function(x) {
  cat("\n------------------------- topGOdata object -------------------------\n")
  cat("\n Description:\n")
  cat("   - ", x@description, "\n")
  
  cat("\n Ontology:\n")
  cat("   - ", x@ontology, "\n")

  ## all genes from the array
  cat("\n", length(x@allGenes) ,"available genes (all genes from the array):\n")
  sym <- x@allGenes[1:min(length(x@allGenes), 5)]
  cat("   - symbol: ", sym, " ...\n")
  if(is.numeric(x@allScores)) {
    score <- x@allScores[1:min(length(x@allGenes), 5)]
    score <- apply(cbind(score, nchar(sym)), 1, function(x) format(x[1], digits = max(x[2] - 2, 1)))
    cat("   - score : ", score, " ...\n")
  }
  cat("   -", sum(x@geneSelectionFun(x@allScores)), " significant genes. \n")
  
  ## feasible genes
  cat("\n", numGenes(x) ,"feasible genes (genes that can be used in the analysis):\n")
  sym <- genes(x)[1:min(numGenes(x), 5)]
  cat("   - symbol: ", sym, " ...\n")
  if(is.numeric(x@allScores)) {
    score <- geneScore(x)[1:min(numGenes(x), 5)]
    score <- apply(cbind(score, nchar(sym)), 1, function(x) format(x[1], digits = max(x[2] - 2, 1)))
    cat("   - score : ", score, " ...\n")
  }
  cat("   -", numSigGenes(x), " significant genes. \n")

  cat("\n GO graph:\n")
  cat("   - a graph with", edgemode(x@graph), "edges\n")
  cat("   - number of nodes =", numNodes(x@graph), "\n")
  cat("   - number of edges =", numEdges(x@graph), "\n")
  
  ##cat("\n Signif. genes sellection:\n\n")
  ##print(x@geneSelectionFun)
  cat("\n------------------------- topGOdata object -------------------------\n\n")
}



## TODO ....
## function to return the name of GO terms which have the no. of annotated
## genes in some interval


######################################################################

######################## topGOresult methods ########################

setMethod("initialize", "topGOresult",
          function(.Object, description = character(),
                   score, testName, testClass) {
            .Object@description <- description
            .Object@score <- score
            .Object@testName <- testName
            .Object@testClass <- testClass

            .Object
          })


#################### the accessor functions ####################

setMethod("description", "topGOresult", function(object) object@description)

if(!isGeneric("score"))
  setGeneric("score", function(object) standardGeneric("score"))
                               
setMethod("score", "topGOresult", function(object) object@score)


if(!isGeneric("testName"))
  setGeneric("testName", function(object) standardGeneric("testName"))
                               
setMethod("testName", "topGOresult", function(object) object@testName)


if(!isGeneric("testClass"))
  setGeneric("testClass", function(object) standardGeneric("testClass"))
                               
setMethod("testClass", "topGOresult", function(object) object@testClass)

#################### the replacement functions ####################

setMethod("description<-", "topGOresult",
          function(object, value) {object@description <- value; object})

if(!isGeneric("score<-"))
  setGeneric("score<-", function(object, value) standardGeneric("score<-"))
                               
setMethod("score<-", "topGOresult", function(object, value) {object@score <- value; object})


if(!isGeneric("testName<-"))
  setGeneric("testName<-", function(object, value) standardGeneric("testName<-"))
                               
setMethod("testName<-", "topGOresult", function(object, value) {object@testName <- value; object})


if(!isGeneric("testClass<-"))
  setGeneric("testClass<-", function(object, value) standardGeneric("testClass<-"))
                               
setMethod("testClass<-", "topGOresult", function(object, value) {object@testClass <- value; object})




######################################################################

######################## groupInfo methods ########################

setMethod("initialize", "groupStats",
          function(.Object, testStatistic, name, allMembers, groupMembers, testStatPar) {
            .Object@name <- name
            .Object@allMembers <- allMembers
            .Object@members <- groupMembers
            .Object@testStatistic <- testStatistic
            .Object@testStatPar <- testStatPar
            
            .Object
          })


#################### the accessor functions ####################

if(!isGeneric("Name"))
  setGeneric("Name", function(object) standardGeneric("Name"))
                               
setMethod("Name", "groupStats", function(object) object@name)

if(!isGeneric("allMembers"))
  setGeneric("allMembers", function(object) standardGeneric("allMembers"))
                               
setMethod("allMembers", "groupStats", function(object) object@allMembers)

if(!isGeneric("members"))
  setGeneric("members", function(object) standardGeneric("members"))
                               
setMethod("members", "groupStats", function(object) object@members)

if(!isGeneric("testStatistic"))
  setGeneric("testStatistic", function(object) standardGeneric("testStatistic"))
                               
setMethod("testStatistic", "groupStats", function(object) object@testStatistic)

if(!isGeneric("testStatPar"))
  setGeneric("testStatPar", function(object) standardGeneric("testStatPar"))
                               
setMethod("testStatPar", "groupStats", function(object) object@testStatPar)

#################### the replacement functions ####################

if(!isGeneric("Name<-"))
  setGeneric("Name<-", function(object, value) standardGeneric("Name<-"))
                               
setMethod("Name<-", "groupStats", function(object, value) {object@Name <- value; object})

if(!isGeneric("allMembers<-"))
  setGeneric("allMembers<-", function(object, value) standardGeneric("allMembers<-"))
                               
setMethod("allMembers<-", "groupStats", function(object, value) {object@allMembers <- value; object})

if(!isGeneric("members<-"))
  setGeneric("members<-", function(object, value) standardGeneric("members<-"))
                               
setMethod("members<-", "groupStats", function(object, value) {object@members <- value; object})


#################### other functions ####################

if(!isGeneric("numMembers"))
  setGeneric("numMembers", function(object) standardGeneric("numMembers"))
                               
setMethod("numMembers", "groupStats", function(object) length(object@members))

if(!isGeneric("numAllMembers"))
  setGeneric("numAllMembers", function(object) standardGeneric("numAllMembers"))
                               
setMethod("numAllMembers", "groupStats", function(object) length(object@allMembers))


## MAIN function -- it should return the "p-value" of the Test Satatistic
if(!isGeneric("runTest"))
  setGeneric("runTest", function(object) standardGeneric("runTest"))
                               
setMethod("runTest", "groupStats",
          function(object) object@testStatistic(object))


## function to update/(build a new object)
if(!isGeneric("updateGroup"))
  setGeneric("updateGroup", function(object, name, members, ...) standardGeneric("updateGroup"))

setMethod("updateGroup",
          signature(object = "groupStats", name = "character", members = "character"),
          function(object, name, members) {
            object@name <- name
            object@members <- members

            return(object)
          })


#################### classicCount methods ####################

#################### constructor ####################
setMethod("initialize", "classicCount",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   sigMembers = character(),
                   ...) {
            .Object <- callNextMethod(.Object, testStatistic, name, allMembers,
                                      groupMembers, testStatPar = list(...))
            .Object@significant <- which(allMembers %in% sigMembers)

            .Object
          })


if(!isGeneric("sigMembers<-"))
  setGeneric("sigMembers<-", function(object, value) standardGeneric("sigMembers<-"))
                               
setMethod("sigMembers<-", "classicCount",
          function(object, value) {
            object@significant <- which(object@allMembers %in% value)
            object
          })

if(!isGeneric("sigAllMembers")) 
  setGeneric("sigAllMembers", function(object) standardGeneric("sigAllMembers"))
                               
setMethod("sigAllMembers", "classicCount",
          function(object) object@allMembers[object@significant])

if(!isGeneric("numSigAll")) 
  setGeneric("numSigAll", function(object) standardGeneric("numSigAll"))
                               
setMethod("numSigAll", "classicCount",
          function(object) length(object@significant))

if(!isGeneric("sigMembers")) 
  setGeneric("sigMembers", function(object) standardGeneric("sigMembers"))
                               
setMethod("sigMembers", "classicCount",
          function(object) intersect(sigAllMembers(object), members(object)))

if(!isGeneric("numSigMembers")) 
  setGeneric("numSigMembers", function(object) standardGeneric("numSigMembers"))
                               
setMethod("numSigMembers", "classicCount",
          function(object) sum(members(object) %in% sigAllMembers(object)))          

if(!isGeneric("contTable")) 
  setGeneric("contTable", function(object) standardGeneric("contTable"))

setMethod("contTable", "classicCount",
          function(object) {

            numHits <- numSigMembers(object)
            numSig <- numSigAll(object)

            sig <- c(numHits, numSig - numHits)
            not.sig <- c(numMembers(object) - numHits,
                         numAllMembers(object) - numMembers(object) - numSig + numHits)
            
            contMat <- cbind(sig = sig, notSig = not.sig)
            row.names(contMat) <- c("anno", "notAnno")

            ##if(!is.null(.LOG.FILE)) {
            ## cat("\n\n (Cont mat)   ", file = .LOG.FILE, append = TRUE)
            ## write.table(contMat, file = .LOG.FILE, quote = FALSE, sep = "\t",
            ##            append = TRUE)## col.names = colnames(contMat), row.names = TRUE)
            ##}
            
            return(contMat)
          })



#################### classicScore methods ####################

setMethod("initialize", "classicScore",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   score = numeric(),
                   alternative = "less",
                   ...) {

            .alternative <- switch(alternative,
                                   greater = TRUE,
                                   less = FALSE,
                                   stop("alternative should be greater or less"))
            
            ## first we order the members aording to the score
            index <- order(score, decreasing = .alternative)
            if(length(allMembers) != length(score))
              warning("score length don't match.")
            
            .Object <- callNextMethod(.Object, testStatistic, name, allMembers[index],
                                      groupMembers, testStatPar = list(...))
            
            .Object@score <- as.numeric(score)[index]
            .Object@.alternative <- .alternative
            .Object
          })

if(!isGeneric("alternative"))
  setGeneric("alternative", function(object) standardGeneric("alternative"))

setMethod("alternative", "classicScore", function(object) object@.alternative)

## methods to get the score 
if(!isGeneric("allScore"))
  setGeneric("allScore", function(object, use.names) standardGeneric("allScore"))
                               
setMethod("allScore",
          signature(object = "classicScore", use.names = "missing"),
          function(object) object@score)

setMethod("allScore",
          signature(object = "classicScore", use.names = "logical"),
          function(object, use.names = FALSE) {
            if(use.names) {
              ret.val <- object@score
              names(ret.val) <- allMembers(object)
              return(ret.val)
            }

            object@score
          })


if(!isGeneric("membersScore"))
  setGeneric("membersScore", function(object) standardGeneric("membersScore"))
                               
setMethod("membersScore", "classicScore",
          function(object) {
            index <- allMembers(object) %in% members(object)
            ss <- object@score[index]
            names(ss) <- allMembers(object)[index]

            ## just to be safe of the order
            return(ss[members(object)])
          })


## methods to assign the score
## the value should be a named vector, the names should be allMembers
if(!isGeneric("score<-"))
  setGeneric("score<-", function(object, value) standardGeneric("score<-"))
                               
setMethod("score<-", "classicScore",
          function(object, value) {
            ## some checking
            if(length(names(value)) == 0)
              stop("no names associated with the score")
            if(!all(names(value) %in% object@allMembers))
              warning("The new score names do not match with allMembers.")
            if(!all(object@members %in% names(value)))
              stop("You need to build a new object!")
            
            value <- sort(value, object@.alternative)
            object@allMembers <- names(value)
            object@score <- as.numeric(value)

            return(object)
          })
          
## rank the members of the group according to their score
if(!isGeneric("rankMembers"))
  setGeneric("rankMembers", function(object) standardGeneric("rankMembers"))
                               
setMethod("rankMembers", "classicScore",
          function(object) {
            return(which(object@allMembers %in% object@members))
          })



######################################################################





#################### elimCount methods ####################

#################### constructor ####################
setMethod("initialize", "elimCount",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   sigMembers = character(),
                   elim = character(),
                   cutOff = 0.01,
                   ...) {
            .Object <- callNextMethod(.Object, testStatistic, name,
                                      allMembers, groupMembers,
                                      sigMembers, testStatPar = list(...))

            .Object@elim <- which(.Object@members %in% elim)
            .Object@cutOff <- cutOff

            .Object
          })

if(!isGeneric("elim<-"))
  setGeneric("elim<-", function(object, value) standardGeneric("elim<-"))
                               
setMethod("elim<-", "elimCount",
          function(object, value) {
            object@elim <- which(object@members %in% value)
            object
          })

if(!isGeneric("cutOff<-"))
  setGeneric("cutOff<-", function(object, value) standardGeneric("cutOff<-"))
                               
setMethod("cutOff<-", "elimCount",
          function(object, value)  {object@cutOff <- value; object})


if(!isGeneric("cutOff")) 
  setGeneric("cutOff", function(object) standardGeneric("cutOff"))
                               
setMethod("cutOff", "elimCount", function(object) object@cutOff)

if(!isGeneric("elim")) 
  setGeneric("elim", function(object) standardGeneric("elim"))
                               
setMethod("elim", "elimCount",
          function(object) object@members[object@elim])

## account for the eliminated members
#TODO:

setMethod("numMembers", "elimCount",
          function(object) length(object@members) - length(object@elim))

setMethod("numAllMembers", "elimCount",
          function(object) length(object@allMembers) - length(object@elim))

setMethod("sigAllMembers", "elimCount",
          function(object) setdiff(object@allMembers[object@significant], elim(object)))
          
setMethod("numSigAll", "elimCount",
          function(object) length(sigAllMembers(object)))

setMethod("sigMembers", "elimCount",
          function(object) intersect(sigAllMembers(object), members(object)))

setMethod("numSigMembers", "elimCount",
          function(object) length(sigMembers(object)))





#################### elimScore methods ####################

#################### constructor ####################
setMethod("initialize", "elimScore",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   score = numeric(),
                   alternative = "less",
                   elim = character(),
                   cutOff = 0.01,
                   ...) {
            .Object <- callNextMethod(.Object, testStatistic, name,
                                      allMembers, groupMembers, score,
                                      alternative, testStatPar = list(...))

            .Object@elim <- which(.Object@members %in% elim)
            .Object@cutOff <- cutOff

            .Object
          })

setMethod("elim<-", "elimScore",
          function(object, value) {
            object@elim <- which(object@members %in% value)
            object
          })

setMethod("cutOff<-", "elimScore",
          function(object, value)  {object@cutOff <- value; object})


setMethod("cutOff", "elimScore", function(object) object@cutOff)
                               
setMethod("elim", "elimScore", function(object) object@members[object@elim])

setMethod("members", "elimScore",
          function(object) {
            if(length(object@elim) == 0)
              return(object@members)

            return(object@members[-object@elim])
          })

setMethod("allMembers", "elimScore",
          function(object) object@allMembers[!(object@allMembers %in% elim(object))])

setMethod("numMembers", "elimScore",
          function(object) length(object@members) - length(object@elim))

setMethod("numAllMembers", "elimScore",
          function(object) length(object@allMembers) - length(object@elim))

setMethod("allScore",
          signature(object = "elimScore", use.names = "missing"),
          function(object) {
            ret.val <- object@score
            index <- !(object@allMembers %in% elim(object))

            return(ret.val[index])
          })

setMethod("allScore",
          signature(object = "elimScore", use.names = "logical"),
          function(object, use.names = FALSE) {
            ret.val <- object@score
            index <- !(object@allMembers %in% elim(object))
            
            if(use.names)
              names(ret.val) <- object@allMembers
            
            return(ret.val[index])
          })

setMethod("membersScore", "elimScore",
          function(object) {
            ## since members(object) has the eliminated members removed
            index <- object@allMembers %in% members(object)

            ss <- object@score[index]
            names(ss) <- object@allMembers[index]

            ## just to be safe of the order
            return(ss[members(object)])
          })


setMethod("rankMembers", "elimScore",
          function(object) {
            ## NEED to use methods allMembers and members since they accont for removed members
            return(which(allMembers(object) %in% members(object)))
          })





#################### weightCount methods ####################

#################### constructor ####################
setMethod("initialize", "weightCount",
          function(.Object,
                   testStatistic,
                   name = character(),
                   allMembers = character(),
                   groupMembers = character(),
                   sigMembers = character(),
                   weights = numeric(),
                   sigRatio = "ratio",
                   penalise = "more",
                   ...) {

            .Object <- callNextMethod(.Object, testStatistic,
                                      paste(name, sigRatio, sep = " : "),
                                      allMembers, groupMembers,
                                      sigMembers, testStatPar = list(...))

            if(length(weights) == 0)
              .Object@weights <- numeric()
            else 
              if(is.null(names(weights)) &&
                 length(intersect(names(weights), .Object@members)) != length(.Object@members))
                stop("The weight vector must be a named vector.")
            
            .Object@weights <- as.numeric(weights[.Object@members]) ## to remove the names

            ## probably the easiest way to store which sigRatio function is used
            .Object@sigRatio <- switch(sigRatio,
                                       ratio = .sigRatio.ratio,
                                       log = .sigRatio.log,
                                       "01" = .sigRatio.01,
                                       stop("Please give a valid sigRation name"))

            if(penalise == "more") 
              ## need to think about this more
              .Object@penalise <- function(a, b) return(1 / (abs(log(a * b) / 2) + 1))
            else
              .Object@penalise <- function(a, b) return(1)
          
            .Object@roundFun <- floor
            
            .Object
          })

if(!isGeneric("penalise")) 
  setGeneric("penalise", function(object, a, b) standardGeneric("penalise"))

setMethod("penalise",
          signature(object = "weightCount", a = "numeric", b = "numeric"),
          function(object, a, b) object@penalise(a, b))


setMethod("updateGroup",
          signature(object = "weightCount", name = "character", members = "character"),
          function(object, name, members, weights) {
            object <- callNextMethod(object, name, members)
            if(!missing(weights) && length(members) > 0) {
              if(length(intersect(names(weights), members)) != length(members))
                stop("weights and members vectors don't agree.")

              object@weights <- as.numeric(weights[members])
            }
            
            return(object)
          })

if(!isGeneric("Weights")) 
  setGeneric("Weights", function(object, use.names) standardGeneric("Weights"))

setMethod("Weights",
          signature(object = "weightCount", use.names = "missing"), 
          function(object) object@weights)

setMethod("Weights",
          signature(object = "weightCount", use.names = "logical"),
          function(object, use.names = FALSE) {
            if(use.names) {
              ret.val <- object@weights
              names(ret.val) <- members(object)
              return(ret.val)
            }
            
            object@weights
          })


if(!isGeneric("Weights<-"))
  setGeneric("Weights<-", function(object, value) standardGeneric("Weights<-"))
                               
setMethod("Weights<-", "weightCount",
          function(object, value) {
            object@weights <- as.numeric(value[object@members])
            object
          })


#if(!isGeneric("sigRatio")) 
#  setGeneric("sigRatio", function(object) standardGeneric("sigRatio"))
#                               
#setMethod("sigRatio", "weightCount", function(object) object@sigRatio)

if(!isGeneric("sigRatio<-"))
  setGeneric("sigRatio<-", function(object, value) standardGeneric("sigRatio<-"))
                               
setMethod("sigRatio<-", "weightCount",
          function(object, value)  {
            object@sigRatio <- switch(value,
                                      log = .sigRatio.log,
                                      ratio = .sigRatio.ratio,
                                      "01" = .sigRatio.01,
                                      stop("Please give a valid sigRation name"))
            object
          })


if(!isGeneric("getSigRatio")) 
  setGeneric("getSigRatio", function(object, a, b) standardGeneric("getSigRatio"))
                               
setMethod("getSigRatio", "weightCount",
          function(object, a, b) object@sigRatio(a, b))


## account for the weights of the members
#TODO:

setMethod("numMembers", "weightCount",
          function(object) object@roundFun(sum(object@weights)))

setMethod("numAllMembers", "weightCount",
          function(object) {
            a <- sum(object@weights) + (length(object@allMembers) - length(object@members))
            return(object@roundFun(a))
          })

setMethod("numSigAll", "weightCount",
          function(object) {
            ## which group members are sig.
            index <- members(object) %in% sigAllMembers(object)
            ## how many sig. members in total
            num.ones <- length(object@significant) - sum(index)

            return(object@roundFun(sum(object@weights[index]) + num.ones))
          })

setMethod("numSigMembers", "weightCount",
          function(object) {
            ## sum the weights of the group members which are sig.
            return(object@roundFun(sum(object@weights[object@members %in% sigAllMembers(object)])))
          })
