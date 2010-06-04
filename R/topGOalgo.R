######################################################################
## file containing all the algorithms for scoring GO terms.


## High-level function for runing the GO algorithms


.algoComp <- rbind(c(1, 0, 1, 1, 1, 0),
                   c(1, 0, 1, 1, 1, 0),
                   c(1, 0, 0, 0, 0, 0),
                   c(1, 0, 1, 1, 1, 0),
                   c(1, 0, 0, 0, 0, 0))
rownames(.algoComp) <- c("classic", "elim", "weight", "weight01", "parentchild")
colnames(.algoComp) <- c("fisher", "z", "ks", "t", "globaltest", "category")

.testNames <- c("GOFisherTest" , "GOKSTest", "GOtTest", "GOglobalTest")
names(.testNames) <- c("fisher", "ks", "t", "globaltest")

.algoClass <- c("classic", "elim", "weight", "weight01", "parentchild")
names(.algoClass) <- c("classic", "elim", "weight", "weight01", "parentchild")


## functions to extract the information from the .algoComp
whichAlgorithms <- function() {
  rownames(.algoComp)
}

whichTests <- function() {
  colnames(.algoComp)[colSums(.algoComp) > 0]
}

## function runTest(topGOdata, "elim", "KS", ...)
## ... are parameters for the test statistic

if(!isGeneric("runTest")) 
  setGeneric("runTest", function(object, algorithm, statistic, ...) standardGeneric("runTest"))

setMethod("runTest",
          signature(object = "topGOdata", algorithm = "character", statistic = "character"),
          function(object, algorithm, statistic, ...) { ## ... parameters for the test statistic

            statistic <- tolower(statistic)
            algorithm <- tolower(algorithm)
            ## we check if the algorithm support the given test statistic
            if(.algoComp[algorithm, statistic] == 0)
              stop("The algorithm doesn't support the given test statistic")

            algorithm <- .algoClass[algorithm]
            ## strong asumtions! 
            if(algorithm == "parentchild")
              testType <- "pC"
            else {
              testType <- as.character(findMethodSignatures(.testNames[statistic]))
              testType <- sub("classic", algorithm, testType) 
            }
            
            if(!isClass(testType))
              stop("Error in defining the test statistic Class")
            
            test.stat <- new(testType, testStatistic = get(.testNames[statistic]), name = statistic, ...)

            return(getSigGroups(object, test.stat))
          })

##resultKS <- runTest(GOdata, algorithm = "classic", statistic = "KS")

## dispatch method for weight01 (or topGO) algorithm - the default algorithm
setMethod("runTest",
          signature(object = "topGOdata", algorithm = "missing", statistic = "character"),
          function(object, statistic, ...) { ## ... parameters for the test statistic
            return(runTest(object, "weight01", statistic, ...))
          })



######################################################################
######################################################################
## This function is use for dispatching each algorithm
## probably more checking could be done in the end!

if(!isGeneric("getSigGroups")) 
  setGeneric("getSigGroups", function(object, test.stat, ...) standardGeneric("getSigGroups"))


########################## CLASSIC algorithm ##########################
## dispatch method for classic algorithm
setMethod("getSigGroups",
          signature(object = "topGOdata", test.stat = "classicCount"),
          function(object, test.stat) { ## ... parameters for each algorithm
            
            ## update the test.stat object
            allMembers(test.stat) <- genes(object)
            sigMembers(test.stat) <- sigGenes(object)

            ## first take aside nodes that don't have any sig members
            x <- termStat(object)
            index <- x$Significant == 0
            not.sigTerms <- rownames(x)[index]
            .sigTerms <- rownames(x)[!index]

            cat("\n\t\t\t -- Classic Algorithm -- \n")
            cat(paste("\n\t\t the algorithm is scoring ", length(.sigTerms), " nontrivial nodes\n", sep =""))
            cat("\t\t parameters: \n")
            cat("\t\t\t test statistic: ", Name(test.stat), "\n")

            ## check if there is at least one nontrivial GO term !!!
            if(length(.sigTerms) > 0) {
              ## apply the algorithm
              algoRes <- .sigGroups.classic(genesInTerm(object, .sigTerms), test.stat)
            } else {
              algoRes <- numeric()
              warning("No enrichment can pe performed - there are no feasible GO terms!")
            }
            
            aux <- rep(1, length(not.sigTerms))
            names(aux) <- not.sigTerms
            algoRes <- c(algoRes, aux)
  
            ## STORE THE RESULTS
            .whichAlgorithm <- "classic"
            attr(.whichAlgorithm, "testClass") <- as.character(class(test.stat))
            return(new("topGOresult",
                       description = paste(description(object), "\nOntology:", ontology(object), sep = " "),
                       score = algoRes, testName = Name(test.stat),
                       algorithm = .whichAlgorithm,
                       geneData = c(.getGeneData(object), SigTerms = length(.sigTerms))))
          })

setMethod("getSigGroups",
          signature(object = "topGOdata", test.stat = "classicScore"),
          function(object, test.stat) {

            ## update the test.stat object if there is the case
            if(length(allScore(test.stat)) == 0) {
              allMembers(test.stat) <- genes(object)
              score(test.stat) <- geneScore(object, use.names = TRUE)
            }
            
            GOlist <- genesInTerm(object)
            cat("\n\t\t\t -- Classic Algorithm -- \n")
            cat(paste("\n\t\t the algorithm is scoring ", length(GOlist), " nontrivial nodes\n", sep =""))
            cat("\t\t parameters: \n")
            cat("\t\t\t test statistic: ", Name(test.stat), "\n")
            cat("\t\t\t score order: ", ifelse(scoreOrder(test.stat), "decreasing", "increasing"), "\n")
            
            algoRes <- .sigGroups.classic(GOlist, test.stat)
            
            ## STORE THE RESULTS
            .whichAlgorithm <- "classic"
            attr(.whichAlgorithm, "testClass") <- as.character(class(test.stat))
            return(new("topGOresult",
                       description = paste(description(object), "\nOntology:", ontology(object), sep = " "),
                       score = algoRes, testName = Name(test.stat),
                       algorithm = .whichAlgorithm,
                       geneData = c(.getGeneData(object), SigTerms = length(GOlist))))
          })


setMethod("getSigGroups",
          signature(object = "topGOdata", test.stat = "classicExpr"),
          function(object, test.stat) {

            if(length(allMembers(test.stat)) == 0)
              stop("In function getSigGroups: no expression data found")

            GOlist <- genesInTerm(object)
            cat("\n\t\t\t -- Classic Algorithm -- \n")
            cat(paste("\n\t\t the algorithm is scoring ", length(GOlist), " nontrivial nodes\n", sep =""))
            cat("\t\t parameters: \n")
            cat("\t\t\t test statistic: ", Name(test.stat), "\n")
            
            algoRes <- .sigGroups.classic(GOlist, test.stat)

            .whichAlgorithm <- "classic"
            attr(.whichAlgorithm, "testClass") <- as.character(class(test.stat))
            return(new("topGOresult",
                       description = paste(description(object), "\nOntology:", ontology(object), sep = " "),
                       score = algoRes, testName = Name(test.stat),
                       algorithm = .whichAlgorithm,
                       geneData = c(.getGeneData(object), SigTerms = length(GOlist))))
          })




########################## ELIM algorithm ##########################
## dispatch method for elim algorithm
## to set the test statistic
## test.stat <- new("elimCount", testStatistic = GOFisherTest, cutOff = 0.05,  name = "Fisher test")
setMethod("getSigGroups",
          signature(object = "topGOdata", test.stat = "elimCount"),
          function(object, test.stat, ...) { ## ... parameters for each algorithm
            
            ## update the test.stat object
            allMembers(test.stat) <- genes(object)
            sigMembers(test.stat) <- sigGenes(object)
            
            ## first take aside nodes that don't have any sig members
            x <- termStat(object)
            index <- x$Significant == 0
            not.sigTerms <- rownames(x)[index]
            .sigTerms <- rownames(x)[!index]

            cat("\n\t\t\t -- Elim Algorithm -- \n")
            cat(paste("\n\t\t the algorithm is scoring ", length(.sigTerms), " nontrivial nodes\n", sep =""))
            cat("\t\t parameters: \n")
            cat("\t\t\t test statistic: ", Name(test.stat), "\n")
            cat("\t\t\t cutOff: ", cutOff(test.stat), "\n")

            ## check if there is at least one nontrivial GO term !!!
            if(length(.sigTerms) > 0) {
              ## restrict the graph
              g <- subGraph(.sigTerms, graph(object))
              ## apply the algorithm
              algoRes <- .sigGroups.elim(g, test.stat)
            } else {
              algoRes <- numeric()
              warning("No enrichment can pe performed - there are no feasible GO terms!")
            }
              
            aux <- rep(1, length(not.sigTerms))
            names(aux) <- not.sigTerms
            algoRes <- c(algoRes, aux)
  
            .whichAlgorithm <- "elim"
            attr(.whichAlgorithm, "testClass") <- as.character(class(test.stat))
            return(new("topGOresult",
                       description = paste(description(object), "\nOntology:", ontology(object), sep = " "),
                       score = algoRes,
                       testName = paste(Name(test.stat), cutOff(test.stat), sep = " : "), 
                       algorithm = .whichAlgorithm,
                       geneData = c(.getGeneData(object), SigTerms = length(.sigTerms))))
          })


setMethod("getSigGroups",
          signature(object = "topGOdata", test.stat = "elimScore"),
          function(object, test.stat, ...) { ## ... parameters for each algorithm

            ## update the test.stat object if there is the case
            if(length(allScore(test.stat)) == 0) {
              allMembers(test.stat) <- genes(object)
              score(test.stat) <- geneScore(object, use.names = TRUE)
            }

            GOlist <- genesInTerm(object)
            cat("\n\t\t\t -- Elim Algorithm -- \n")
            cat(paste("\n\t\t the algorithm is scoring ", length(GOlist), " nontrivial nodes\n", sep =""))
            cat("\t\t parameters: \n")
            cat("\t\t\t test statistic: ", Name(test.stat), "\n")
            cat("\t\t\t cutOff: ", cutOff(test.stat), "\n")
            cat("\t\t\t score order: ", ifelse(scoreOrder(test.stat), "decreasing", "increasing"), "\n")
            
            ## apply the algorithm
            algoRes <- .sigGroups.elim(graph(object), test.stat)
            
            .whichAlgorithm <- "elim"
            attr(.whichAlgorithm, "testClass") <- as.character(class(test.stat))
            return(new("topGOresult",
                       description = paste(description(object), "\nOntology:", ontology(object), sep = " "),
                       score = algoRes,
                       testName = paste(Name(test.stat), cutOff(test.stat), sep = " : "), 
                       algorithm = .whichAlgorithm,
                       geneData = c(.getGeneData(object), SigTerms = length(GOlist))))

          })


setMethod("getSigGroups",
          signature(object = "topGOdata", test.stat = "elimExpr"),
          function(object, test.stat, ...) { ## ... parameters for each algorithm
            
            ## update the test.stat object if there is the case
            if(length(allMembers(test.stat)) == 0) 
              stop("No expression data found")
            
            GOlist <- genesInTerm(object)
            cat("\n\t\t\t -- Elim Algorithm -- \n")
            cat(paste("\n\t\t the algorithm is scoring ", length(GOlist), " nontrivial nodes\n", sep =""))
            cat("\t\t parameters: \n")
            cat("\t\t\t test statistic: ", Name(test.stat), "\n")
            cat("\t\t\t cutOff: ", cutOff(test.stat), "\n")
            
            ## apply the algorithm
            algoRes <- .sigGroups.elim(graph(object), test.stat)
            
            .whichAlgorithm <- "elim"
            attr(.whichAlgorithm, "testClass") <- as.character(class(test.stat))
            return(new("topGOresult",
                       description = paste(description(object), "\nOntology:", ontology(object), sep = " "),
                       score = algoRes,
                       testName = paste(Name(test.stat), cutOff(test.stat), sep = " : "), 
                       algorithm = .whichAlgorithm,
                       geneData = c(.getGeneData(object), SigTerms = length(GOlist))))
          })




########################## weight01(topGO) algorithm ##########################
## dispatch method for weight01 algorithm
## to set the test statistic
## test.stat <- new("weight01Count", testStatistic = GOFisherTest, name = "Fisher test")
setMethod("getSigGroups",
          signature(object = "topGOdata", test.stat = "weight01Count"),
          function(object, test.stat, ...) { ## ... parameters for each algorithm
            
            ## update the test.stat object
            allMembers(test.stat) <- genes(object)
            sigMembers(test.stat) <- sigGenes(object)
            
            ## first take aside nodes that don't have any sig members
            x <- termStat(object)
            index <- x$Significant == 0
            not.sigTerms <- rownames(x)[index]
            .sigTerms <- rownames(x)[!index]

            cat("\n\t\t\t -- Weight01 Algorithm -- \n")
            cat(paste("\n\t\t the algorithm is scoring ", length(.sigTerms), " nontrivial nodes\n", sep =""))
            cat("\t\t parameters: \n")
            cat("\t\t\t test statistic: ", Name(test.stat), "\n")

            ## check if there is at least one nontrivial GO term !!!
            if(length(.sigTerms) > 0) {
              ## restrict the graph
              g <- subGraph(.sigTerms, graph(object))
              ## apply the algorithm
              algoRes <- .sigGroups.weight01(g, test.stat)
            } else {
              algoRes <- numeric()
              warning("No enrichment can pe performed - there are no feasible GO terms!")
            }
              
            aux <- rep(1, length(not.sigTerms))
            names(aux) <- not.sigTerms
            algoRes <- c(algoRes, aux)
  
            .whichAlgorithm <- "weight01"
            attr(.whichAlgorithm, "testClass") <- as.character(class(test.stat))
            return(new("topGOresult",
                       description = paste(description(object), "\nOntology:", ontology(object), sep = " "),
                       score = algoRes, testName = Name(test.stat),
                       algorithm = .whichAlgorithm,
                       geneData = c(.getGeneData(object), SigTerms = length(.sigTerms))))

          })


setMethod("getSigGroups",
          signature(object = "topGOdata", test.stat = "weight01Score"),
          function(object, test.stat, ...) { ## ... parameters for each algorithm

            ## update the test.stat object if there is the case
            if(length(allScore(test.stat)) == 0) {
              allMembers(test.stat) <- genes(object)
              score(test.stat) <- geneScore(object, use.names = TRUE)
            }

            GOlist <- genesInTerm(object)
            cat("\n\t\t\t -- Weight01 Algorithm -- \n")
            cat(paste("\n\t\t the algorithm is scoring ", length(GOlist), " nontrivial nodes\n", sep =""))
            cat("\t\t parameters: \n")
            cat("\t\t\t test statistic: ", Name(test.stat), "\n")
            cat("\t\t\t score order: ", ifelse(scoreOrder(test.stat), "decreasing", "increasing"), "\n")
            
            ## apply the algorithm
            algoRes <- .sigGroups.weight01(graph(object), test.stat)
            
            .whichAlgorithm <- "weight01"
            attr(.whichAlgorithm, "testClass") <- as.character(class(test.stat))
            return(new("topGOresult",
                       description = paste(description(object), "\nOntology:", ontology(object), sep = " "),
                       score = algoRes, testName = Name(test.stat),
                       algorithm = .whichAlgorithm,
                       geneData = c(.getGeneData(object), SigTerms = length(GOlist))))
          })


setMethod("getSigGroups",
          signature(object = "topGOdata", test.stat = "weight01Expr"),
          function(object, test.stat, ...) { ## ... parameters for each algorithm
            
            ## update the test.stat object if there is the case
            if(length(allMembers(test.stat)) == 0) 
              stop("No expression data found")
            
            GOlist <- genesInTerm(object)
            cat("\n\t\t\t -- Weight01 Algorithm -- \n")
            cat(paste("\n\t\t the algorithm is scoring ", length(GOlist), " nontrivial nodes\n", sep =""))
            cat("\t\t parameters: \n")
            cat("\t\t\t test statistic: ", Name(test.stat), "\n")
            
            ## apply the algorithm
            algoRes <- .sigGroups.weight01(graph(object), test.stat)
            
            .whichAlgorithm <- "weight01"
            attr(.whichAlgorithm, "testClass") <- as.character(class(test.stat))
            return(new("topGOresult",
                       description = paste(description(object), "\nOntology:", ontology(object), sep = " "),
                       score = algoRes, testName = Name(test.stat),
                       algorithm = .whichAlgorithm,
                       geneData = c(.getGeneData(object), SigTerms = length(GOlist))))
          })


########################## WEIGHT algorithm ##########################
## dispatch method for the weight algorithm
## to set the test statistic
## test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", ...)
setMethod("getSigGroups",
          signature(object = "topGOdata", test.stat = "weightCount"),
          function(object, test.stat, ...) { ## ... parameters for each algorithm
            
            ## update the test.stat object
            allMembers(test.stat) <- genes(object)
            sigMembers(test.stat) <- sigGenes(object)
            
            ## first take aside nodes that don't have any sig members
            x <- termStat(object)
            index <- x$Significant == 0
            not.sigTerms <- rownames(x)[index]
            .sigTerms <- rownames(x)[!index]
            
            cat("\n\t\t\t -- Weight Algorithm -- \n")
            cat(paste("\n\t\t The algorithm is scoring ", length(.sigTerms), " nontrivial nodes\n", sep =""))
            cat("\t\t parameters: \n")
            cat("\t\t\t test statistic: ", Name(test.stat), "\n")
            
            ## check if there is at least one nontrivial GO term !!!
            if(length(.sigTerms) > 0) {
              ## restrict the graph
              g <- subGraph(.sigTerms, graph(object))
              ## apply the algorithm
              algoRes <- .sigGroups.weight(g, test.stat)
            } else {
              algoRes <- numeric()
              warning("No enrichment can pe performed - there are no feasible GO terms!")
            }
            
            aux <- rep(1, length(not.sigTerms))
            names(aux) <- not.sigTerms
            algoRes <- c(algoRes, aux)

            .whichAlgorithm <- "weight"
            attr(.whichAlgorithm, "testClass") <- as.character(class(test.stat))
            return(new("topGOresult",
                       description = paste(description(object), "\nOntology:", ontology(object), sep = " "),
                       score = algoRes, testName = Name(test.stat),
                       algorithm = .whichAlgorithm,
                       geneData = c(.getGeneData(object), SigTerms = length(.sigTerms))))
          })




########################## Parent-Child algorithm ##########################
## dispatch method for parent-child algorithm
##
## to set the test statistic
## test.stat <- new("pC", testStatistic = GOFisherTest, joinFun = "intersect", name = "Fisher test")
setMethod("getSigGroups",
          signature(object = "topGOdata", test.stat = "pC"),
          function(object, test.stat, ...) { ## ... parameters for each algorithm
            
            ## first take aside nodes that don't have any sig members
            x <- termStat(object)
            index <- x$Significant == 0
            not.sigTerms <- rownames(x)[index]
            .sigTerms <- rownames(x)[!index]

            cat("\n\t\t\t -- Parent-Child Algorithm -- \n")
            cat(paste("\n\t\t the algorithm is scoring ", length(.sigTerms), " nontrivial nodes\n", sep =""))
            cat("\t\t parameters: \n")
            cat("\t\t\t test statistic: ", Name(test.stat), "\n")

            ## check if there is at least one nontrivial GO term !!!
            if(length(.sigTerms) > 0) {
              ## restrict the graph
              g <- subGraph(.sigTerms, graph(object))
              ## apply the algorithm
              algoRes <- .sigGroups.parentChild(g, test.stat, sigGenes(object))
            } else {
              algoRes <- numeric()
              warning("No enrichment can pe performed - there are no feasible GO terms!")
            }
                        
            aux <- rep(1, length(not.sigTerms))
            names(aux) <- not.sigTerms
            algoRes <- c(algoRes, aux)
  
            .whichAlgorithm <- "parentchild"
            attr(.whichAlgorithm, "testClass") <- as.character(class(test.stat))
            return(new("topGOresult",
                       description = paste(description(object), "\nOntology:", ontology(object), sep = " "),
                       score = algoRes, testName = Name(test.stat),
                       algorithm = .whichAlgorithm,
                       geneData = c(.getGeneData(object), SigTerms = length(.sigTerms))))
          })


setMethod("getSigGroups",
          signature(object = "topGOdata", test.stat = "parentChild"),
          function(object, test.stat, ...) { ## ... parameters for each algorithm
            
            ## first take aside nodes that don't have any sig members
            x <- termStat(object)
            index <- x$Significant == 0
            not.sigTerms <- rownames(x)[index]
            .sigTerms <- rownames(x)[!index]

            cat("\n\t\t\t -- Parent-Child Algorithm -- \n")
            cat(paste("\n\t\t the algorithm is scoring ", length(.sigTerms), " nontrivial nodes\n", sep =""))
            cat("\t\t parameters: \n")
            cat("\t\t\t test statistic: ", Name(test.stat), "\n")
            cat("\t\t\t join function: ", test.stat@joinFun, "\n")

            ## check if there is at least one nontrivial GO term !!!
            if(length(.sigTerms) > 0) {
              ## restrict the graph
              g <- subGraph(.sigTerms, graph(object))
              ## apply the algorithm
              algoRes <- .sigGroups.parentChild(g, test.stat, sigGenes(object))
            } else {
              algoRes <- numeric()
              warning("No enrichment can pe performed - there are no feasible GO terms!")
            }
            
            aux <- rep(1, length(not.sigTerms))
            names(aux) <- not.sigTerms
            algoRes <- c(algoRes, aux)

            .whichAlgorithm <- "parentchild"
            attr(.whichAlgorithm, "testClass") <- as.character(class(test.stat))
            return(new("topGOresult",
                       description = paste(description(object), "\nOntology:", ontology(object), sep = " "),
                       score = algoRes, testName = Name(test.stat),
                       algorithm = .whichAlgorithm,
                       geneData = c(.getGeneData(object), SigTerms = length(.sigTerms))))
          })



########################## CLASSIC algorithm ##########################
.sigGroups.classic <- function(GOlist, test.stat) {
  
  sigList <- sapply(GOlist,
                    function(term) {
                      ## we can't geve names .......
                      ## Name(gi)
                      members(test.stat) <- term
                      return(runTest(test.stat))
                    })
  
  return(sigList)
}


########################## ELIM algorithm ##########################
.sigGroups.elim <- function(goDAG, test.stat) {

  goDAG.r2l <- reverseArch(goDAG)
  nodeLevel <- buildLevels(goDAG.r2l, leafs2root = FALSE)
  levelsLookUp <- nodeLevel$level2nodes
  
  ##adjs.cutOff <- cutOff(test.stat) / numNodes(goDAG)
  adjs.cutOff <- cutOff(test.stat)
  #cat(paste("\n\t\t Parameters:\t cutOff = ", adjs.cutOff, "\n", sep =""))
  
  ## we use a lookup table to search for nodes that have were significant
  sigNodes.LookUP <- new.env(hash = T, parent = emptyenv())

  ## hash table for the genes that we eliminate for each node
  ## we store the genes that we want to eliminate
  elimGenes.LookUP <- new.env(hash = T, parent = emptyenv())
  
  ## hash table to store the result
  sigList <- new.env(hash = T, parent = emptyenv())


  for(i in nodeLevel$noOfLevels:1) {
    currNodes.names <- get(as.character(i), envir = levelsLookUp, mode = 'character')
    ##children.currNodes <- adj(goDAG.r2l, currNodes.names)
    currAnno <- .genesInNode(goDAG, currNodes.names)

    .num.elimGenes <- length(unique(unlist(as.list(elimGenes.LookUP))))
    cat(paste("\n\t Level ", i, ":\t", length(currNodes.names),
              " nodes to be scored\t(", .num.elimGenes, " eliminated genes)\n", sep =""))
    
    for(termID in currNodes.names) {
      ## just in case the test.stat is modified by some methods
      group.test <- updateGroup(test.stat, name = termID, members = currAnno[[termID]])
      
      ## remove the genes from the term, if there are
      if(!exists(termID, envir = elimGenes.LookUP, mode = "character"))
        elim(group.test) <- character(0)
      else
        elim(group.test) <- get(termID, envir = elimGenes.LookUP, mode = 'character')

      ## run the test and store the result (the p-value)
      termSig <- runTest(group.test)
      assign(termID, termSig, envir = sigList)
      
      ## if we have a significant GO term 
      if(termSig <= adjs.cutOff) {
        ## we mark it
        assign(termID, termSig, envir = sigNodes.LookUP)

        ## we set the genes that we want to eliminate from the all ancestors
        elimGenesID <- currAnno[[termID]]

        ## we look for the ancestors
        ancestors <- setdiff(nodesInInducedGraph(goDAG, termID), termID)
        
        oldElimGenesID <- mget(ancestors, envir = elimGenes.LookUP,
                               mode = 'character', ifnotfound = list(NA))

        ## add the new genesID to the ancestors
        newElimGenesID <- lapply(oldElimGenesID,
                                 function(termGenes) {
                                   aux <- union(termGenes, elimGenesID)
                                   return(aux[!is.na(aux)])
                                 })

        
        ## update the lookUp table
        if(length(newElimGenesID) > 0)
          multiassign(names(newElimGenesID), newElimGenesID, envir = elimGenes.LookUP)
      }
    }
  }
  
  return(unlist(as.list(sigList)))
}


########################## WEIGHT algorithm ##########################
.sigGroups.weight <- function(goDAG, test.stat) {

  ## SOME FUNCTIONS THAT WE NEED
  
  ## make a join between x and y unsing the names as key
  ## combFun should be a vectorial function
  ## we can use functions like: pmin, *, pmax, ...
  x.comb.y <- function(x, y, combFun = get('*')) {
    
    ## check if we weight the same set: always the case when we
    ## set the genes in downNodes.LookUP
    nx <- names(x)
    ny <- names(y)
    if((length(nx) == length(ny)) && all(nx == ny))
      return(combFun(x, y))
    
    allNames <- union(nx, ny)
    commonNames <- intersect(nx, ny)
    
    z <- numeric(length(allNames))
    names(z) <- allNames
    
    z[nx] <- x
    z[ny] <- y
    z[commonNames] <- combFun(x[commonNames], y[commonNames])
    
    return(z)
  }


  setGeneWeights <- function(termID, geneWeights, LookUP.table) {
    oldWeigths <- get(termID, envir = LookUP.table)
    assign(termID, x.comb.y(geneWeights, oldWeigths), envir = LookUP.table)
  }

  
  ## this is a recursive function in termChildren
  computeTermSig <- function(group.test, termChildren) {

    termID <- Name(group.test)
    ## look to see if there are some genes that must be weighted
    Weights(group.test) <- get(termID, envir = upNodes.LookUP)    

    ## run the test (the p-value)
    termSig <- runTest(group.test)

    ## if we came from some recursive call
    if(length(termChildren) == 0) {
      assign(termID, termSig, envir = sigList)
      return()
    }
    
    ## since we konw that length(termChildren) > 0
    w <- pW <- numeric(length(termChildren))
    names(w) <- names(pW) <- termChildren
    
    for(child in termChildren) {
      ## look for the child significance
      ## we should have the child significance in sigList
      if(!exists(child, envir = sigList, inherits = FALSE))
        stop('the child of node u, should had been processed')
      childSig <- get(child, envir = sigList)
      
      w[child] <- getSigRatio(group.test, a = childSig, b = termSig)
      pW[child] <- penalise(group.test, a = childSig, b = termSig)
    }
    
    ## if w[child] > 1 than that child is more significant
    sig.termChildren <- names(w[w > 1])
    
    ## CASE 1:  if we don't have significant children
    ##          in this case we down-weight the genes in each child
    if(length(sig.termChildren) == 0) {
      for(child in termChildren) {
        ## the genes in this child
        gene.child <- currAnno[[child]]
        
        ## weights for the genes in this child
        gene.childWeights <- rep(pW[child] * w[child], length(gene.child))
        names(gene.childWeights) <- gene.child

        ## put the gene weights in the downNodes.LookUP hash-table
        if(!exists(child, envir = downNodes.LookUP, inherits = FALSE))
          gCW <- gene.childWeights
        else {
          ## if for child exists some genes that are weighted
          ## they are 'combine' using "pmin" function
          oldWeigths <- get(child, envir = downNodes.LookUP)
          gCW <- x.comb.y(gene.childWeights, oldWeigths, combFun = pmin)
        }

        assign(child, gCW, envir = downNodes.LookUP)
      }
      ## since we have only less significant genes that means
      ## that we don't recompute the significance for node 'u'
      assign(termID, termSig, envir = sigList)
      return()
    }
    
    ## CASE 2:   'child' is more significant that 'u'
    .nn <- members(group.test)
    geneWeights <- numeric(length(.nn)) + 1
    names(geneWeights) <- .nn

    for(child in sig.termChildren) {
      ## the genes in this child
      gene.child <- currAnno[[child]]
      
      gene.childWeights <- rep(pW[child] / w[child], length(gene.child))
      names(gene.childWeights) <- gene.child

      geneWeights <- x.comb.y(gene.childWeights, geneWeights)
    }

    ## we set the new weights only for node 'termID'
    setGeneWeights(termID, geneWeights, upNodes.LookUP)

    if(is.null(.gene.weights))
      assign('.gene.weights', geneWeights, envir = .main.envir)
    else
      assign('.gene.weights', x.comb.y(.gene.weights, geneWeights), envir = .main.envir)
    
    computeTermSig(group.test, setdiff(termChildren, sig.termChildren))
  }
  
  ##----------------------------------------------------------------------

  goDAG.r2l <- reverseArch(goDAG)
  ## get the levels list 
  nodeLevel <- buildLevels(goDAG.r2l, leafs2root = FALSE)
  levelsLookUp <- nodeLevel$level2nodes

  ## the annotations of the genes to GO
  currAnno <- .genesInNode(goDAG, nodes(goDAG))

  ## we use a lookup table to search for nodes that we have weighted the genes 
  aux <- lapply(currAnno,
                function(x) {
                  gw <- numeric(length(x)) + 1
                  names(gw) <- x
                  return(gw)
                })
  upNodes.LookUP <- l2e(aux, new.env(hash = T, parent = emptyenv()))
  downNodes.LookUP <-  new.env(hash = T, parent = emptyenv())

  ## we also use a hash table to store the result
  sigList <- new.env(hash = T, parent = emptyenv())
  
  ## we will use frequently the children
  allChildren <- adj(goDAG.r2l, nodes(goDAG))
  

  ## START
  .main.envir <- environment()
  
  for(i in nodeLevel$noOfLevels:1) {
    currNodes.names <- get(as.character(i), envir = levelsLookUp, mode = 'character')

    ## some messages
    cat(paste("\n\t Level ", i, ":\t", length(currNodes.names), " nodes to be scored.\n", sep =""))

    for(termID in currNodes.names) {
      ## take the IDs of the gene in the current node
      ## we store the termID in the "name" slot of the test.stat
      group.test <- updateGroup(test.stat, name = termID, members = currAnno[[termID]])

      .gene.weights <- NULL
      ## look at all children of node u and compute the significance
      computeTermSig(group.test, allChildren[[termID]])

      if(!is.null(.gene.weights)) {
        ## get the upper induced subgraph from 'u' (including 'u')
        nodesUpSubgraph <- nodesInInducedGraph(goDAG, termID)
        ## set the weights for all the nodes
        
        lapply(setdiff(nodesUpSubgraph, termID), setGeneWeights, .gene.weights, upNodes.LookUP)
      }
    }
  }

  
  ## recompute the significance of the nodes in the downNodes.LookUP
  for(termID in ls(downNodes.LookUP)) {
    newWeights <- x.comb.y(get(termID, envir = upNodes.LookUP),
                           get(termID, envir = downNodes.LookUP))
    
    termSig <- runTest(updateGroup(test.stat, name = termID,
                                   members = currAnno[[termID]], weights = newWeights))
    
    assign(termID, termSig, envir = sigList)
  }
  
  return(unlist(as.list(sigList)))
}






########################## WEIGHT algorithm ##########################
.sigGroups.weight01 <- function(goDAG, test.stat) {

  ## SOME FUNCTIONS THAT WE NEED
  
  removeGenes <- function(termID, genesID, LookUP.table) {
    if(!exists(termID, envir = LookUP.table, mode = "character"))
      assign(termID, genesID, envir = LookUP.table)
    else
      assign(termID, union(genesID, get(termID, envir = LookUP.table)), envir = LookUP.table)
  }
  
  
  ## this is a recursive function in termChildren
  computeTermSig <- function(group.test, termChildren) {

    termID <- Name(group.test)

    ## remove the genes from the term, if there are
    if(!exists(termID, envir = upNodes.LookUP, mode = "character"))
      elim(group.test) <- character(0)
    else
      elim(group.test) <- get(termID, envir = upNodes.LookUP, mode = 'character')

    ## run the test (the p-value)
    termSig <- runTest(group.test)

    ## if we came from some recursive call
    if(length(termChildren) == 0) {
      assign(termID, termSig, envir = sigList)
      return()
    }
    
    ## since we konw that length(termChildren) > 0
    w <- numeric(length(termChildren))
    names(w) <- termChildren
    
    for(child in termChildren) {
      ## look for the child significance
      if(!exists(child, envir = sigList, inherits = FALSE))
        stop('the child of node u, should had been processed')
      ## we should have the child significance in sigList
      childSig <- get(child, envir = sigList)
      
      w[child] <- .sigRatio.01(a = childSig, b = termSig)
    }
    
    ## if w[child] > 1 than that child is more significant
    sig.termChildren <- names(w[w > 1])
    
    ## CASE 1:  if we don't have significant children 
    if(length(sig.termChildren) == 0) {
      ## since we have only less significant nodes that means
      ## that we don't recompute the significance for node 'u'
      assign(termID, termSig, envir = sigList)
      return()
    }
    
    ## CASE 2:   'child' is more significant that 'u'
    geneToRemove <- unique(unlist(currAnno[sig.termChildren]))

    ## we remove the genes only for node 'termID'
    removeGenes(termID, geneToRemove, upNodes.LookUP)

    if(is.null(.genesToRemove))
      assign('.genesToRemove', geneToRemove, envir = .main.envir)
    else
      assign('.genesToRemove', union(.genesToRemove, geneToRemove), envir = .main.envir)
    
    computeTermSig(group.test, setdiff(termChildren, sig.termChildren))
  }
  
  ##----------------------------------------------------------------------

  goDAG.r2l <- reverseArch(goDAG)
  ## get the levels list 
  nodeLevel <- buildLevels(goDAG.r2l, leafs2root = FALSE)
  levelsLookUp <- nodeLevel$level2nodes

  ## the annotations of the genes to GO
  currAnno <- .genesInNode(goDAG, nodes(goDAG))

  ## hash table for the genes that we eliminate for each node
  ## we store the genes that we want to eliminate
  upNodes.LookUP <- new.env(hash = T, parent = emptyenv())

  ## we also use a hash table to store the result
  sigList <- new.env(hash = T, parent = emptyenv())
  
  ## we will use frequently the children
  allChildren <- adj(goDAG.r2l, nodes(goDAG))
  

  ## START
  .main.envir <- environment()
  
  for(i in nodeLevel$noOfLevels:1) {
    currNodes.names <- get(as.character(i), envir = levelsLookUp, mode = 'character')

    ## some messages
    .num.elimGenes <- length(unique(unlist(as.list(upNodes.LookUP))))
    cat(paste("\n\t Level ", i, ":\t", length(currNodes.names),
              " nodes to be scored\t(", .num.elimGenes, " eliminated genes)\n", sep =""))
    
    for(termID in currNodes.names) {
      ## take the IDs of the gene in the current node
      ## we store the termID in the "name" slot of the test.stat
      group.test <- updateGroup(test.stat, name = termID, members = currAnno[[termID]])

      .genesToRemove <- NULL
      ## compute the significance of node u by looking at all children 
      computeTermSig(group.test, allChildren[[termID]])

      if(!is.null(.genesToRemove)) {
        ## get the nodes in the upper induced subgraph from 'u' (including 'u')
        ## and eliminate the genes for all these nodes
        lapply(setdiff(nodesInInducedGraph(goDAG, termID), termID),
               removeGenes, .genesToRemove, upNodes.LookUP)
      }
    }
  }
  
  return(unlist(as.list(sigList)))
}





########################## Parent-Child algorithm ##########################
.sigGroups.parentChild <- function(goDAG, test.stat, sigGenes) {

  nodeLevel <- buildLevels(goDAG, leafs2root = TRUE)
  levelsLookUp <- nodeLevel$level2nodes
  
  ## hash table to store the result
  sigList <- new.env(hash = T, parent = emptyenv())

  ## we don't score the root node 
  for(i in nodeLevel$noOfLevels:2) {
    currNodes.names <- get(as.character(i), envir = levelsLookUp, mode = 'character')
    parents.currNodes <- adj(goDAG, currNodes.names)
    currAnno <- .genesInNode(goDAG, currNodes.names)

    cat(paste("\n\t Level ", i, ":\t", length(currNodes.names),
              " nodes to be scored.\n", sep =""))
    
    for(termID in currNodes.names) {
      ## just in case the test.stat is modified by some methods
      parents <- .genesInNode(goDAG, parents.currNodes[[termID]])
      group.test <- updateGroup(test.stat, members = currAnno[[termID]],
                                parents = parents, sigMembers = sigGenes)
      
      ## run the test and store the result (the p-value)
      termSig <- runTest(group.test)
      assign(termID, termSig, envir = sigList)
    }
  }

  ## the root will have p-value 1
  assign(get("1", envir = levelsLookUp, mode = 'character'), 1, envir = sigList)
  
  return(unlist(as.list(sigList)))
}

