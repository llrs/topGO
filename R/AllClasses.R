#################### topGOdata class ####################
## 
## The node atributes are environments containing the genes/probes
## annotated to the respective node
##
## If genes is a numeric vector than this should represent the gene's
## score. If it is factor it should discriminate the genes in
## interesting genes and the rest

## TODO: it will be a good idea to replace the allGenes and allScore with an
## ExpressionSet class. In this way we can use tests like global test, globalAncova....
##  -- ALL variables sarting with . are just for internal class usage (private)


setClass("topGOdata",
         representation = representation(
           ## some description of the data/experiment
           description = "character",
           ## which of the ontology to be used: BP, CC, MF, (BC, BM, CM, ...) 
           ontology = "character",
           ## the vector containing the genes/probes used in the experiment
           allGenes = "character",
           ## the vector containing the score or the  of the  genes/probes 
           allScores = "ANY",
           ## if each has a score, geneSelectionFun: function(x, ...)
           geneSelectionFun = "function",
           ## which genes are annotated to the specific ontology
           feasible = "logical",
           ## the GO ontology: graph structure and annotations: leaves to root
           graph = "graphNEL",
           ## expression matrix                                     #!
           expressionMatrix = "ANY",                                #!
           ## phenotype information (groups that shall be compared) #!
           phenotype = "ANY"))                                      #!



######################################################################

## probably will add more infos with time here
setClass("topGOresult",
         representation = representation(
           ## some description of the data/experiment
           description = "character",
           ## the p-values of the GO terms (named vector)
           score = "numeric",
           ## which test statistic was used 
           testName = "character",
           ## which class was used 
           testClass = "character"))



######################## groupInfo class ########################
## A virtual class containing basic group (GO term) data:
##     gene names, genes scores, etc...

setClass("groupStats",
         representation = representation(
           ## the name of the group (GO ID for example)
           name = "character",
           ## the names of all the mebers (gene ID)
           allMembers = "character",
           ## the names of the mebers in the group (gene ID)
           members = "character",
           ## function containg a statistical test takeing the
           ## parameter the object itself 
           testStatistic = "function",
           testStatPar = "list",
           "VIRTUAL"))

#################### classicCount class ####################
## A class that extends the virtual class groupStats by adding 
## a slot representing the significant members.

setClass("classicCount", contains = "groupStats", 
         representation = representation(
           ## the index of which members are significant  
           significant = "integer"))


#################### classicScore class ####################
## A class that extends the virtual class groupStats by adding 
## a slot representing the score of each gene. (used for KS test)

setClass("classicScore", contains = "groupStats",
         representation = representation(
           ## the score for each member, the most important
           ## member has the highest score
           score = "numeric",
           ## scoreOrder = TRUE (decreasing) the max(score) is considering the best score
           ## scoreOrder = FALSE (increasing) the min(score) is considre the best score
           scoreOrder = "logical"))


######################################################################


#################### elimCount class ####################

setClass("elimCount", contains = "classicCount", 
         representation = representation(
           ## the index of which of the group members should be removed
           elim = "integer",
           cutOff = "numeric"))


#################### elimScore class ####################

setClass("elimScore", contains = "classicScore",
         representation = representation(
           ## the score for each member, the most important
           ## member has the highest score
           elim = "integer",
           cutOff = "numeric")) 



######################################################################


#################### weightCount class ####################

setClass("weightCount", contains = "classicCount", 
         representation = representation(
           weights = "numeric",
           sigRatio = "function",
           ## if more sig nodes ar penalized more (should be part of the sigRatio) 
           penalise = "function",  ## penalise = c("more", "less")
           roundFun = "function"))


#################### weightScore class ####################

setClass("weightScore", contains = "classicScore",
         representation = representation(
           ## the score for each member, the most important
           ## member has the highest score
           weights = "numeric"))

