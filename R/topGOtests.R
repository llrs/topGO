## put here all tests statistic

if(!isGeneric("GOFisherTest"))
  setGeneric("GOFisherTest", function(object) standardGeneric("GOFisherTest"))
                               
setMethod("GOFisherTest", "classicCount",
          function(object) {

            contMat <- contTable(object)
            
            if(all(contMat == 0))
              p.value <- 1
            else
              p.value <- fisher.test(contMat, alternative = 'greater')$p.value

            return(p.value)
          })



if(!isGeneric("GOKSTest"))
  setGeneric("GOKSTest", function(object) standardGeneric("GOKSTest"))
                               
setMethod("GOKSTest", "classicScore",
          function(object) {

            N <- numAllMembers(object)
            na <- numMembers(object)

            ## if the group is empty ... (should not happen, but you never know!)
            if(na == 0 || na == N)
              return(1)

            x.a <- rankMembers(object)
            x.b <- setdiff(1:N, x.a)

            return(ks.test(x.a, x.b, alternative = "greater")$p.value)
          })



if(!isGeneric("GOtTest"))
  setGeneric("GOtTest", function(object) standardGeneric("GOtTest"))
                               
setMethod("GOtTest", "classicScore",
          function(object) {

            ## we should try to do this smarter!!! 
            var.equal <- FALSE
            ## check if there are any parameters set
            if(length(testStatPar(object)) != 0)
              if("var.equal" %in% names(testStatPar(object)))
                var.equal <- testStatPar(object)[["var.equal"]]
                                        
            nG <- numMembers(object)
            nNotG <- numAllMembers(object) - nG
            
            ## if the group complementary is empty 
            if(nNotG == 0 || nG == 0)
              return(1)

            x.G <- membersScore(object)
            x.NotG <- allScore(object, TRUE)[setdiff(allMembers(object), members(object))]

            aa <- alternative(object)
            ## alternative - TRUE if the max(score) is the best score), FALSE if the min(score)...
            if(nNotG == 1) {
              .stat <- ifelse(aa, sum(x.G >= x.NotG), sum(x.G <= x.NotG))
              p.value <- (.stat + 1) / (nG + nNotG + 1)
            }
            else if(nG == 1) {
              .stat <- ifelse(aa, sum(x.NotG >= x.G), sum(x.NotG <= x.G))
              p.value <- (.stat + 1) / (nG + nNotG + 1)
            }
            else 
              p.value <- t.test(x = x.G, y = x.NotG, var.equal = var.equal,
                                alternative = ifelse(aa, "greater", "less"))$p.value
            
            return(p.value)
          })



if(!isGeneric("GOglobalTest"))
  setGeneric("GOglobalTest", function(object) standardGeneric("GOglobalTest"))
                               
setMethod("GOglobalTest", "classicExpr",
          function(object) {
            
            if(numMembers(object) == 0)
              return(1)

            return(p.value(globaltest(X = membersExpr(object), Y = pType(object))))
          })


