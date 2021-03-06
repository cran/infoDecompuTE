##' Construct the Matrix from Information Decomposition and Compute the
##' Efficiency Factors of Treatment effects
##' 
##' Perform the information decomposition for either the block or treatment
##' effects within a single stratum and Compute the Efficiency Factors for
##' every treatment effect within a single stratum.
##' 
##' The main purpose of this function is to construct a list of resultant
##' matrices associated with each source of variation after the information
##' decomposition and to compute the canonical or average efficiency factors
##' for each treatment effects in each stratum of ANOVA table.
##' 
##' The canonical efficiency factors are generated when the user input the
##' treatment contrasts, otherwise the average efficiency factors, which is the
##' harmonic mean of the canonical efficiency factors, are generated.
##' 
##' @param z a matrix containing the orthogonal projector of a stratum
##' generated by \code{\link{makeOrthProjectors}}.
##' @param T a list of contrast matrices generated by
##' \code{\link{makeContrMat}}.
##' @param N a matrix containing the design matrix generated by
##' \code{\link{makeOverDesMat}}.
##' @param Rep a matrix containing the treatment replication number and is
##' generated by \code{\link{getTrtRep}}.
##' @param trt.Sca a numeric vector for computing a coefficients of the fixed
##' effect parameter in EMS and is generated by \code{\link{getTrtRep}}.
##' @return A list of matrices and numeric vectors containing the efficiency
##' factors of every treatment effect.
##' @author Kevin Chang
##' @examples
##' 
##' design1 <- local({ 
##'   Ani = as.factor(LETTERS[c(1,2,3,4,
##'                             5,6,7,8)])
##'   Trt = as.factor(letters[c(1,1,1,1,
##'                             2,2,2,2)])
##'   data.frame(Ani, Trt, stringsAsFactors = TRUE )
##' })
##' 
##' blk.str = "Ani"
##'     
##' rT = terms(as.formula(paste("~", blk.str, sep = "")), keep.order = TRUE) 
##' blkTerm = attr(rT,"term.labels")
##'      
##' Z = makeBlkDesMat(design1, blkTerm)
##' 
##' trt.str = "Trt"              
##' fT <- terms(as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  #fixed terms
##' 
##' trtTerm <- attr(fT, "term.labels")
##' effectsMatrix <- attr(fT, "factor")        
##' 
##' T <- makeContrMat(design1, trtTerm, effectsMatrix, contr.vec = NA)
##' 
##' N =  makeOverDesMat(design1, trtTerm)
##' 
##' Replist = getTrtRep(design1, trtTerm)   
##'  
##' Rep <- Replist$Rep
##' trt.Sca <- Replist$Sca
##'     
##' effFactors = lapply(makeOrthProjectors(Z), function(z) 
##'       getEffFactor(z, T, N, Rep, trt.Sca))
##' 
##' 
##' @export getEffFactor
getEffFactor <- 
  function(z, T, N, Rep, trt.Sca) {
  #z is Q, T is C and N is X in the paper
  if (!is.matrix(z)) 
    return(z)
  
  nEffect <- length(T)
  PNTginvATNP <- effFactors <- vector("list", nEffect)
  
  names(PNTginvATNP) <- names(effFactors) <- names(T)
  
  #First project into the first treatment vector subspace of stratum "z"
  PNTginvATNP[[1]] <- z %*% N %*% T[[1]] %*% 
    invInfMat(C = z, N = N, T = T[[1]]) %*% 
    T[[1]] %*% t(N) %*% t(z)

  
  if (!all(PNTginvATNP[[1]] < 1e-06)) {
    effFactors[[1]] <- vector("list", nEffect)
    names(effFactors[[1]]) <- names(T)
    for (i in 1:nEffect) {
      
      r.adjust <- MASS::ginv(sqrt(diag(Rep[, i])))
      #browser()
      
      # eigenvalues of the information matrix
      va <- Re(eigen(r.adjust %*% T[[i]] %*% t(N) %*% 
                       PNTginvATNP[[1]] %*% 
                       N %*% T[[i]] %*% r.adjust)$va)
      
      va <- va[which(va > 1e-07)]
      
      trt.coef <- Re(eigen(T[[i]] %*% t(N) %*% 
                            PNTginvATNP[[1]] %*% N 
                           )$va)[1:length(va)]/trt.Sca[i]

      if (isTRUE(all.equal(as.numeric(outer(va, va, "/")), 
                           rep(1, length(va) * length(va))))) {
        
        # harmonic means of the canonical efficiency factors 
        # to give the average efficiency factor
        effFactors[[1]][[i]] <- c(1/mean(1/va), trt.coef[1])
        
      } else {
        effFactors[[1]][[i]] <- c(1/mean(1/va), trt.coef)
      }
    }
  }
  
  newZ <-  (z %*% t(z)) - PNTginvATNP[[1]]
  
  if (nEffect != 1) {
    for (i in 2:nEffect) {
      
      PNTginvATNP[[i]] <- newZ %*% N %*% T[[i]] %*% 
        invInfMat(C = newZ, N = N, T = T[[i]]) %*% T[[i]] %*% 
        t(N) %*% t(newZ)
      
      if (all(PNTginvATNP[[i]] < 1e-06)) 
        next
      
      effFactors[[i]] <- vector("list", nEffect)
      names(effFactors[[i]]) <- names(T)
      
      for (j in 1:nEffect) {
        
        r.adjust <- MASS::ginv(sqrt(diag(Rep[, j])))
        
        va <- Re(eigen(r.adjust %*% T[[j]] %*% t(N) %*% 
                         PNTginvATNP[[i]] %*% N %*% r.adjust)$va)
        
        va <- va[which(va > 1e-07)]
        
        trt.coef <- Re(eigen(T[[j]] %*% t(N) %*% PNTginvATNP[[i]] %*% 
                              N %*% T[[j]])$va)[1:length(va)]/trt.Sca[j]
        
        if(isTRUE(all.equal(as.numeric(outer(va, va, "/")), 
                             rep(1, length(va) * length(va))))) {
          
          # harmonic means of the canonical efficiency factors 
          # to give the average efficiency factor
          
          effFactors[[i]][[j]] <- c(1/mean(1/va), trt.coef[1])
          
        } else {
          effFactors[[i]][[j]] <- c(1/mean(1/va), trt.coef)
        }
        
      }
      newZ <- (newZ %*% t(newZ)) - PNTginvATNP[[i]]
      
    }
  }
      
  PNTginvATNP$Residual <- newZ
  elementToRemove <- numeric(0)
  
  for (i in 1:length(PNTginvATNP)) {
    if (all(abs(PNTginvATNP[[i]]) < 1e-06)) #Removing the matrix with elements all equal to zero
      elementToRemove <- c(elementToRemove, i)
  }
  
  if (length(elementToRemove) > 0) 
    PNTginvATNP <- PNTginvATNP[-elementToRemove]
  
  
  return(list(PNTginvATNP, effFactors))
} 
