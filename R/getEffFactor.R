getEffFactor <- 
  function(z, T, N, Rep, trt.Sca) {
  
  if (!is.matrix(z)) 
    return(z)
  
  nEffect <- length(T)
  PNTginvATNP <- effFactors <- vector("list", nEffect)
  
  names(PNTginvATNP) <- names(effFactors) <- names(T)
  
  PNTginvATNP[[1]] <- z %*% N %*% T[[1]] %*% invInfMat(C = z, N = N, T = T[[1]]) %*% T[[1]] %*% t(N) %*% t(z)
  
  # browser()
  
  if (!all(PNTginvATNP[[1]] < 1e-06)) {
    effFactors[[1]] <- vector("list", nEffect)
    names(effFactors[[1]]) <- names(T)
    for (i in 1:nEffect) {
      r.adjust <- ginv(sqrt(diag(Rep[, i])))
      # eigenvalues of the information matrix
      va <- Re(eigen(r.adjust %*% T[[i]] %*% t(N) %*% PNTginvATNP[[1]] %*% N %*% T[[i]] %*% r.adjust)$va)
      
      va = va[which(va > 1e-07)]
      
      trt.coef = Re(eigen(T[[i]] %*% t(N) %*% PNTginvATNP[[1]] %*% N %*% T[[i]])$va)[1:length(va)]/trt.Sca[i]
      
      if (isTRUE(all.equal(as.numeric(outer(va, va, "/")), rep(1, length(va) * length(va))))) {
        
        # harmonic means of the canonical efficiency factors to give the average efficiency factor
        effFactors[[1]][[i]] <- c(1/mean(1/va), trt.coef[1])
        
      } else {
        effFactors[[1]][[i]] <- c(1/mean(1/va), trt.coef)
      }
    }
  }
  
  newZ <- (z %*% t(z)) - PNTginvATNP[[1]]
  
  if (nEffect != 1) {
    for (i in 2:nEffect) {
      
      PNTginvATNP[[i]] <- newZ %*% N %*% T[[i]] %*% invInfMat(C = newZ, N = N, T = T[[i]]) %*% T[[i]] %*% 
        t(N) %*% t(newZ)
      
      # PNTginvATNP[[i]] = newZ %*% N %*% t(T[[i]]) %*% ginv(t(T[[i]]) %*% t(N) %*% z %*% N %*% T[[i]]) %*% T[[i]]
      # %*% t(N) %*% t(newZ)
      
      if (all(PNTginvATNP[[i]] < 1e-06)) 
        next
      
      effFactors[[i]] <- vector("list", nEffect)
      names(effFactors[[i]]) <- names(T)
      for (j in 1:nEffect) {
        r.adjust <- ginv(sqrt(diag(Rep[, j])))
        va <- Re(eigen(r.adjust %*% T[[j]] %*% t(N) %*% PNTginvATNP[[i]] %*% N %*% T[[j]] %*% r.adjust)$va)
        
        va = va[which(va > 1e-07)]
        
        trt.coef = Re(eigen(T[[j]] %*% t(N) %*% PNTginvATNP[[i]] %*% N %*% T[[j]])$va)[1:length(va)]/trt.Sca[j]
        
        if (isTRUE(all.equal(as.numeric(outer(va, va, "/")), rep(1, length(va) * length(va))))) {
          
          # harmonic means of the canonical efficiency factors to give the average efficiency factor
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
    if (all(PNTginvATNP[[i]] < 1e-06)) 
      elementToRemove <- c(elementToRemove, i)
  }
  
  if (length(elementToRemove) > 0) 
    PNTginvATNP <- PNTginvATNP[-elementToRemove]
  
  
  
  return(list(PNTginvATNP, effFactors))
} 
