getVCs.onePhase <-
function(design.df, random.terms, fixed.terms, var.comp = NA, trt.contr = NA, table.legend = FALSE){
 
  library(MASS)   #need this package to perform ginv and factions functions
 
#######################################################################################  
  #All the pre-initiated functions   
  projMat <- function(X) X %*% ginv(t(X) %*% X) %*% t(X)  #projection matrix
  mK <- function(n) matrix(1/n, nrow=n, ncol=n)           #Averaging matrix
  mJ <- function(n) diag(n) - mK(n)                       #J matrix 
  tr <- function(X) sum(diag(X))                          #Trace operation
  m1 <- function(n) matrix(rep(1,n), nrow = n, ncol = 1)  #M1 matrix
  mI <- function(n) diag(n)                               #identity matrix
  
#######################################################################################  
  #check for factor names 
  isFactorNameNumeric <- function(levels) !as.logical( length(grep("[A-Z]|[a-z]", levels)) )

#######################################################################################  
  #make design matrix
  makeDesignMatrix <- function(nRows, design.df, col){    
    if(grepl(":", col)){ 
      factor <-as.factor( apply(design.df[,unlist(strsplit(col, ":"))], 1, function(x) paste(x, collapse =".")))
    }else{
      factor <- as.factor(design.df[,col])
    }  
        
    facName <- col
    nCols <- nlevels(factor)
   
    Z <- matrix(0, nrow=nRows, ncol=nCols)
    Z[cbind(1:nRows, match(c(factor), 1:nCols))] <- 1
    if(isFactorNameNumeric(levels(factor))){ 
      colNames <- paste(facName, 1:nCols,sep="")
    }else{                                    
      colNames <- levels(factor)
    }
    
    dimnames(Z) <- list(1:nRows, colNames)
    return(Z)
  } 
 
#######################################################################################  
 #make block design matrix 
 makeBlkDesMatrix <- function(design.df, blkTerm){

   # design.df = data.frame containing design
   # blkTerm = block terms

   n <- length(blkTerm)
   nRows <- nrow(design.df)
   Z <- list(NULL)
   Z[[1]] <- diag(nrow(design.df))
   
   for(i in 2:(n+1)){   
    Z[[i]] <- makeDesignMatrix(nRows=nRows, design.df=design.df, col=blkTerm[i-1]) 
   }   
        
   names(Z) <- c("e", blkTerm)  
   return(Z)
  }   
 
#######################################################################################  
 #make the block projection matrix
  makeBlockProjectors <- function(BlkDesignMatrixList, initial = diag(nrow(BlkDesignMatrixList[[1]])) - mK(nrow(BlkDesignMatrixList[[1]]))){
   
    n <- nrow(BlkDesignMatrixList$e)
    Q <- lapply(BlkDesignMatrixList, function(z) projMat(z))
    
    Q <- Q[sort(1:length(Q), decreasing=TRUE)]
    
    elementToRemove = numeric(0)

    cusumQ <- P <- NULL
    P[[1]] <- Q[[1]] %*% initial     
                                                                                          
    if(all(P[[1]] <1e-6)){
     elementToRemove = 1                       #revord the elements that are all zeros
     P[[1]] <- matrix(0, nrow = n, ncol = n)   #make the projectors which has less than 1e-6 to zero, to avoid rounding error
    } 
        
    cusumQ[[1]] <- initial - P[[1]]    
    
    for(i in 2:(length(Q))){
      P[[i]] <- Q[[i]] %*% cusumQ[[i-1]]

      if((all(P[[i]] <1e-6))|| tr(P[[i]]) <= 0){
        elementToRemove = c(elementToRemove,i )   #revord the elements that are all zeros
        P[[i]] <- matrix(0, nrow = n, ncol = n)   #make the projectors which has less than 1e-6 to zero, to avoid rounding error
      }                 
      cusumQ[[i]] <- cusumQ[[i-1]] - P[[i]]
    }     
       
    P <- P[sort(1:length(P), decreasing=TRUE)]
    names(P) <- names(BlkDesignMatrixList)   
   
    if(length(elementToRemove)>0)
      P= P[-(length(Q) - elementToRemove+1)]
                                    
    return(P)                   
  }


#######################################################################################  
  #transfrom each treatment contrast metrix to C matrix
transContrToT = function(fact, contr){
  T = matrix(0, nrow = nlevels(fact), ncol= nlevels(fact))
  if(is.matrix(contr)){
    rownames(contr) = fact
    for(i in 1:ncol(contr)){         
      x = contr[levels(fact),i]
      T = T + (x %*% t(x))/as.numeric(t(x) %*% x)
    }
  
  } else{
     names(contr) = fact
     x = contr[levels(fact)]
      T = T + (x %*% t(x))/as.numeric(t(x) %*% x)      
  }                                                     
  return(T)
}                                                                        

  
#######################################################################################  
   #make the treatment C matrix, this function allow the input of the treatment contrasts.
makeTreatProjectors <- function(design.df, trtCols, effectsMatrix, trt.contr){       
   
    if(length(trtCols) == 1){
      nLevels = nlevels(design.df[,trtCols])
      names(nLevels) =  trtCols   
    }else if(any(grepl(":", trtCols))){
      uniqueTrtCols = unique(unlist(strsplit(trtCols, "\\:")))  
      nLevels <- sapply(design.df[,uniqueTrtCols], function(x) nlevels(as.factor(x)))
    }else{
      nLevels <- sapply(design.df[,trtCols], function(x) nlevels(as.factor(x)))
    }  
  
   nLevels = nLevels[rownames(effectsMatrix)]
   nEffects = ncol(effectsMatrix)  
   
   
   if(all(is.na(trt.contr))){
     indMatrix <- function(x,n){
        if(x == 1) X <- mJ(n)
        else if(x == 2)  X <- mI(n)     
        else       X <- mK(n)
        return(X)
     } 
     
     X <- as.list(rep(1,nEffects))
     names(X)  <- colnames(effectsMatrix)
  
     for(i in 1:nrow(effectsMatrix)){      
        matList <- lapply(effectsMatrix[i,], function(y) indMatrix(y, nLevels[i]))     
        for(j in 1:nEffects) X[[j]] <- X[[j]] %x% matList[[j]]    
     }
   }else{
      indMatrix1 <- function(x,n, trtContr){
        if(x == 1) X <- trtContr
        else if(x == 2)  X <- mI(n)     
        else       X <- mK(n)
        return(X)
     } 
     names(trt.contr) = trtCols
     if(!any(sapply(trt.contr, class)=="list")){
        X <- as.list(rep(1,nEffects))
       names(X)  <- colnames(effectsMatrix)
       names(trt.contr)<- colnames(effectsMatrix)  
       
       for(i in 1:nrow(effectsMatrix)){
       
          trtContr = transContrToT(design.df[,rownames(effectsMatrix)[i]], 
          trt.contr[[rownames(effectsMatrix)[i]]])        
          matList <- lapply(effectsMatrix[i,], function(y) 
                                            indMatrix1(y, nLevels[i], trtContr))     
          for(j in 1:nEffects) X[[j]] <- X[[j]] %x% matList[[j]]    
       }
     } else {
      
      X = list()
      count = 1
      
      totalLength = length(trt.contr[-which(sapply(trt.contr, class)=="list")]) + 
       sum(sapply( trt.contr[which(sapply(trt.contr, class)=="list")], length))
       
      newEffectsMatrix = matrix(0, ncol =totalLength  , nrow = nrow(effectsMatrix))
      rownames(newEffectsMatrix)  =  rownames(effectsMatrix)
      
      for(i in 1:length(trt.contr)){
        if(is.list(trt.contr[[i]])){
          tmpX = rep(1,length(trt.contr[[i]]))
          names(tmpX) = paste(names(trt.contr[i]), names(trt.contr[[i]]),sep=".")
          
          tmpCounter =  count:(count+length(trt.contr[[i]])-1)
         
          X = c(X, tmpX)
         
          newEffectsMatrix[names(trt.contr[i]),tmpCounter] = 1          
          count = count+length(trt.contr[[i]])
        }else{
          tmpX = 1
          names(tmpX) = names(trt.contr[i])
          X = c(X, tmpX)
                    
          newEffectsMatrix[,count] = effectsMatrix[,i] 
          
          tmpCounter =  count = count + 1 
        }
      }
      
      colnames(newEffectsMatrix) = names(X)
      
      for(i in 1:nrow(effectsMatrix)){
       if(is.list(trt.contr[[rownames(effectsMatrix)[i]]])){
            trtContr = lapply(trt.contr[[rownames(effectsMatrix)[i]]], function(x)
                    transContrToT(design.df[,rownames(effectsMatrix)[i]], x))        
          
             matList <- vector(mode = "list", length = length(trt.contr))
            
            for( j in 1:ncol(newEffectsMatrix)){            
              if(newEffectsMatrix[i,j] == 1){  
                matList[[j]] <- trtContr[[match(unlist(strsplit(colnames(newEffectsMatrix)[j], "\\."))[2],
                                                          names(trtContr))]]                                   
              }else if(newEffectsMatrix[i,j] == 2){
                 matList[[j]] <-   mI(nLevels[i])                         
               }else{
                 matList[[j]] <- mK(nLevels[i])
              } 
            }
          }else{
            trtContr = transContrToT(design.df[,rownames(effectsMatrix)[i]], 
                        trt.contr[[rownames(effectsMatrix)[i]]])        
            matList <- lapply(newEffectsMatrix[i,], function(y) 
                    indMatrix1(y, nLevels[i], trtContr))     
          }
          for(j in 1:length(X)) X[[j]] <- X[[j]] %x% matList[[j]]       
        
    }
   }
   }
   return(X)  
  }
  
      
        
#######################################################################################  
  #get the incidence matrix N
  getIncidenceMatrix <- function(design.df, trtCols){
   
   if(length(trtCols) == 1){
      incident = design.df[,trtCols]
      nLevels = levels(design.df[,trtCols])    
    }else if(any(grepl(":", trtCols))){ 
      uniqueTrtCols = unique(unlist(strsplit(trtCols, "\\:")))  

      incident = as.factor(apply(design.df[,uniqueTrtCols], 1, function(x) paste(x, collapse =".")))
      nLevels = sort(levels(interaction(design.df[,uniqueTrtCols])))

    }else{
     incident = as.factor(apply(design.df[,trtCols], 1, function(x) paste(x, collapse =".")))
     nLevels = sort(levels(interaction(design.df[,trtCols])))
    }       
   
   N <- matrix(0, nrow=nrow(design.df), ncol=length(nLevels))
   N[cbind(1:nrow(design.df), match(incident, nLevels))] <- 1
  
   return(N)
  }

#######################################################################################  
  #get a list of replication for each treatment
  getReplicationList <- function(design.df, trtCols){
    
    if(length(trtCols) == 1){
      return(mean(table(design.df[,trtCols])))
      
    }else if(any(grepl(":", trtCols))){ 
     
      trtColsList = lapply(strsplit(trtCols, "\\:"), function(x) design.df[,x])
      repList = sapply(trtColsList, function(y) if(is.factor(y)) {mean(table(y))} else{ mean(table(apply(y, 1, function(x) paste(x, collapse ="."))))} )
        
      levelList =  sapply(trtColsList, function(y) if(is.factor(y)) {nlevels(y)} else{ nlevels(as.factor(apply(y, 1, function(x) paste(x, collapse ="."))))} )
        
      repList = repList * levelList/nlevels(interaction(design.df[,unique(unlist(strsplit(trtCols, "\\:")))  ]))
     
      return(repList)      
    
    }else{
       repList = apply(design.df[,trtCols], 2, function(x)mean(table(x)))
       levelList = apply(design.df[,trtCols], 2, function(x) nlevels(as.factor(x)))
       repList = repList* levelList/ nlevels(interaction(design.df[,trtCols]))
     return(repList)      
        
    }                                                                 
  }    

#######################################################################################  
  #get a list of coeffient for each treatment
  getTrtCoef <- function(design.df, trtCols){
    
    if(length(trtCols) == 1){
      return(length(design.df[,trtCols])/nlevels(design.df[,trtCols]))
    }else if(any(grepl(":", trtCols))){ 
     
      trtColsList = lapply(strsplit(trtCols, "\\:"), function(x) design.df[,x])
      repList = sapply(trtColsList, function(y) if(is.factor(y)) {mean(table(y))} else{ mean(table(apply(y, 1, function(x) paste(x, collapse ="."))))} )
       
      return(repList)      
    }else{
       repList = apply(design.df[,trtCols], 2, function(x)mean(table(x)))
      return(repList)      
    }                                                                 
  }    

#######################################################################################  
  # Compute g-inverse of information matrix
  invInfMat <- function(C,N,T){

   ei <- eigen(t(T) %*% t(N) %*% C %*% N %*% T)
   nn <- length(ei$values)
   L <- matrix(0, nrow=nn, ncol=nn)
   for(i in 1:(nn)) {
    if( Re(ei$values[i]) <1e-6)next
    L <- L + (1/Re(ei$values[i]))*Re(ei$vectors[,i])%*%t(Re(ei$vectors[,i]))
   }
   return(L)
  }
#######################################################################################  
  #get eigenvalue using single value decomposition 
  eigenValue <- function(C,N,T){
   fractions(eigen(t(T) %*% t(N) %*% C %*% N %*% T)$va[1])
  }
  

#######################################################################################  
  #pre- and post-multiply NTginvATN by block projection matrices for each stratum 
  blkProkMat = function(z, T, N, Rep){
      if(!is.matrix(z)) return(z)
      
      nEffect = length(T)
      PNTginvATNP = T
      
      PNTginvATNP[[1]] = z %*% N %*% T[[1]] %*% invInfMat(C=z, N=N, T=T[[1]]) %*% T[[1]] %*% t(N) %*% t(z) 
     
      newZ = (z %*% t(z)) - PNTginvATNP[[1]]
      
      if(nEffect !=1){            
        for(i in 2:nEffect){                      
          PNTginvATNP[[i]] = newZ %*% N %*% t(T[[i]]) %*% invInfMat(C=newZ, N=N, T=T[[i]]) %*% T[[i]] %*% t(N) %*% t(newZ)    
              newZ = (newZ %*% t(newZ)) - PNTginvATNP[[i]]
        }
      }

      PNTginvATNP$Residual = newZ
      elementToRemove = numeric(0)
      for(i in 1:length(PNTginvATNP)){        
        if(all(PNTginvATNP[[i]] <1e-6))
          elementToRemove = c(elementToRemove, i)
      }
    
      if(length(elementToRemove)>0)        
        PNTginvATNP= PNTginvATNP[-elementToRemove] 
                                     
      return(PNTginvATNP) 
    }
    
######################################################################################  
  #pre- and post-multiply NTginvATN by block projection matrices produces the effFactors for each stratum  
  trtProkMat = function(z, T, N, Rep){
      if(!is.matrix(z)) return(z)
      
      nEffect = length(T)
      PNTginvATNP = effFactors = vector("list", nEffect)
      
      names(PNTginvATNP) = names(effFactors) = names(T)
      
      PNTginvATNP[[1]] = z %*% N %*% T[[1]] %*% invInfMat(C=z, N=N, T=T[[1]]) %*% T[[1]] %*% t(N) %*% t(z) 
      if(!all(PNTginvATNP[[1]] <1e-6)){    
        effFactors[[1]] = vector("list", nEffect)
        names(effFactors[[1]]) = names(T)
        for(i in 1:nEffect){
          va = Re(eigen(T[[i]] %*% t(N) %*% PNTginvATNP[[1]] %*% N %*% T[[i]])$va)
          effFactors[[1]][[i]] = mean(va[which(va>1e-6)]/Rep[names(T[i])])
        } 
      }
      
      newZ = (z %*% t(z)) - PNTginvATNP[[1]]
      
      if(nEffect !=1){            
        for(i in 2:nEffect){                      
          PNTginvATNP[[i]] = newZ %*% N %*% t(T[[i]]) %*% invInfMat(C=newZ, N=N, T=T[[i]]) %*% T[[i]] %*% t(N) %*% t(newZ)    
          if(all(PNTginvATNP[[i]] <1e-6)) next
   
          effFactors[[i]] = vector("list", nEffect)
           names(effFactors[[i]]) = names(T)
          for(j in 1:nEffect){
             va = Re(eigen(T[[j]] %*% t(N) %*% PNTginvATNP[[i]] %*% N %*% T[[j]])$va)
            effFactors[[i]][[j]] =mean(va[which(va>1e-6)]/Rep[names(T[j])])
          }
          newZ = (newZ %*% t(newZ)) - PNTginvATNP[[i]]
        }
      } 
                                     
      return(effFactors) 
    }
   
#########################################################################################  
#Main methods starts here->
#Extract the fixed and random terms   
  
  rT = terms(as.formula(paste("~", random.terms, sep = "")), keep.order = TRUE) #random terms 
  fT = terms(as.formula(paste("~", fixed.terms, sep = "")), keep.order = TRUE)  #fixed terms 
  
#########################################################################################  
#Preparing the block structures    
  write("1. Preparing the block structure.", "")

  blkTerm = attr(rT,"term.labels")
  Z = makeBlkDesMatrix(design.df, rev(attr(rT,"term.labels")))
  
  Pb <- makeBlockProjectors(Z)
  if(names(Pb)[1] == "e")
     names(Pb)[1] = "Within"

#########################################################################################  
#Prepating the treatment structures
  write("2. Preparing the treatment structure.", "")


  trtTerm = attr(fT,"term.labels")
  effectsMatrix = attr(fT,"factor")
  
                                                 
  T =  makeTreatProjectors(design.df, trtTerm, effectsMatrix, trt.contr)
  N =  getIncidenceMatrix(design.df, trtTerm)
  Rep = getReplicationList(design.df, trtTerm)
  trt.Coef = getTrtCoef(design.df, trtTerm)  
  
  if(any(grepl("\\.", names(T)))){  
    names(Rep) = trtTerm
    names(trt.Coef) = trtTerm
    
    Rep = Rep[sapply(strsplit(names(T), "\\."), function(x) x[1])]
    trt.Coef = trt.Coef[sapply(strsplit(names(T), "\\."), function(x) x[1])]
  }
  
  names(Rep) = names(T)
  names(trt.Coef) = names(T)
#########################################################################################  
#Start calculating the VCs
#1-phase experiment
write("3. Start calculating the variance components.", "")
write("4. Pre- and post-multiply NTginvATN by block projection matrices.", "")

  #pre- and post-multiply NTginvATN by block projection matrices  
   PNTginvATNP <- lapply(Pb,function(z) blkProkMat(z, T, N))
  
   #Now construct variance matrices
   PNTginvATNP<- PNTginvATNP[sort(1:length(PNTginvATNP), decreasing=TRUE)]
      
  #Now construct variance matrices
  if(all(is.na(var.comp))){ 
    V <- lapply(Z, function(x) x %*% t(x)) 
  }else{
    V <- lapply(makeBlkDesMatrix(design.df, var.comp), function(x) x %*% t(x))
  }
  
  
  VC <-character(length(V))
   names(VC)  = names(V)
   VC = t(data.frame(VC))
 ##############################################################################################################  
write("5. Get coefficients of each source of variation for the random effects,", "")
   
    for(i in 1:length(PNTginvATNP)){
      tmp <- matrix(0, nrow=length(names(PNTginvATNP[[i]])), ncol=length(V), 
                    dimnames=list(names(PNTginvATNP[[i]]), names(V)))        
      for(j in 1:(length(names(PNTginvATNP[[i]])))){
        for(z in 1:(length(V))){
          tmp[j,z] <- tr(PNTginvATNP[[i]][[j]] %*% V[[z]])                                   
        }         
      } 
      
      
       if(nrow(tmp) == 1 && rownames(tmp) == "Residual"){
         tmp = c(tmp[1], tmp[-1]/tmp[1])  
         VC = rbind(VC, attr(fractions(tmp),"fracs"))
          
         rownames(VC)[nrow(VC)] = paste("Between", names(PNTginvATNP[i]), sep = " ")
 
       }else{
        VC = rbind(VC, character(length = length(V)))
          
        if(names(PNTginvATNP[i]) == "Within"){
          rownames(VC)[nrow(VC)] =  paste(names(PNTginvATNP[i]), sep = " ")
        }else{
          rownames(VC)[nrow(VC)] = paste("Between", names(PNTginvATNP[i]), sep = " ")
        }
        
        rownames(tmp) = paste("  ", rownames(tmp), sep = " ")
        
        tmp = t(apply(tmp, 1, function(x) attr(fractions(c(x[1], x[-1]/x[1])),"fracs")))
        VC = rbind(VC, tmp)
      }           
      
    }
   colnames(VC)[1] = "DF"
   VC = noquote(VC[-1,])   
    
   if(table.legend){
    Legend = paste(paste( letters[1:(length(colnames(VC))-1)],colnames(VC)[-1], sep = " = "))
    colnames(VC)[-1] = letters[1:(length(colnames(VC))-1)]
    VC = list(VC, Legend)
   }


##############################################################################################################  
write("   and for fixed effects.", "")

  effFactors = lapply(Pb, function(z) trtProkMat(z, T, N, Rep))
  effFactors <- effFactors[sort(1:length(effFactors), decreasing=TRUE)]

   trt  = numeric(length(trt.Coef) + length(Rep))
   names(trt) = c( names(T), paste("eff",  names(T), sep = ".") )
 
   for(i in 1:length(effFactors)){
     trt = rbind(trt, character(length = length(T)*2))
      if( names(effFactors[i]) =="Within"){
        rownames(trt)[nrow(trt)] = paste(names(effFactors[i]), sep = " ")
      }else{
        rownames(trt)[nrow(trt)] = paste("Between", names(effFactors[i]), sep = " ")
      }            
           
      for(j in 1:length(effFactors[[i]])){
        if(is.null(effFactors[[i]][[j]])) next
        trt.temp = attr(fractions(c(trt.Coef*unlist(effFactors[[i]][[j]]),unlist(effFactors[[i]][[j]]))),"fracs")
        trt = rbind(trt,trt.temp)         
        rownames(trt)[nrow(trt)] =  paste("  ", names(effFactors[[i]][j]), sep = " ")  
      }
    
   }
    
    trt = trt[-1,]
  
    trt = noquote(ifelse(trt == "NaN", "", trt))
    
    if(table.legend){
      Legend = paste(paste( letters[1:(length(colnames(trt)))],colnames(trt), sep = " = "))
      colnames(trt) = letters[1:(length(colnames(trt)))]
      trt = list(trt = trt, Legend = Legend)
    }


   write(paste(rep("#", 80), sep ="", collapse = ""), "")

   return(list(random = VC, fixed = trt))
}

