makeContrMat <- function(design.df, effectNames, effectsMatrix, contr.vec) {
    # Square contrast matrix with the specific contrast (trtContr)
    indMatrix1 <- function(x, n, trtContr) {
        if (x == 1) 
            X <- trtContr else if (x == 2) 
            X <- identityMat(n) else X <- K(n)
        return(X)
    }
    
    # Square contrast matrix without defining the specific contrast
    indMatrix <- function(x, n) {
        if (x == 1) 
            X <- J(n) else if (x == 2) 
            X <- identityMat(n) else X <- K(n)
        return(X)
    }
    
    # transform each treatment contrast matrix to C matrix
    transContrToT <- function(fact, contr) {
        T <- matrix(0, nrow = nlevels(fact), ncol = nlevels(fact))
        if (is.matrix(contr)) {
            rownames(contr) <- fact
            for (i in 1:ncol(contr)) {
                x <- contr[levels(fact), i]
                T <- T + (x %*% t(x))/as.numeric(t(x) %*% x)
            }
            
        } else {
            names(contr) <- fact
            x <- contr[levels(fact)]
            T <- T + (x %*% t(x))/as.numeric(t(x) %*% x)
        }
        return(T)
    }
    
	#browser()
	
    # obtaining the numbers of the levels from the design
    if (any(grepl("[[:punct:]]", effectNames))) {
        uniqueTrtCols <- unique(unlist(strsplit(effectNames, "[[:punct:]]+")))
        #browser()
        nLevels <- sapply(design.df[, uniqueTrtCols], function(x) nlevels(as.factor(x)))
        
    } else if (length(effectNames) == 1) {
        uniqueTrtCols = effectNames
        nLevels <- nlevels(design.df[, effectNames])
        names(nLevels) <- effectNames
        
    } else {
        uniqueTrtCols = effectNames
        nLevels <- sapply(design.df[, effectNames], function(x) nlevels(as.factor(x)))
    }
    
    nLevels <- nLevels[rownames(effectsMatrix)]
    nEffects <- ncol(effectsMatrix)
    
    
    # Without the contrasts specifically defined
    if (all(is.na(contr.vec))) {
        X <- as.list(rep(1, nEffects))
        names(X) <- effectNames
        
        for (i in 1:nrow(effectsMatrix)) {
            matList <- lapply(effectsMatrix[i, ], function(y) indMatrix(y, nLevels[i]))
            for (j in 1:nEffects) X[[j]] <- X[[j]] %x% matList[[j]]
        }
        
    } else {
        new.contr.vec <- vector(length = length(uniqueTrtCols), mode = "list")
        names(new.contr.vec) <- uniqueTrtCols
   		
        for (i in 1:length(new.contr.vec)) {
            new.contr.vec[i] <- ifelse(is.null(contr.vec[[names(new.contr.vec)[i]]]), 
                NA, contr.vec[names(new.contr.vec)[i]])
        }
        
                 
        X <- list()
        count <- 1
        
        contr.vec <- new.contr.vec
        #browser()
		
		effectList =  vector(length = length(effectNames), mode = "list")
		names(effectList) = effectNames
		
		 # Construct the interactions based on the contrasts
        #if (any(grepl("[[:punct:]]", effectNames))) {
        #    interTerms <- grep("[[:punct:]]", effectNames)
            
            for (i in effectNames) {
                mainTerms <- unlist(strsplit(i, "[[:punct:]]"))
                mainTerms.list <- vector(length = length(mainTerms), mode = "list")
                names(mainTerms.list) <- mainTerms
                
                for (j in 1:length(mainTerms)) {
					if (is.null(names(contr.vec[[mainTerms[j]]]))) 
						mainTerms.list[[j]] <- 0 
					else 
						mainTerms.list[[j]] <- names(contr.vec[[mainTerms[j]]])
                }
                
                
                effectList[[i]] <- vector(length = nlevels(interaction(mainTerms.list)), 
                  mode = "list")
                # names(contr.vec)[[i]] = effectNames[i]
                
                names(effectList[[i]]) <- levels(interaction(mainTerms.list))
            }
        #}
	
		totalLength <- length(effectList[sapply(effectList, class) != "list"])
		
		if(any(sapply(effectList, class) == "list")){
	
			 totalLength = totalLength + sum(sapply(effectList[which(sapply(effectList, class) == "list")], length))
		}
        
        newEffectsMatrix <- matrix(0, ncol = totalLength, nrow = nrow(effectsMatrix))
        rownames(newEffectsMatrix) <- rownames(effectsMatrix)
   
		
        # construct a new effects matrix to use for the Kronecker products
        for (i in 1:length(effectList)) {
             
            if (is.list(effectList[[i]])) {
                tmpX <- rep(1, length(effectList[[i]]))
                names(tmpX) <- paste(names(effectList[i]), names(effectList[[i]]), sep = ".")
                
                tmpCounter <- count:(count + length(effectList[[i]]) - 1)
                
                X <- c(X, tmpX)
                rowNames <- unlist(strsplit(names(effectList[i]), "[[:punct:]]"))
                newEffectsMatrix[rowNames, tmpCounter] <- 1
                count <- count + length(effectList[[i]])
                
            } else {
                tmpX <- 1
                names(tmpX) <- names(effectList[i])
                X <- c(X, tmpX)
                
                newEffectsMatrix[, count] <- effectsMatrix[, i]
                
                tmpCounter <- count <- count + 1
            }
        }
        
        colnames(newEffectsMatrix) <- names(X)
        #browser()
		
        for (i in 1:nrow(newEffectsMatrix)) {
            if (is.list(contr.vec[[rownames(effectsMatrix)[i]]])) {
                trtContr <- lapply(contr.vec[[rownames(effectsMatrix)[i]]], function(x) transContrToT(design.df[, 
                  rownames(effectsMatrix)[i]], x))
                
                matList <- vector(mode = "list", length = totalLength)
                
                for (j in 1:ncol(newEffectsMatrix)) {
                 
                  if (newEffectsMatrix[i, j] == 1) {
                    trtNames <- match(unlist(strsplit(colnames(newEffectsMatrix)[j], 
                      "\\.")), names(trtContr))
                    
                    matList[[j]] <- trtContr[[trtNames[!(is.na(trtNames))]]]
                  } else if (newEffectsMatrix[i, j] == 2) {
                    matList[[j]] <- identityMat(nLevels[i])
                  } else {
                    matList[[j]] <- K(nLevels[i])
                  }
                }
            } else {
                if (all(is.na(contr.vec[[i]]))) {
                  matList <- lapply(newEffectsMatrix[i, ], function(y) indMatrix(y, 
                    nLevels[i]))
                  
                } else {
                  trtContr <- transContrToT(design.df[, rownames(effectsMatrix)[i]], 
                    contr.vec[[rownames(effectsMatrix)[i]]])
                  matList <- lapply(newEffectsMatrix[i, ], function(y) indMatrix1(y, 
                    nLevels[i], trtContr))
                  
                }
            }
               
            for (j in 1:length(X)) X[[j]] <- X[[j]] %x% matList[[j]]
            
        }
        
        
        
    }
    
	
	 
    names(X) <- gsub("\\.0", "", names(X))
    return(X)
} 
