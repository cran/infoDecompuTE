summaryAovTwoPhase <- function(design.df, blk.str1, blk.str2, trt.str, var.comp = NA, 
    blk.contr = NA, trt.contr = NA, table.legend = FALSE, response = NA, latex = FALSE, 
    fixed.names = NA, decimal = FALSE, digits = 2, list.sep = TRUE) {
    
  design.df <- data.frame(sapply(design.df, function(x) gsub("[[:punct:]]", "", as.character(x))))
    
	#browser()
	newTerms = adjustMissingLevels(design.df, trt.str)
		
	design.df = newTerms$design.df
	trt.str = newTerms$str.for
	
	newTerms = adjustMissingLevels(design.df, blk.str1)
		
	design.df = newTerms$design.df
	blk.str1 = newTerms$str.for
	
    # Extract the fixed and random terms
    
    rT1 <- stats::terms(stats::as.formula(paste("~", blk.str1, sep = "")), keep.order = TRUE)  #random terms phase 1
    rT2 <- stats::terms(stats::as.formula(paste("~", blk.str2, sep = "")), keep.order = TRUE)  #random terms phase 2
    fT <- stats::terms(stats::as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  #fixed terms
    
    
    # Preparing the block structures write('1.  Preparing the block structure.', '')
    
    blkTerm1 <- attr(rT1, "term.labels")
    blkTerm2 <- attr(rT2, "term.labels")
    
    # check for complete confounding and changing the names of the block factors
    # browser()
    if (any(grepl("\\:", blkTerm2))) {
        
        check.blkTerm2 <- unique(unlist(strsplit(blkTerm2, "\\:")))
    } else {
        check.blkTerm2 <- blkTerm2
    }
    
    if (any(grepl("\\:", blkTerm1))) {
        check.blkTerm1 <- unique(unlist(strsplit(blkTerm1, "\\:")))
    } else {
        check.blkTerm1 <- blkTerm1
    }
    
    # Check for Complete confounding.
    for (i in 1:length(check.blkTerm2)) {
        if (length(check.blkTerm1) == 1) {
            if (all(as.numeric(as.factor(design.df[, check.blkTerm1])) == as.numeric(as.factor(design.df[, 
                check.blkTerm2[i]])))) {
                cat("Note: Complete confounding between ", check.blkTerm1, " and ", 
                  check.blkTerm2[i], "!\n", sep = "")
                
                colnames(design.df)[which(colnames(design.df) == check.blkTerm1)] <- paste(colnames(design.df)[which(colnames(design.df) == 
                  check.blkTerm1)], "CCW", sep = "")
                
                blk.str1 <- paste(blkTerm1[which(blkTerm1 == check.blkTerm1)], "CCW", 
                  sep = "")
                
                if (!is.na(blk.contr)) {
                  names(blk.contr)[which(names(blk.contr) == check.blkTerm1)] <- paste(names(blk.contr)[which(names(blk.contr) == 
                    check.blkTerm1)], "CCW", sep = "")
                }
                check.blkTerm1 <- paste(check.blkTerm1, "CCW", sep = "")
                
            }
        } else {
            check.temp <- apply(design.df[, check.blkTerm1], 2, function(x) all(as.numeric(as.factor(x)) == 
                as.numeric(as.factor(design.df[, check.blkTerm2[i]]))))
            if (any(check.temp)) {
                cat("Note: Complete confounding between ", check.blkTerm1[which(check.temp)], 
                  " and ", check.blkTerm2[i], "!\n", sep = "")
                colnames(design.df)[which(colnames(design.df) == check.blkTerm1[which(check.temp)])] <- paste(colnames(design.df)[which(colnames(design.df) == 
                  check.blkTerm1[which(check.temp)])], "CCW", sep = "")
                
                blk.str1 <- gsub(check.blkTerm1[which(check.temp)], paste(check.blkTerm1[which(check.temp)], 
                  "CCW", sep = ""), blk.str1)
                if (!is.na(blk.contr)) {
                  names(blk.contr)[which(names(blk.contr) == check.blkTerm1[which(check.temp)])] <- paste(names(blk.contr)[which(names(blk.contr) == 
                    check.blkTerm1[which(check.temp)])], "CCW", sep = "")
                }
                check.blkTerm1[which(check.blkTerm1 == check.blkTerm1[which(check.temp)])] <- paste(check.blkTerm1[which(check.blkTerm1 == 
                  check.blkTerm1[which(check.temp)])], "CCW", sep = "")
            }
        }
    }
    
	rT1 <- stats::terms(stats::as.formula(paste("~", blk.str1, sep = "")), keep.order = TRUE)  #random terms phase 1
	blkTerm1 <- attr(rT1, "term.labels")
     
	
	#browser()
	
    Z1 <- makeBlkDesMat(design.df, rev(blkTerm1))
    Z2 <- makeBlkDesMat(design.df, rev(blkTerm2))
    
    # write('2. Defining the block structures of second Phase.', '')
    Pb <- makeOrthProjectors(Z2)
    
	if(length(blkTerm2) > 1)
		blkTerm2 = adjustEffectNames(effectsMatrix = attr(rT2, "factors"), effectNames = blkTerm2)

    
	
    if (names(Pb)[1] == "e") {
        names(Pb)[1] <- paste("Within", paste(unique(unlist(strsplit(names(Pb)[-1], 
            "[[:punct:]]+"))), collapse = "."))
        
        names(Pb)[-1] <- names(Z2)[-1] <-rev(blkTerm2)
		 
    } else {
        names(Pb) <- names(Z2) <- rev(blkTerm2)
    }
    
    # write('3. Defining the block structures of first phase within second Phase.', '')
    effectsMatrix <- attr(rT1, "factor")
    
	effectsMatrix[nrow(effectsMatrix), effectsMatrix[nrow(effectsMatrix),]==2] <- 1
	
	if(length(blkTerm1) > 1)
		blkTerm1 = adjustEffectNames(effectsMatrix = effectsMatrix, effectNames = blkTerm1)
	
	if (names(Z1)[1] == "e") {
		names(Z1)[-1] <- rev(blkTerm1)
		 
    } else {
        names(Z1) <- rev(blkTerm1)
    }
 	
	#browser()   
    T <- makeContrMat(design.df, effectNames = blkTerm1, effectsMatrix = effectsMatrix, 
        contr.vec = blk.contr)
    N <- makeOverDesMat(design.df, blkTerm1)
    

    
    #browser()
    
    res <- paste("Within", paste(unique(unlist(strsplit(names(T), "[[:punct:]]+"))), 
        collapse = "."))
    
    Pb1 <- lapply(Pb, function(z) infoDecompMat(z, T, N))
    
    # t.name = unique(unlist(strsplit(names(T), '[[:punct:]]')))
    t.name <- names(T)
    
    if (length(Pb) > 1) {
        pb1.names <- lapply((Pb1[-1]), names)
    } else {
        pb1.names <- lapply(Pb1, names)
    }
    
    # browser()
    for (i in 1:length(pb1.names)) {
        comp <- t.name %in% pb1.names[[i]]
        if (any(comp)) {
            break
        } else if (i == length(pb1.names)) {
            names(Pb1)[1] <- ""
        }
    }
    
    for (i in 1:length(Pb1)) {
        names(Pb1[[i]])[which(names(Pb1[[i]]) == "Residual")] <- res
    }
    
    # Preparing the treatment structures write('4.  Preparing the treatment structure.',
    # '')
    
   trtTerm <- attr(fT, "term.labels")
    effectsMatrix <- attr(fT, "factor")
     #browser()
    
	 effectsMatrix[nrow(effectsMatrix), effectsMatrix[nrow(effectsMatrix),]==2] <- 1
	 
	#browser()
	if(length(trtTerm) > 1)
		trtTerm = adjustEffectNames(effectsMatrix, trtTerm)
		
  	T <- makeContrMat(design.df = design.df, effectNames = trtTerm, effectsMatrix = effectsMatrix, 
  	                  contr.vec = trt.contr)
  	N <- makeOverDesMat(design.df = design.df, effectNames = trtTerm)
  	Replist <- getTrtRep(design.df, trtTerm)
	Rep <- Replist$Rep
  	trt.Sca <- Replist$Sca
  	
	
	#When there are treatment contrasts defined 
  	if (any(grepl("\\.", names(T)))) {
  	  colnames(Rep) <- trtTerm
  	  names(trt.Sca) <- trtTerm 
  	  Rep <- Rep[, sapply(strsplit(names(T), "\\."), function(x) x[1])]
  	  trt.Sca <- trt.Sca[sapply(strsplit(names(T), "\\."), function(x) x[1])]
  	} else {
  	  colnames(Rep) <- trtTerm
  	  names(trt.Sca) <- trtTerm      
	}
	 
    
	#browser()
	
    # Start calculating the VCs 2-phase experiment write('5. Start calculating the
    # variance components.', '') write('6. Pre- and post-multiply NTginvATN by block
    # projection matrices.', '')
    
	effFactors <- lapply(Pb1, function(y) lapply(y, function(z) getEffFactor(z, T, N, 
	                                                                         Rep, trt.Sca)))
	effFactors <- effFactors[sort(1:length(effFactors), decreasing = TRUE)]   
    
    
    v.mat <- getVMat.twoPhase(Z.Phase1 = Z1, Z.Phase2 = Z2, design.df = design.df, var.comp = var.comp)
    
    if (all(is.na(var.comp))) {
        names(v.mat)[-1] <- c(rev(blkTerm1), rev(blkTerm2))
    }
    
    
    ANOVA <- getCoefVC.twoPhase(Pb = effFactors, design.df = design.df, v.mat = v.mat, 
        response, table.legend, decimal = decimal, digits = digits)
    
    
    effFactors <- lapply(Pb1, function(y) lapply(y, function(z) getEffFactor(z, T, N, 
        Rep, trt.Sca)))
    effFactors <- effFactors[sort(1:length(effFactors), decreasing = TRUE)]
    
    
    if (latex) {
        EF <- getFixedEF.twoPhase(effFactors = effFactors, trt.Sca = trt.Sca, T = T, Rep = Rep, 
                                table.legend, decimal = decimal, digits = digits, list.sep = FALSE)
        return(toLatexTable(ANOVA = ANOVA, EF = EF, fixed.names = fixed.names))
    } else {
      EF <- getFixedEF.twoPhase(effFactors = effFactors, trt.Sca = trt.Sca, T = T, Rep = Rep, 
                                table.legend, decimal = decimal, digits = digits, list.sep = list.sep)
        return(list(ANOVA = ANOVA, Fixed = EF))
    }
    
    
} 
