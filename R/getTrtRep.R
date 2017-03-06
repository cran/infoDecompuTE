
getTrtRep <- function(design.df, trtTerm) {
    
    if (length(trtTerm) == 1 && !any(grepl("[[:punct:]]", trtTerm))) {
        
        return(list(Rep = as.matrix(table(design.df[, trtTerm])), Sca = 1))
    } else if (any(grepl("[[:punct:]]", trtTerm))) {
        
        level <- t(sapply(strsplit(sort(levels(interaction(design.df[, unique(unlist(strsplit(trtTerm, 
            "[[:punct:]]+")))]))), "\\."), rbind))
         	 	 
		 		
		 #level <- t(sapply(strsplit( sort(unique(apply(design.df[, unique(unlist(strsplit(trtTerm, #"[[:punct:]]")))],1,  function(x) paste(x, collapse = ".")))), "\\."), rbind))
  
        colnames(level) <- unique(unlist(strsplit(trtTerm, "[[:punct:]]+")))
        
        inter <- trtTerm[grepl("[[:punct:]]", trtTerm)]
        
        for (i in 1:length(inter)) {
            level <- cbind(level, apply(level[, unique(unlist(strsplit(inter[i], "[[:punct:]]+")))], 
                1, function(x) paste(x, collapse = ".")))
            colnames(level)[ncol(level)] <- inter[i]
        }
        
        trtTermList <- lapply(strsplit(trtTerm, "[[:punct:]]+"), function(x) design.df[, 
            x])
        names(trtTermList) <- trtTerm
        
        repList <- lapply(trtTermList, function(y) if (is.factor(y)) {
            table(y)
        } else {
            table(apply(y, 1, function(x) paste(x, collapse = ".")))
        })
        
        repMat <- level
        
        for (i in 1:length(repList)) {
            level.temp <- level[, names(repList)[i]]
            repMat <- cbind(repMat, repList[[i]][level.temp])
        }
        
        repMat <- repMat[, -(1:ncol(level))]
        
        if (is.matrix(repMat)) {
            repMat <- apply(repMat, 2, function(x) ifelse(is.na(x), 0, as.numeric(x)))
            
            colnames(repMat) <- names(repList)
            rownames(repMat) <- NULL
            
            
            levelList <- sapply(trtTermList, function(y) if (is.factor(y)) {
                nlevels(y)
            } else {
                nlevels(as.factor(apply(y, 1, function(x) paste(x, collapse = "."))))
            })/apply(repMat, 2, function(x) sum(x != 0))
            
            repList <- repMat %*% diag(levelList)
            
        } else {
            repMat <- ifelse(is.na(repMat), 0, as.numeric(repMat))
            
            levelList <- sapply(trtTermList, function(y) if (is.factor(y)) {
                nlevels(y)
            } else {
                nlevels(as.factor(apply(y, 1, function(x) paste(x, collapse = "."))))
            })/sum(repMat != 0)
            repList <- as.matrix(repMat * levelList)
        }
        
        		
        return(list(Rep = repList, Sca = levelList))
        
    } else {
        
        level <- t(sapply(strsplit(sort(levels(interaction(design.df[, trtTerm]))), 
            "\\."), rbind))
		
		#level <- t(sapply(strsplit( sort(unique(apply(design.df[, trtTerm],1,  
		#		function(x) paste(x, collapse = ".")))), "\\."), rbind))		
			
        colnames(level) <- trtTerm
        
        repList <- lapply(design.df[, trtTerm], table)
        
        repMat <- level
        
        for (i in 1:length(repList)) {
            level.temp <- level[, names(repList)[i]]
            repMat <- cbind(repMat, repList[[i]][level.temp])
        }
        
        repMat <- repMat[, -(1:ncol(level))]
        repMat <- apply(repMat, 2, function(x) ifelse(is.na(x), 0, as.numeric(x)))
        colnames(repMat) <- names(repList)
        rownames(repMat) <- NULL
        
        levelList <- apply(design.df[, trtTerm], 2, function(x) nlevels(as.factor(x)))/apply(repMat, 
            2, function(x) sum(x != 0))
        
        repList <- repMat %*% diag(levelList)
        
		
        return(list(Rep = repList, Sca = levelList))
        
    }
} 
