
getTrtCoef <- function(design.df, trtTerm) {
    
    if (length(trtTerm) == 1 && !any(grepl("[[:punct:]]", trtTerm))) {
        return(length(design.df[, trtTerm])/nlevels(design.df[, trtTerm]))
    } else if (any(grepl("[[:punct:]]", trtTerm))) {
        
        trtTermList <- lapply(strsplit(trtTerm, "[[:punct:]]+"), function(x) design.df[, 
            x])
        repList <- sapply(trtTermList, function(y) if (is.factor(y)) {
            mean(table(y))
        } else {
            mean(table(apply(y, 1, function(x) paste(x, collapse = "."))))
        })
        
        return(repList)
    } else {
        repList <- apply(design.df[, trtTerm], 2, function(x) mean(table(x)))
        return(repList)
    }
} 
