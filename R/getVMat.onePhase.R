getVMat.onePhase <- function(Z.Phase1, design.df, var.comp = NA) {
    
    v.mat <- lapply(Z.Phase1, function(x) x %*% t(x))
    # Now construct variance matrices
    if (all(is.na(var.comp))) {
        return(v.mat)
        
    } else {
        
        if (names(v.mat)[1] == "e") {
            match.names <- match(var.comp, names(v.mat))
            
            if (any(is.na(match.names))) 
                match.names <- match.names[!is.na(match.names)]
            
            return(v.mat[c(1, match.names)])
        } else {
            match.names <- match(var.comp, names(v.mat))
            
            if (any(is.na(match.names))) 
                match.names <- match.names[!is.na(match.names)]
            
            return(v.mat[match.names])
        }
        
    }
    
} 
