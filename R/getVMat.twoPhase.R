# Now construct variance matrices
getVMat.twoPhase <- function(Z.Phase1, Z.Phase2, design.df, var.comp = NA) {
    #browser()
    v.mat <- lapply(c(Z.Phase1, Z.Phase2[-1]), function(x) x %*% t(x))
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
