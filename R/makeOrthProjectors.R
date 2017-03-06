
makeOrthProjectors <- function(BlkDesList) {
    
    initial <- diag(nrow(BlkDesList[[1]])) - K(nrow(BlkDesList[[1]]))
    n <- nrow(BlkDesList$e)
    Q <- lapply(BlkDesList, function(z) projMat(z))
    
    Q <- Q[sort(1:length(Q), decreasing = TRUE)]
    
    elementToRemove <- numeric(0)
    
    cusumQ <- P <- NULL
    P[[1]] <- Q[[1]] %*% initial
    
    if (all(P[[1]] < 1e-06)) {
        elementToRemove <- 1  #revord the elements that are all zeros
        P[[1]] <- matrix(0, nrow = n, ncol = n)  #make the projectors which has less than 1e-6 to zero, to avoid rounding error
    }
    
    cusumQ[[1]] <- initial - P[[1]]
    
    for (i in 2:(length(Q))) {
        P[[i]] <- Q[[i]] %*% cusumQ[[i - 1]]
        
        if ((all(P[[i]] < 1e-06)) || tr(P[[i]]) <= 0) {
            elementToRemove <- c(elementToRemove, i)  #revord the elements that are all zeros
            P[[i]] <- matrix(0, nrow = n, ncol = n)  #make the projectors which has less than 1e-6 to zero, to avoid rounding error
        }
        cusumQ[[i]] <- cusumQ[[i - 1]] - P[[i]]
    }
    
    P <- P[sort(1:length(P), decreasing = TRUE)]
    names(P) <- names(BlkDesList)
    
    if (length(elementToRemove) > 0) 
        P <- P[-(length(Q) - elementToRemove + 1)]
    
    return(P)
} 
