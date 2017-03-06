getCoefVC.twoPhase <- function(Pb, design.df, v.mat, response, table.legend, decimal, 
    digits) {
    
    if (all(is.na(response))) {
        response <- rep(NA, nrow(design.df))
    }
    
    MS <- lapply(Pb, function(x) lapply(x, function(q) lapply(q[[1]], function(w) t(response) %*% 
        w %*% response)))
    
    V <- v.mat
    VC <- rep("1", length(V) + 2)
    names(VC) <- c("DF", c(names(V), "MS"))
    VC <- t(data.frame(VC))
    noindent <- FALSE
    # browser()
    for (i in 1:length(Pb)) {
        VC <- rbind(VC, character(length = length(V) + 2))
        
        if (grepl("Within", names(Pb[i]))) {
            rownames(VC)[nrow(VC)] <- paste(names(Pb[i]), sep = " ")
        } else if (names(Pb[i]) == "") {
            VC <- VC[-(nrow(VC)), ]
            noindent <- TRUE
        } else {
            rownames(VC)[nrow(VC)] <- paste("Between", names(Pb[i]), sep = " ")
        }
        
        for (k in 1:(length(Pb[[i]]))) {
            tmp <- matrix(0, nrow = length(names(Pb[[i]][[k]][[1]])), ncol = (length(V) + 
                1), dimnames = list(names(Pb[[i]][[k]][[1]]), c(names(V), "MS")))
            
            for (j in 1:((length(names(Pb[[i]][[k]][[1]]))))) {
                for (z in 1:length(V)) {
                  tmp[j, z] <- tr(Pb[[i]][[k]][[1]][[j]] %*% V[[z]])
                }
                tmp[j, (length(V) + 1)] <- as.numeric(MS[[i]][[k]][[j]])
            }
            
            if (nrow(tmp) == 1 && grepl("Residual", rownames(tmp))) {
                

                if (decimal) {
                  tmp <- t(apply(tmp, 1, function(x) c(round(c(x[1], x[-length(x)]/x[1]), 
                                                             digits = digits), as.character(round(x[length(x)]/x[1], digits = 3)))))
                } else {
                  tmp <- t(apply(tmp, 1, function(x) c(attr(fractions(c(round(x[1]), x[-length(x)]/x[1])), 
                                                            "fracs"), as.character(round(x[length(x)]/x[1], digits = 3)))))
                }
                
                VC <- rbind(VC, tmp)
                if (grepl("Within", names(Pb[[i]][k])) && length(Pb[[i]]) > 1) {
                  if (noindent) {
                    rownames(VC)[nrow(VC)] <- names(Pb[[i]][k])
                  } else {
                    rownames(VC)[nrow(VC)] <- paste("  ", names(Pb[[i]][k]), sep = " ")
                  }
                } else if (grepl("Within", names(Pb[[i]][k])) && length(Pb[[i]]) == 1) {
                  rownames(VC)[nrow(VC)] <- paste("Between", names(Pb[i]), sep = " ")
                  VC <- VC[-(nrow(VC) - 1), ]
                } else {
                  if (noindent) {
                    rownames(VC)[nrow(VC)] <- paste("Between", names(Pb[[i]][k]), sep = " ")
                  } else {
                    rownames(VC)[nrow(VC)] <- paste("   Between", names(Pb[[i]][k]), 
                      sep = " ")
                  }
                }
            } else {
                
                VC <- rbind(VC, character(length = length(V) + 2))
                
                if (grepl("Within", names(Pb[[i]][k])) && length(Pb[[i]]) == 1) {
                  
                  noindent <- TRUE
                }
                
                if (noindent) {
                  if (grepl("Within", names(Pb[[i]][k])) && length(Pb[[i]]) == 1) {
                    VC <- VC[-(nrow(VC)), ]
                  } else if (grepl("Within", names(Pb[[i]][k]))) {
                    rownames(VC)[nrow(VC)] <- names(Pb[[i]][k])
                  } else {
                    rownames(VC)[nrow(VC)] <- paste("Between", names(Pb[[i]][k]), sep = " ")
                  }
                  rownames(tmp) <- paste("  ", rownames(tmp), sep = " ")
                } else {
                  if (grepl("Within", names(Pb[[i]][k]))) {
                    rownames(VC)[nrow(VC)] <- paste("  ", names(Pb[[i]][k]), sep = " ")
                  } else {
                    rownames(VC)[nrow(VC)] <- paste("   Between", names(Pb[[i]][k]), 
                      sep = " ")
                  }
                  rownames(tmp) <- paste("     ", rownames(tmp), sep = " ")
                }
                
                
                
                if (decimal) {
                  tmp <- t(apply(tmp, 1, function(x) c(round(c(x[1], x[-length(x)]/x[1]), 
                    digits = digits), as.character(round(x[length(x)]/x[1], digits = 3)))))
                } else {
                  tmp <- t(apply(tmp, 1, function(x) c(attr(fractions(c(round(x[1]), x[-length(x)]/x[1])), 
                    "fracs"), as.character(round(x[length(x)]/x[1], digits = 3)))))
                }
                
                
                VC <- rbind(VC, tmp)
            }
        }
        noindent <- FALSE
    }
    
    if (all(is.na(response))) {
        VC <- VC[, -ncol(VC)]
    }
    
    if (length((which(apply(apply(VC[, -1], 2, function(x) x == VC[, "e"]), 2, all)))) > 
        1) {
        VC <- noquote(VC[-1, -2])
    } else {
        VC <- noquote(VC[-1, ])
    }
    
    if (table.legend) {
        Legend <- paste(paste(letters[1:(length(colnames(VC)) - 1)], colnames(VC)[-1], 
            sep = " = "))
        colnames(VC)[-1] <- letters[1:(length(colnames(VC)) - 1)]
        VC <- list(VC = VC, Legend = Legend)
    }
    
    return(VC)
} 
