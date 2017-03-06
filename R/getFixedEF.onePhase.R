getFixedEF.onePhase <- function(effFactors, trt.Sca, T, Rep, table.legend, decimal, 
    digits, list.sep) {
    
    trt <- numeric(length(trt.Sca) + ncol(Rep))
    names(trt) <- c(names(T), paste("eff", names(T), sep = "."))
    
    for (i in 1:length(effFactors)) {
        trt <- rbind(trt, character(length = length(T) * 2))
        if (grepl("Within", names(effFactors[i]))) {
            rownames(trt)[nrow(trt)] <- paste(names(effFactors[i]), sep = " ")
        } else {
            rownames(trt)[nrow(trt)] <- paste("Between", names(effFactors[i]), sep = " ")
        }
        
        for (j in 1:length(effFactors[[i]][[2]])) {
            if (is.null(effFactors[[i]][[2]][[j]])) 
                next
            
					
			effCoefList = effFactors[[i]][[2]][[j]]
				
			#browser()
				
			if (decimal) {
				char.trt.eff <- round(sapply(effCoefList, function(x) x[1]), digits = digits)
				char.trt <- sapply(effCoefList, function(x) 
					ifelse(length(x)>2,  paste(round(x[2:length(x)], digits = digits), collapse = ","), 
									round(x[2], digits = digits)))
			} else {
				char.trt <- sapply(effCoefList, function(x) 
					ifelse(length(x)>2,  paste(attr(fractions(x[2:length(x)]), "fracs"), collapse = ","), 
										attr(fractions(x[2]), "fracs")))
									
				char.trt.eff <- attr(fractions(sapply(effCoefList, function(x) x[1])), "fracs")
			}
			
            trt.temp <- c(char.trt, char.trt.eff)
            
            trt <- rbind(trt, trt.temp)
            rownames(trt)[nrow(trt)] <- paste("  ", names(effFactors[[i]][[2]][j]), sep = " ")
        }
        
    }
    
    trt <- trt[-1, ]
    
    trt <- noquote(ifelse(trt == "NaN", "", trt))
    trt <- noquote(ifelse(trt == "0", "", trt))
	  
    if(list.sep){
      trt.Fixed = trt[,-grep("^eff", colnames(trt))]
      trt.EF = trt[,grep("^eff", colnames(trt))]
      
      if(length(grep("^eff", colnames(trt))) ==1){
        trt.Fixed = noquote(matrix(trt.Fixed))
        trt.EF = noquote(matrix(trt.EF))
        
        rownames(trt.Fixed) = rownames(trt.EF) = rownames(trt)
        colnames(trt.Fixed) =  colnames(trt)[1]
        colnames(trt.EF) = paste("eff.", colnames(trt)[1], sep = "")
      }
      trt = list(EF = trt.EF, Coef = trt.Fixed)
    }
    
    
    if (table.legend) {
      if(list.sep){
        Legend.EF <- paste(paste(letters[1:(length(colnames(trt.EF)))], colnames(trt.EF), sep = " = "))
        colnames(trt.EF) <- letters[1:(length(colnames(trt.EF)))]
        Legend.Coef <- paste(paste(letters[1:(length(colnames(trt.Fixed)))], colnames(trt.Fixed), sep = " = "))
        colnames(trt.Fixed) <- letters[1:(length(colnames(trt.Fixed)))]
        
        trt <- list(EF = trt.EF, Legend.EF = Legend.EF, Coef = trt.Fixed, Legend.Coef = Legend.Coef)
      } else{
        Legend <- paste(paste(letters[1:(length(colnames(trt)))], colnames(trt), sep = " = "))
        colnames(trt) <- letters[1:(length(colnames(trt)))]
        trt <- list(trt = trt, Legend = Legend)
      }
    }
    
    return(trt)
    
} 
