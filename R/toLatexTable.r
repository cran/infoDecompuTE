
toLatexTable <- function(ANOVA, EF, fixed.names) {
    
    
    matchRowNames <- rownames(ANOVA)
    random.ColNames <- colnames(ANOVA)
    fixed.ColNames <- colnames(EF)
    
    #browser()
    
    if (is.na(fixed.names[1])) 
        fixed.names <- c("\\tau", "\\gamma", "\\rho", "\\phi", "\\delta", "\\omega")
    
    fixed.ColNames <- gsub("\\.", "", fixed.ColNames)
    
    match.fixed.names <- unique(unlist(strsplit(fixed.ColNames[1:(length(fixed.ColNames)/2)], 
        "[[:punct:]]")))
    fixed.names <- fixed.names[1:length(match.fixed.names)]
    match.fixed.names <- c(match.fixed.names, "*", "(", ")")
    fixed.names <- c(fixed.names, "*", "(", ")")
    
    fixed.names <- sapply(strsplit(gsub("([[:punct:]])", "\\.\\1\\.", fixed.ColNames[1:(length(fixed.ColNames)/2)]), 
        "\\."), function(x) paste(fixed.names[match(x, match.fixed.names)], collapse = ""))
    
    fixed.names <- gsub("\\*", "", fixed.names)
    
    # browser()
    
    if (length(grep("MS", random.ColNames)) == 1) {
        ColNamesLen <- length(random.ColNames) - 1
    } else {
        ColNamesLen <- length(random.ColNames)
    }
    
    
    if (any(random.ColNames == "e")) {
        startColNames <- 3
    } else {
        startColNames <- 2
    }
    
    
    # check for interaction effects
    random.names <- sapply(strsplit(gsub("([[:punct:]])", "\\.\\1\\.", random.ColNames[startColNames:ColNamesLen]), 
        "\\."), function(x) paste(substr(x, 1, 1), collapse = ""))
    
    random.names <- gsub("\\*", "", random.names)
    
    #random.names <- gsub("\\:", "", random.names)
    
    if (ColNamesLen > 2 && startColNames == 3) {
        random.names <- c("\\sigma^2", paste("\\sigma_{", random.names, "}^2", sep = ""))
    } else if (startColNames == 2) {
        random.names <- c(paste("\\sigma_{", random.names, "}^2", sep = ""))
    } else {
        random.names <- c("\\sigma^2")
    }
    
    ANOVA <- ifelse(ANOVA == "0", "", ANOVA)
    finalTable <- cbind(ANOVA, EF[match(matchRowNames, rownames(EF)), ])
    tempEF <- EF
    
    # avoid repeat in the fixed componenets <- check!!!!!
    for (i in 1:nrow(ANOVA)) {
        index <- which(matchRowNames[i] == rownames(tempEF))[1]
        
        if (is.na(index)) 
            next
        
        finalTable[i, (ncol(ANOVA) + 1):ncol(finalTable)] <- tempEF[index, ]
        if (nrow(tempEF) == 2) 
            next
        
        tempEF <- tempEF[-index, ]
        # rownames(EF)[grep(matchRowNames[i], rownames(EF))[1]] <- ''
    }
    
    finalTable <- ifelse(is.na(finalTable), "", finalTable)
    
    output <- "\\begin{table}[ht]\n\\centering\n \\caption{Theoretical ANOVA table}\n"
    
    if (length(grep("MS", random.ColNames)) == 1) {
        output <- c(output, paste("\\begin{tabular}[t]{lrll", paste(rep("l", length(fixed.names)), 
            collapse = ""), "} \n", sep = ""))
    } else {
        output <- c(output, paste("\\begin{tabular}[t]{lrl", paste(rep("l", length(fixed.names)), 
            collapse = ""), "} \n", sep = ""))
    }
    
    output <- c(output, "\\toprule \n")
    
    
    firstRow <- paste("\\multicolumn{1}{l}{\\textbf{Source of Variation}} & \\multicolumn{1}{l}{\\textbf{DF}} & \\multicolumn{1}{l}{\\textbf{EMS}}&", 
        sep = "")
    
    if (length(grep("MS", random.ColNames)) == 1) 
        firstRow <- rbind(firstRow, "\\multicolumn{1}{l}{\\textbf{MS}} &")
    
    firstRow <- rbind(firstRow, paste(paste("\\multicolumn{1}{l}{$\\bm{E_{", fixed.names, 
        sep = "", collapse = "}}$}&"), "}}$}\\\\ \n", sep = ""))
    
    # browser()
    
    output <- c(output, firstRow)
    output <- c(output, "\\midrule \n")
    for (i in 1:length(matchRowNames)) {
        
        # row names
        SV <- matchRowNames[i]
        SV <- gsub("   ", "\\\\quad ", SV)
        
        # DF
        DF <- paste("$", finalTable[i, 1], "$", sep = "")
        DF <- ifelse(DF == "$$", "", DF)
        
        # Random VC
        coef.VC <- finalTable[i, 2:(1 + length(random.names))]
        random.VC <- ""
        for (j in 1:length(coef.VC)) {
            if (coef.VC[j] == "") 
                next
            if (coef.VC[j] == 1) 
                coef.VC[j] <- ""
            
            random.VC <- paste(random.VC, coef.VC[j], random.names[j], sep = "")
            
            random.VC <- paste(random.VC, "+", sep = "")
        }
        
        # Fixed VC
        if (length(grep("MS", random.ColNames)) == 1) {
            coef.VC <- finalTable[i, (3 + length(random.names)):(2 + length(random.names) + 
                length(fixed.names))]
        } else {
            coef.VC <- finalTable[i, (2 + length(random.names)):(1 + length(random.names) + 
                length(fixed.names))]
        }
        
        fixed.VC <- ""
        for (j in 1:length(coef.VC)) {
            if (coef.VC[j] == "") 
                next
            
            if (coef.VC[j] == 1) 
                coef.VC[j] <- ""
            
			if( grepl(",", coef.VC[j])){
				fixed.VC <- paste(fixed.VC, "(", coef.VC[j], ")\\theta_{", fixed.names[j], "}", 
					sep = "")
			}else {			
				fixed.VC <- paste(fixed.VC, coef.VC[j], "\\theta_{", fixed.names[j], "}", 
					sep = "")
            }
			
            fixed.VC <- paste(fixed.VC, "+", sep = "")
        }
        
        total.VC <- ifelse(random.VC == "", "", paste("$", random.VC, "+", fixed.VC, 
            "$", sep = ""))
        total.VC <- ifelse(fixed.VC == "", paste("$", random.VC, "$", sep = ""), total.VC)
        total.VC <- gsub("\\+\\+", "\\+", total.VC)
        total.VC <- gsub("\\+\\$", "\\$", total.VC)
        total.VC <- gsub("\\$\\$", "", total.VC)
        
        if (length(grep("MS", random.ColNames)) == 1) {
            eff <- paste("$", finalTable[i, (3 + length(random.names) + length(fixed.names)):ncol(finalTable)], 
                "$", sep = "", collapse = " & ")
        } else {
            eff <- paste("$", finalTable[i, (2 + length(random.names) + length(fixed.names)):ncol(finalTable)], 
                "$", sep = "", collapse = " & ")
        }
        
        
        eff <- gsub("\\$\\$", "", eff)
        
        # browser()
        finalTable[i, 2:(1 + length(random.names))]
        
        if (length(grep("MS", random.ColNames)) == 1) {
            currentRow <- paste(SV, " & ", DF, " & ", total.VC, " & ", finalTable[i, 
                2 + length(random.names)], " &", eff, "\\\\", sep = "")
        } else {
            currentRow <- paste(SV, " & ", DF, " & ", total.VC, " &", eff, "\\\\", sep = "")
        }
        
        if ((grepl("Between", matchRowNames[i + 1]) || grepl("Within", matchRowNames[i + 
            1])) && !grepl("^Within", matchRowNames[i]) && !grepl("^Between", matchRowNames[i])) {
            currentRow <- paste(currentRow, "\\hline \n")
        } else {
            currentRow <- paste(currentRow, "\n")
        }
        
        output <- c(output, currentRow)
    }
    
    output <- c(output, "\\bottomrule \n \\end{tabular} \n \\label{tab:} \n\\end{table} \n")
    output <- c(output, "")
    
    return(cat(output))
} 
