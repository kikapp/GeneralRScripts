library(xtable)

sum.anno.single <- function(temp.df, varname, labelname, n) {
  anote.string <- paste0(temp.df[[labelname]][1], "\nMean: ", round(mean(temp.df[[varname]]),n), "\n",
                         "Median: ",  round(median(temp.df[[varname]]),n), "\n",
                         "SD: ", round(sd(temp.df[[varname]]),n), "\n",
                         "IQR: ", round(quantile(temp.df[[varname]], 0.75) - quantile(temp.df[[varname]], 0.25),n), "\n")
  return(anote.string)
}


# varname <- "log.lbm"
# # labelname <- "statin.use"
# temp.labelname <- "statin.use"
# groupname <- unique(temp.df[[labelname]])[1]
# temp.varname <- varname
# temp.df <- culled.data
# n <- 3
sum.anno <- function(temp.df, varname, labelname, n) {
  options(digits=10)
  header.str <- paste(c("", unique(temp.df[[labelname]])))#, sep = "\t")
  mean.str <- paste(c("Mean:",laply(unique(temp.df[[labelname]]), function(groupname, temp.df, temp.varname, temp.labelname, n) {
    temp.temp <- temp.df[temp.df[[temp.labelname]] == groupname,]
    return(round(mean( temp.temp[[temp.varname]] ), n))
  }, temp.df, varname, labelname, n)))#, sep = "\t")
  median.str <- paste(c("Median:", laply(unique(temp.df[[labelname]]), function(groupname, temp.df, temp.varname, temp.labelname, n) {
    temp.temp <- temp.df[temp.df[[temp.labelname]] == groupname,]
    round(median( temp.temp[[temp.varname]] ), n)
  }, temp.df, varname, labelname, n)))#, sep = "\t")
  sd.str <- paste(c("SD:",laply(unique(temp.df[[labelname]]), function(groupname, temp.df, temp.varname, temp.labelname, n) {
    temp.temp <- temp.df[temp.df[[temp.labelname]] == groupname,]
    round(sd( temp.temp[[temp.varname]] ), n)
  }, temp.df, varname, labelname, n)))#, sep = "\t")
  iqr.str <- paste(c("IQR:",laply(unique(temp.df[[labelname]]), function(groupname, temp.df, temp.varname, temp.labelname, n) {
    temp.temp <- temp.df[temp.df[[temp.labelname]] == groupname,]
    round(quantile( temp.temp[[temp.varname]], 0.75) - quantile( temp.temp[[temp.varname]], 0.25), n)
  }, temp.df, varname, labelname, n)))#, sep = "\t")
  
  max.l <- max(nchar(c(header.str, mean.str, median.str, sd.str, iqr.str)))
  st.ar <- header.str
  str.list <- laply(list(header.str, mean.str, median.str, sd.str, iqr.str), function(st.ar, max.l) {
    padding <- max.l - nchar(st.ar) + 1
    for (i in 1:length(padding) ) {
      st.ar[i] <-  paste0(c(st.ar[i], rep(" ", times=padding[i])), collapse = "")
    }
    return(paste0(c(st.ar, "\n"), collapse = ""))
  }, max.l)
  to.return <- paste0(str.list, collapse = "")
  return(to.return)
}
# sum.anno(culled.data, "WBLEAN_BASE", "statin.use", 5)                      

n.sum <-  function(dlist, p.test = T) {
	temp <- list()
	tablename <- "Sample Size"

	for ( level in rownames( summary(dlist) ) ) {
		freq <- nrow(dlist[[level]])
		temp[[level]] <- freq
	}
	
	to.return <- temp[[1]]
	for (level in 2:nrow( summary(temp) ) ) {
		to.return <- cbind(to.return, temp[[level]])
	}
	colnames(to.return) <- rep(rownames( summary(dlist)), rep(1, nrow(summary(dlist) ) ) ) 
	rownames(to.return) <- rep(tablename)

	to.return <- cbind(rownames(to.return), to.return)
	to.return <- rbind(c("",rep(rownames( summary(dlist)), rep(1, nrow(summary(dlist) ) ) )), 
					   to.return,
					   rep("", ncol(to.return)))
  if (p.test) {
    to.return <- cbind(to.return, matrix("", ncol=1, nrow=3))
    to.return[1, ncol(to.return)] = "p-value"
  }
	return( to.return )
}

# list.of.varnames <- cont.vars.to.summarize
# list.of.data <- dlist
# data.list <- dlist
# var.info <- list.of.varnames[[1]]
aggregate.cont.sum <- function(list.of.data, list.of.varnames, to.round, p.val = TRUE) {
  
  list.of.tables <- llply(list.of.varnames, function(var.info, data.list, p.val) {
    return(cont.sum(data.list, var.info, to.round, p.val))
  }, list.of.data, p.val)
  
  to.return <- list.of.tables[[1]]
  
  if ( length(list.of.tables)   > 1 ) {
    for (temp.tab in 2:nrow(summary(list.of.tables))) {
      to.return <- rbind(to.return, list.of.tables[[temp.tab]])
    }	
  }
  #Remove duplicate column names
  rows.to.change <- grep(c(names(list.of.data)[1]), to.return[,2])[-1]
  if ( length(rows.to.change) > 0) {
    to.return[rows.to.change,c(2:length(to.return[1,]))] <- ""
  }
  
  return(to.return)						
}

# list.of.data<-culled.groups.full
# list.of.varnames<-cont.vars.to.summarize
# to.round<-tab.prec
# list.of.data <- census_list
# list.of.varnames <- cont_var_list
aggregate.cont.sum.median.iqr <- function(list.of.data, list.of.varnames, to.round) {
  list.of.tables <- llply(list.of.varnames, function(var.info, data.list, to.round) {
    return(cont.sum.median.iqr(data.list, var.info, to.round))
  }, list.of.data, to.round)
  
  to.return <- list.of.tables[[1]]
  
  if ( length(list.of.tables)   > 1 ) {
    for (temp.tab in 2:nrow(summary(list.of.tables))) {
      to.return <- rbind(to.return, list.of.tables[[temp.tab]])
    }  
  }
  #Remove duplicate column names
  rows.to.change <- grep(c(names(list.of.data)[1]), to.return[,2])[-1]
  if ( length(rows.to.change) > 0) {
    to.return[rows.to.change,c(2:length(to.return[1,]))] <- ""
  }
  
  return(to.return)						
}  
  
# varname.array <- cont_var_list[[1]]
  
cont.sum <- function(dlist, varname.array, to.round, p.val = TRUE) {
	temp <- list()
  temp.for.kruskal <- data.frame(samples=NA, levels = NA)
	varname <- varname.array[1]
	tablename <- varname.array[2]
# 	print(varname)
  
	for ( level in rownames( summary(dlist) ) ) {
		x <- as.numeric(dlist[[level]][[varname]])
		Mean = round(mean(x, na.rm = TRUE),to.round)
		SD = round(sd(x, na.rm = TRUE),to.round)
		N = length(x) - sum( is.na(x) ) 
		Missing = round(mean( is.na(x) ) * 100, 2) 
		if (Missing < 0.01) {
			Missing <- "< 0.01"
		}
		temp[[level]] <- rbind(paste(as.character(Mean), " (", as.character(SD), ")", sep=""),  Missing, N) 
    temp.for.kruskal <- rbind(temp.for.kruskal, data.frame( samples = x, levels = rep(level, length(x)) ) )
	}
	
  if (p.val) {
    pval <- kruskal.test(formula = samples ~ factor(levels), data = temp.for.kruskal)$p.value
    if (pval < 0.001) {
      pval <- "< 0.001"
    } else {
      pval <- paste0(" ", as.character(round(pval, digits = 3)))
    }
  }
	to.return <- temp[[1]]
	if (nrow( summary(temp) ) > 1) {
	  for (level in 2:nrow( summary(temp) ) ) {
	    to.return <- cbind(to.return, temp[[level]])
	  }
	}
	to.return <- rbind( c("", rownames( summary(dlist)) ), 
	                    c(tablename, rep("Mean (SD)", ncol(to.return) )),
	                    c("", to.return[1,]),
	                    c("Missing (%)", to.return[2,]),
	                    c("N", to.return[3,]) )
  
	if (p.val) {
	  to.return <- cbind( to.return, matrix("", ncol=1, nrow = nrow(to.return)))
	  to.return[2,ncol(to.return)] <- pval
	}
	return( to.return )
}


# dlist<-culled.groups.full
# varname.array<-cont.vars.to.summarize[[1]]
# to.round<-tab.prec

# to.round<-0
# dlist <- culled.groups.full
# cont.sum.median.iqr(culled.groups.full,cont.vars.to.summarize[[1]], tab.prec)
cont.sum.median.iqr <- function(dlist, varname.array, to.round) {
  temp <- list()
  temp.for.kruskal <- data.frame(samples=NA, levels = NA)
  varname <- varname.array[1]
  tablename <- varname.array[2]
  
  for ( level in rownames( summary(dlist) ) ) {
    x <- dlist[[level]][[varname]]
    Median = round(median(x, na.rm = TRUE),to.round)
    IQR = round(quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE),to.round)
    N = length(x) - sum( is.na(x) ) 
    Missing = round(mean( is.na(x) ) * 100, 2) 
    if (Missing < 0.01) {
      Missing <- "< 0.01"
    }
    temp[[level]] <- rbind(paste(as.character(Median), " (", as.character(IQR), ")", sep=""),  Missing, N) 
    temp.for.kruskal <- rbind(temp.for.kruskal, data.frame( samples = x, levels = rep(level, length(x)) ) )
  }
  
  pval <- kruskal.test(formula = samples ~ factor(levels), data = temp.for.kruskal)$p.value
  if (pval < 0.001) {
    pval <- "< 0.001"
  } else {
    pval <- paste0(" ", as.character(round(pval, digits = 3)))
  }
  
  to.return <- temp[[1]]
  for (level in 2:nrow( summary(temp) ) ) {
    to.return <- cbind(to.return, temp[[level]])
  }
  
  to.return <- rbind( c("", rownames( summary(dlist)) ), 
                      c(tablename, rep("Median (IQR)", ncol(to.return) )),
                      c("", to.return[1,]),
                      c("Missing (%)", to.return[2,]),
                      c("N", to.return[3,]) )
  
  to.return <- cbind( to.return, matrix("", ncol=1, nrow = nrow(to.return)))
  to.return[2,ncol(to.return)] <- pval
  return( to.return )
}
# 
# dlist <- culled.groups
# varname.array <- cat.vars.to.summarize[[8]]
# 
# dlist <- data_list
# varname.array <- table_vars[[14]]

# dlist <- census_list
# varname.array <- cat_var_list[[1]]

cat.sum <- function(dlist, varname.array) {
	temp <- list()
	varname <- varname.array[[1]]
	tablename <- varname.array[[2]]
	row.order <- varname.array[[3]]
	for ( level in rownames( summary(dlist) ) ) {
		#level = "G-G"
		x <- dlist[[level]][[varname]]
		freq <- table(x)
		prop <- round(table(x)/length(x),3)
		table(x)
		leftover.vars <- setdiff(row.order, unique(x))
# 		leftover.vars <- setdiff(union(row.order, unique(x)), intersect(row.order, unique(x)))
		if ( length(leftover.vars) > 0 ) {
			placeholders <- rep(0, length(leftover.vars))
			names(placeholders) <- leftover.vars
			
			freq <- c(freq, placeholders)
			prop <- c(prop, placeholders)
		}
		
    #REMOVE ANY VARIABLES WHICH AREN'T SPECIFIED
    freq <- freq[match(row.order, names(freq))]
		prop <- prop[match(row.order, names(prop))]
		
    if ( length(row.order) > 1 ) {
			new.order <- match(names(freq), row.order)
#         aaply(row.order, 1, function(row.name, ordered.row.name) {
# 										#print(row.name)
# 										return(which(ordered.row.name == row.name))
# 										}, names(freq))
      freq <- freq[new.order]
			prop <- prop[new.order]
		}
		
		temp[[level]] <- matrix(paste0(freq, " (",  prop*100, ")"), ncol = 1 )
    colnames(temp[[level]] ) <- "N (%)"
		rownames(temp[[level]] ) <- names(freq)
	}
	
	to.return <- temp[[1]]
  if (nrow( summary(temp) ) > 1) {
    for (level in 2:nrow( summary(temp) ) ) {
      to.return <- cbind(to.return, temp[[level]])
    }
  }
	colnames(to.return) <- rownames( summary(dlist))
#     rep(rownames( summary(dlist)), rep(2, nrow(summary(dlist) ) ) ) 
	
	#Order
	#Check for "MISSING" row
	if ("MISSING" %in% rownames(to.return) & !("MISSING" %in% row.order) ) {
		missing.row.n <- which(rownames(to.return) == "MISSING")
		Missing <- to.return[missing.row.n,]
		to.return <- to.return[-missing.row.n,]
		to.return <- rbind(to.return, Missing)
	}	
	
	to.return <- cbind(rownames(to.return), to.return)
	to.return <-     rbind(c(tablename, rownames( summary(dlist))), 
	                       to.return,
	                       rep("", ncol(to.return)))
  
  
#     rbind(c(tablename,rep(rownames( summary(dlist)), rep(2, nrow(summary(dlist) ) ) )), 
# 					   to.return,
# 					   rep("", ncol(to.return)))
	return( to.return )
}

# dlist <- culled.groups

#
# varname.array <- cat.vars.to.summarize[[9]]

#cat.sum.perc.only(dlist, varname.array)

cat.sum.perc.only <- function(dlist, varname.array) {
	temp <- list()
	varname <- varname.array[[1]]
  print(varname)
	tablename <- varname.array[[2]]
	row.order <- varname.array[[3]]
	for ( level in rownames( summary(dlist) ) ) {
# 		level = rownames(summary(dlist))[5]
#     names( dlist[[level]])
#     [[varname]]
		x <- as.matrix(dlist[[level]][[varname]])
		#freq <- table(x)
		prop <- round(table(x)/length(x),3)
    

    leftover.vars <- row.order[!row.order %in% as.matrix(unique(x))]
     
    if ( length(leftover.vars) > 0 ) {
			placeholders <- rep(0, length(leftover.vars))
			names(placeholders) <- leftover.vars
			
			#freq <- c(freq, placeholders)
		    prop <- c(prop, placeholders)
		}
		
		if ( length(row.order) > 1 ) {
			new.order <- aaply(row.order, 1, function(row.name, ordered.row.name) {
										#print(row.name)
										#print(ordered.row.name)
										#print(which(ordered.row.name == row.name))
										return(which(ordered.row.name == row.name))
								  		 }, names(prop))
			#freq <- freq[new.order]
			prop <- prop[new.order]
		}
		
		temp[[level]] <- cbind(c("%", prop*100 ) )
	}
	
	to.return <- temp[[1]]
	if (nrow( summary(temp) ) > 1) {
		for (level in 2:nrow( summary(temp) ) ) {
			to.return <- cbind(to.return, temp[[level]])
		}
	}
	
  #colnames(to.return) <- rep(rownames( summary(dlist)), rep(1, nrow(summary(dlist) ) ) ) 
	
	#Order
	#Check for "MISSING" row
	if ("MISSING" %in% rownames(to.return) & !("MISSING" %in% row.order) ) {
		missing.row.n <- which(rownames(to.return) == "MISSING")
		Missing <- to.return[missing.row.n,]
		to.return <- to.return[-missing.row.n,]
		to.return <- rbind(to.return, Missing)
	}	
	
	to.return <- cbind(rownames(to.return), to.return)
	to.return[1,1] <- tablename
# 	to.return <- rbind(c(tablename,rep(rownames( summary(dlist)), rep(1, nrow(summary(dlist) ) ) )), 
# 					   to.return,
# 					   rep("", ncol(to.return)))
	return( to.return )
}

# dlist <- data.list
# variable.name <- varname
#varname.array <- table_vars[[6]]
cat.chisq.test <- function(dlist, varname.array) {
  temp <- list()
  varname <- varname.array[[1]]
  #varname <- list.of.varnames[[3]][[1]]
#   c(dlist[[1]][[varname]],dlist[[3]][[varname]])
#   dim(matrix(dlist[[7]][[varname]]))
#   data.subset <- dlist[[1]]
  #Combine observations
#   print(varname)
  #print(summary(to.table))
  to.table <- ldply(dlist, function(data.subset, variable.name) { 
#     print(class(data.subset[variable.name]))
    return(data.subset[variable.name])
    }, varname)
#   print(varname)
#   print(to.table)
  #print(table(to.table[,1],to.table[,2]))
  pval <- chisq.test(table(to.table[,1],to.table[,2]), simulate.p.value = TRUE)$p.value 
  if (ncol(table(to.table[,1],to.table[,2])) == 1 | min(colSums(table(to.table[,1],to.table[,2]))) == 0 ) {
    pval <- "NA"
  } else {
    if (pval < 0.001) {
      pval <- "< 0.001"
    } else {
      pval <- paste0(" ", as.character(round(pval, digits = 3)))
    }
  }
  return(pval)
}
				  
# dlist <- data.list
# list.of.data <- dlist
# list.of.varnames <- cat.vars.to.summarize
# variable.name <- varname
#varname.array <- cat.vars.to.summarize[[32]]

# dlist <- census_list
# list.of.data <- dlist
# list.of.varnames <- table_vars
# variable.name <- cat_var_list[[1]]
# varname.array <- cat.vars.to.summarize[[32]]

aggregate.cat.sum <- function(list.of.data, list.of.varnames, p.test = TRUE) {
	list.of.tables <- llply(list.of.varnames, function(var.info, data.list) {
                print(var.info[[1]])
								return(cat.sum(data.list, var.info))
								}, list.of.data)
	
	to.return <- list.of.tables[[1]]
  
	if (!p.test) {
	  if ( length(list.of.tables)  > 1 ) {
	    for (temp.tab in 2:nrow(summary(list.of.tables))) {
	      temp.new.table <- list.of.tables[[temp.tab]]
	      to.return <- rbind(to.return, temp.new.table)
	    }	
	  }
	}
  
	#GENERATE P-VALUES
	if (p.test) {
	  list.of.pvalues <- llply(list.of.varnames, function(var.info, data.list) {
	    return(cat.chisq.test(data.list, var.info))
	  }, list.of.data)
	  
	  to.return <- cbind(to.return, matrix(rep("", nrow(to.return), ncol=1 ) ) )
	  to.return[1,ncol(to.return)] <- paste0(" ", list.of.pvalues[[1]])

  
	  if ( length(list.of.tables)  > 1 ) {
	    for (temp.tab in 2:nrow(summary(list.of.tables))) {
	      temp.new.table <- list.of.tables[[temp.tab]]
	      temp.new.table <- cbind(temp.new.table, matrix(rep("", nrow(temp.new.table), ncol=1 ) ) )
	      temp.new.table[1,ncol(to.return)] <- paste0(" ", list.of.pvalues[[temp.tab]])
	      to.return <- rbind(to.return, temp.new.table)
	    }	
	  }
	}	
	
	#Remove duplicate column names
	rows.to.change <- grep(c(names(list.of.data)[1]), to.return[,2])[-1]
	if ( length(rows.to.change) > 0) {
	  to.return[rows.to.change,c(2:(length(to.return[1,])-1))] <- ""
	}
	return(to.return)							
}

#list.of.data <- culled.groups
#list.of.varnames <- cat.vars.to.summarize
# var.info <- list.of.varnames[[3]]
# data.list <- list.of.data
# data.subset <- list.of.data

aggregate.cat.sum.perc.only <- function(list.of.data, list.of.varnames) {
  #GENERATE TABLES
	list.of.tables <- llply(list.of.varnames, function(var.info, data.list) {
								return(cat.sum.perc.only(data.list, var.info))
								}, list.of.data)
  #GENERATE P-VALUES
  list.of.pvalues <- llply(list.of.varnames, function(var.info, data.list) {
    return(cat.chisq.test(data.list, var.info))
  }, list.of.data)
	
  to.return <- list.of.tables[[1]]
  
  to.return <- cbind(to.return, matrix(rep("", nrow(to.return), ncol=1 ) ) )
	to.return[1,ncol(to.return)] <- paste0(" ", list.of.pvalues[[1]])
  
	if ( length(list.of.tables)  > 1 ) {
		for (temp.tab in 2:nrow(summary(list.of.tables))) {
      temp.new.table <- list.of.tables[[temp.tab]]
      temp.new.table <- cbind(temp.new.table, matrix(rep("", nrow(temp.new.table), ncol=1 ) ) )
      temp.new.table[1,ncol(to.return)] <- paste0(" ", list.of.pvalues[[temp.tab]])
			to.return <- rbind(to.return, temp.new.table)
		}	
	}
	
	#Remove duplicate column names
	rows.to.change <- grep(c(names(list.of.data)[1]), to.return[,2])[-1]
	if ( length(rows.to.change) > 0) {
		to.return[rows.to.change,c(2:length(to.return[1,]))] <- ""
	}
	return(to.return)						
}

xtableFormattingMethodsLoaded <- TRUE

