

# Returns the number of unique values in an array
cUnique <- function(x, count.na = FALSE) {
  members <- unique(x)
  if (NA %in% members & !count.na) {
    return(length(members) - 1)
  } else {
    return(length(members))
  }
}


# Loads and combines files in a folder that match a specified
# string
# if merge == TRUE, a dataframe of merged files is returned, 
#   else a list with an element corresponding to each file is returned
openFilesInDirectory <- function(directory, match_string, merge = FALSE, delim_str =",") {
  
  file_array <-  paste0(directory, "/", list.files(directory)[grep(pattern=match_string, list.files(directory))])
  
  data_list <- llply(file_array, function(file_path, delim_str) {
    print(file_path)
    to_return <- read.table(file = file_path, header = TRUE, sep = delim_str, stringsAsFactors = FALSE, fill=TRUE, quote="\"" )
    to_return["loaded_file_name"] <- tail(strsplit(file_path, "/")[[1]],1)
    return(to_return)
  }, delim_str)

  if(merge |  length(file_array) == 1) {
    data_list <- ldply(data_list, identity)
  }
  
  return(data_list)
}

# Converts factors into numeric or character arrays
factorConvert <- function(var, to_type = "numeric") {
  
  if (to_type == "numeric" ){
    to_return <- as.numeric(levels(var)[as.numeric(var)])
  }
  
  if (to_type == "character" ){
    to_return <- levels(var)[as.numeric(var)]
  }
  return(to_return)
}

# Shows the extent of NAs in a data frame

showNAs <- function(temp_df) {
  df_names <- names(temp_df)
  to_return <- ldply(df_names, function(var_name, temp_df) {
    to_return <- temp_df[[var_name]]
    to_return[to_return == "MISSING"] <- NA
    return( data.frame(var = var_name, 
                       NA_count = sum(is.na(to_return)),
                       NA_mean = mean(is.na(to_return))))
  }, temp_df)
  return(to_return)
}

# Shows the unique values for each factor or character
# variable in a data frame temp_df
# nchars specifies how many characters to display for each 
#  unique value
showUniqueValues <- function(temp_df, nchars = 5) {
  df_names <- names(temp_df)
  if (nchars == 0) {nchars = 10e6}
  to_return <- ldply(df_names, function(var_name, temp_df) {
    to_return <- temp_df[[var_name]]
    if (class(to_return) == "character" | class(to_return) == "factor") {
      return( data.frame(var = var_name, 
                         values = paste0(substr(unique(to_return), 1, nchars), collapse = ", "),
                         kurtosis = NA))
    }
    if (!(class(to_return) == "character" | class(to_return) == "factor")) {
      return( data.frame(var = var_name, 
                         values = "NA",
                         kurtosis = kurtosis(to_return, na.rm = TRUE)))
    }
  }, temp_df)
  return(to_return)
}

# Shows the first word of each unique value for character or
# factor in a temp_df 
showUniqueValuesFirstWord <- function(temp_df) {
  df_names <- names(temp_df)
  to_return <- ldply(df_names, function(var_name, temp_df) {
    to_return <- temp_df[[var_name]]
    
    if (class(to_return) == "character" | class(to_return) == "factor") {
      return( data.frame(var = var_name, 
                         values = paste0(substr(unique(to_return), 1, 
                                                laply(strsplit(unique(to_return), " "), function(entry) {
                                                  return(nchar(entry[1]))
                                                }) ), 
                                         collapse = ", "),
                         kurtosis = NA))
    }
    
    if (!(class(to_return) == "character" | class(to_return) == "factor")) {
      return( data.frame(var = var_name, 
                         values = "NA",
                         kurtosis = kurtosis(to_return, na.rm = TRUE)))
    }
  }, temp_df)
  return(to_return)
}


##### KK 2013-04-08: Added manual list of k values at which to calculate gap widths
##### From here: https://svn.r-project.org/R-packages/trunk/cluster/R/clusGap.R

#### Originally from orphaned package SLmisc
#### (Version: 1.4.1, 2007-04-12, Maintainer: Matthias Kohl <kohl@sirs-lab.com>)
#### License: GPL (version 2 or later)
####
#### which said
####  "function corresponds to function gap in package SAGx"

## MM: SAGx is now in Bioconductor --- 1.10.1{devel} or 1.11.1{release}
##     had gap() *corrected* to re-cluster using FUNcluster --> see ./gap-SAGx.R.~orig~
##
## MM: Package 'lga' -- has gap() and lga and robust lga [-> UBC]
##    - it uses  boot() nicely  [2012-01: ORPHANED because  Justin Harrington is amiss]
## MM: renamed arguments, and changed almost everything

# x <- sampled_data[, cause_index]
# FUNcluster <- pam1
# k_arr <- seq(2,10,4)
# B <- 10
clusGapManualK <- function (x, FUNcluster, K.max, k_arr = seq(2,10,1), B = 100, verbose = interactive(), ...)
{
  k_length <- length(k_arr)
  K.max <- max(k_arr)
  stopifnot(is.function(FUNcluster), length(dim(x)) == 2, K.max >= 2,
            (n <- nrow(x)) >= 1, (p <- ncol(x)) >= 1)
  if(B != (B. <- as.integer(B)) || (B <- B.) <= 0)
    stop("'B' has to be a positive integer")
  
  if(is.data.frame(x))
    x <- as.matrix(x)
  ii <- seq_len(n)
  W.k <- function(X, kk) {
#     X <- x
#     kk <- k_arr[k]
    clus <- if(kk > 1) FUNcluster(X, kk)$cluster else rep.int(1L, nrow(X))
    ##                 ---------- =  =       -------- kmeans() has 'cluster'; pam() 'clustering'
    0.5* sum(vapply(split(ii, clus),
                    function(I) { xs <- X[I,, drop=FALSE]
                                  sum(dist(xs)/nrow(xs)) }, 0.))
  }
  logW <- E.logW <- SE.sim <- numeric(k_length)
  if(verbose) cat("Clustering k = ", k_arr[1], ", ", k_arr[2], ",..., K.max (= ",K.max,"): .. ", sep='')
  for(k in 1:k_length) 
    logW[k] <- log(W.k(x, k_arr[k]))
  if(verbose) cat("done\n")
  
  ## Scale 'x' into "hypercube" -- we later fill with H0-generated data
  xs <- scale(x, center=TRUE, scale=FALSE)
  m.x <- rep(attr(xs,"scaled:center"), each = n)# for back transforming
  V.sx <- svd(xs)$v
  rng.x1 <- apply(xs %*% V.sx, # = transformed(x)
                  2, range)
  
  logWks <- matrix(0., B, k_length) 
  if(verbose) cat("Bootstrapping, b = 1,2,..., B (= ", B,
                  ")  [one \".\" per sample]:\n", sep="")
  for (b in 1:B) {
#     b <- 2
    ## Generate "H0"-data as "parametric bootstrap sample" :
    z1 <- apply(rng.x1, 2,
                function(M, nn) runif(nn, min=M[1], max=M[2]),
                nn=n)
    z <- tcrossprod(z1, V.sx) + m.x # back transformed
    for(k in 1:k_length) {
      logWks[b,k] <- log(W.k(z,k_arr[k]))
    }
    if(verbose) cat(".", if(b %% 50 == 0) paste(b,"\n"))
  }
  if(verbose && (B %% 50 != 0)) cat("",B,"\n")
  E.logW <- colMeans(logWks)
  SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
  structure(class = "clusGap",
            list(Tab = cbind(k = matrix(k_arr, ncol = 1), logW, E.logW, gap = E.logW - logW, SE.sim),
                 ## K.max == nrow(T)
                 n = n, B = B, FUNcluster=FUNcluster))
}



g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

blankground <- function() {
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        panel.margin = unit(0,"null"),
        plot.margin = rep(unit(0,"null"),4),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank()
        
  )
}

miscRLoaded <- TRUE




