# for export...

# standardCoefs
# posthocPairwiseT
# cohensD
# etaSquared
# cramersV
# unlibrary
# rmAll
# quantileCut
# sortFrame
# tFrame
# importList
# ciMean
# aad
# modeOf
# maxFreq
# wideToLong
# longToWide
# wideToMV
# longRM
# wideRM
# who

# not for export...
# print.whoList
# pooledSD



standardCoefs <- function( x ) {
  
  # read off the useful info
  term.names <- names(x$coefficients)[-1] # all names except "intercept"
  b <- x$coefficients[-1] # grab coefficients except "intercept"
  
  # construct the design matrix
  predictors <- model.matrix(x$terms, data = x$model)
  predictors <- predictors[ ,-1, drop=FALSE] # (we don't want the intercept)
  predictors <- as.data.frame(predictors)  # hack!!
  
  # standard deviations
  sy <- sd(x$model[[1]])       # sd for the outcome  
  sx <- sapply(predictors,sd)  # sd for the predictors
  
  # now compute beta
  beta <- b * sx / sy

  # convert to matrix
  coefficients <- cbind(b,beta)
  colnames(coefficients) <- c("b","beta")
  rownames(coefficients) <- term.names

  return(coefficients)


}


posthocPairwiseT <- function(x,...) {
  
  # only works for a one-way anova!
  outcome <- x$model[[1]] 
  group <- x$model[[2]]
  
  # run pairwise t.test
  out <- pairwise.t.test(outcome, group, ...)
  
  # alter the names to match the original and return
  var.names <- names(x$model)
  out$data.name <- paste(var.names[1],"and",var.names[2])
  return(out)
}


cohensD <- function(x, y = NULL, data = NULL, method = "pooled",  mu = 0 ) {
  
  # split formula if need be...
  if (is(x,"formula")) {
    if (is.null(data)) { data <- sys.frame(-1) }
    outcome <- eval( x[[2]], data )
    group <- eval( x[[3]], data )
    group <- as.factor(group) 
    if (nlevels(group) != 2L) {
      stop("grouping factor must have exactly 2 levels")
    }
    x <- split(outcome,group)
    y <- x[[2]]
    x <- x[[1]]
  }
      
  if (is.null(y)) { d <- (mean(x) - mu) / sd(x) }  # one sample 
  else {  # two sample... 
    mean.diff <- mean(x) - mean(y)
    sd.est <- switch( EXPR = method,             
          "x.sd" = sd(x),
          "y.sd" = sd(y),
          "pooled" = pooledSD(x,y),
          "corrected" = pooledSD(x,y),
          "raw" = pooledSD(x,y,FALSE),
          "paired" = sd(x-y),
          "unequal" = sqrt( (var(x)+var(y))/2 )
    )               
    d <- mean.diff / sd.est
    if( method == "corrected") { 
      n <- length(x) + length(y)
      d <- d * (n - 3)/(n-2.25)
    }
  }
  return(abs(d))
}

                      
pooledSD <- function(x,y,debias = TRUE) {  
  sq.devs <- (c( x-mean(x), y-mean(y) ))^2
  n <- length(sq.devs)
  if(debias) { psd <- sqrt(sum(sq.devs)/(n-2)) }
  else { psd <- sqrt(sum(sq.devs)/n)}
  return(psd)
}


etaSquared<- function( x ) {
  
  # get the ANOVA table
  anova.table <- summary.aov( x )[[1]]
  
  # read off the useful parts of the ANOVA table
  term.names <- rownames(anova.table)
  ss <- anova.table$`Sum Sq`
  
  # some useful conversion
  k <- length(ss) - 1 # number of terms in the model
  ss.resid <- ss[k + 1] # read off the residual ss
  ss.tot <- sum(ss) # calculate the total ss
  ss <- ss[1:k] # drop "residuals" from the ss vector
  term.names <- term.names[1:k] # drop "residuals" from term names
  
  # compute eta-squared and partial eta-squared
  eta.squared <- ss / ss.tot
  partial.eta.squared <- ss / (ss + ss.resid)

  # convert to matrix
  E <- cbind(eta.squared, partial.eta.squared)
  colnames(E) <- c("eta.sq","partial.eta.sq")
  rownames(E) <- term.names

  # return
  return(E)

}


cramersV <- function (...) {
  
  test <- chisq.test(...)
  chi2 <- test$statistic
  N <- sum(test$observed)
  
  # for GOF test, give same answer as phi
  if (test$method =="Chi-squared test for given probabilities"){
    k <- 2
  } 
  else { # otherwise, correct for size of table
    S <- dim(test$observed)
    k <- min(S)
  }
  
  V <- sqrt( chi2 / (N*(k-1)) )
  names(V) <- NULL
  return(V)
  
}


tFrame <- function(x) {
	if (!is(x,"data.frame")) {
		stop("'tFrame' is intended to apply to data frames only")
	}
	x <- t(x) # coerce to matrix and transpose
	x <- as.data.frame(x)
	return( x )
}

sortFrame <- function(x,..., alphabetical = TRUE){

  dots <- as.list(substitute(list(...)))[-1L] # list of quoted sort terms
  if( length(dots) == 0 ){ return(x) } # do nothing if null arguments
  rel.vars <- unlist(lapply(dots,all.vars)) # which variables are referred to
  y <- lapply(x[rel.vars], xtfrm) # numeric frame that sorts equivalently
  if( alphabetical == TRUE ) { # case conversion if necessary...   
    char.vars <- unlist(lapply(x,is,"character")) # find character vars
    char.vars <- names(which(char.vars)) # relevant variable names
    n <- length(y[[1]]) # number of cases
    for( v in char.vars ) { 
      z <- xtfrm(tolower(x[[v]])) # sorts equivalently to lower case version 
      y[[v]] <- z + y[[v]]/(n+1) # original only breaks ties
    }
  }
  
  ord <- with(y, do.call(order, dots)) # the sort order
  return( x[ord,] ) # sort and return

}


importList <- function(x, ask = TRUE ) {
  
  envir = parent.frame() # import to parent environment
  
  vars <- names(x)  # get variable names
  vars <- make.unique(vars) # make sure the names are unique
  vars <- make.names(vars) # convert to legitimate R names 
  
  if( ask ) {
    cat("Names of variables to be created:\n")
    print(vars)
    ans <- NA
    while( ! (ans %in% c("y","n") ) ) {
      ans <- readline("Create these variables? [y/n] ")
    }
    if (ans == "n") { 
      return( invisible(0) )
    }  
  }
  
  for (v in seq_along(vars)) { # for each variable in x:
    assign(x = vars[v], value = x[[v]], envir = envir)
  }
  return( invisible(1) )
  
} 

quantileCut <- function(x, n, ...) {
  p <- seq(0,1,1/n)
  breaks <- quantile( x, p )
  eps <- (max(x)-min(x)) / 1000
  breaks[1] <- breaks[1] - eps
  breaks[n+1] <- breaks[n+1] + eps  
  out <- cut(x, breaks, ...)
  return( out )
}


unlibrary <- function(package) { 
  env.name <- deparse(substitute(package))   # allow input to drop the quote marks 
  env.name <- paste("package:", env.name, sep="")   # add the "package:" bit (TODO: this is hacky)
  detach(name = env.name, character.only = TRUE)   # now detach it
}


rmAll <- function(ask = TRUE) {
 
  # preliminaries
  env <- parent.frame() # evaluate in the parent frame
  object.list <- objects(env) # and the list of objects
  
  
  # ask the user if they mean it... 
  if ( ask ) {
      
      # don't bother if already empty
      if (length(object.list) == 0) {
        print("Workspace is already empty")
        return( invisible(1) )
      } 
    
      # first, display all the objects...
      cat("Current contents of workspace:\n")
      print( object.list )

      # then ask user for a decision...
      full.prompt <- paste( "Remove all objects? [y/n] ",sep = " ")
      response <- NA
      while( !(response %in% c("y","n")) ) {
        response <- readline( full.prompt )
      } 
  
      # bail out if necessary
      if( response == "n" ) { 
        return( invisible(0) ) 
      }
        
  }

  # remove everything and return
  rm( list = object.list, envir = env ) # ... remove all objects
  return( invisible(1) ) # ... return with invisible flag
  
}


ciMean <- function(x, conf = .95, na.rm = FALSE) {
	
  # remove missing data if requested
  if (na.rm) { x <- x[!is.na(x)] }
  
  # calculate confidence interval using normal approximation
  quantiles <- c( (1-conf)/2 , (1+conf)/2 ) # quantiles of t 
  n <- length(x) # sample size
  CI <- mean(x) + qt(p = quantiles, df = n-1) * sd(x) / sqrt(n) # normal CI
  
  # assign the names attribute
  names(CI) <- paste(100*quantiles,'%',sep="")
  
  # return the confidence interval
  return(CI)
 
}


aad <- function(x, na.rm = FALSE) { 
	if (na.rm) { x <- x[!is.na(x)] }
  y <- mean( abs(x - mean(x)) )
  return(y)
}
  


modeOf <- function(x, na.rm = TRUE) {
  
  na.freq <- 0                                        
  if (na.rm == FALSE) { na.freq <- sum( is.na(x) ) }  # count the NAs if needed
  x <- x[!is.na(x)]                                   # delete NAs  
  obs.val <- unique(x)                                # find unique values
  valFreq <- function(x, y){ sum(y == x) }
  freq <- unlist((lapply( obs.val, valFreq, x )))     # apply for all unique values
  max.freq <- max(freq)                               # modal frequency
  if (na.rm == FALSE & na.freq > max.freq) {
    modal.values <- NA                                # mode is NA if appropriate...
  } else {
    modal.cases <- freq == max.freq                   # otherwise find modal cases
    modal.values <- obs.val[modal.cases]              # and corresponding values
  }
  return( modal.values )
  
}
  
maxFreq <- function(x, na.rm = TRUE) {

  na.freq <- 0                                        
  if (na.rm == FALSE) { na.freq <- sum( is.na(x) ) }  # count the NAs if needed
  x <- x[!is.na(x)]                                   # delete NAs  
  obs.val <- unique(x)                                # find unique values
  valFreq <- function(x, y){ sum(y == x) }
  freq <- unlist((lapply( obs.val, valFreq, x )))     # apply for all unique values
  max.freq <- max(freq, na.freq)                      # modal frequency    
  return( max.freq )                                 
  
}
   
wideToLong <- function(data, rms = NULL, ...) {

  # get the repeated measures info
  if( is.null( rms ) ) {
    rms <- wideRM( data , ...)
  }
  
  # loop over the repeated measures outcomes
  factor.names <- colnames(rms$within)[-(1:2)]
  cells <- unique(rms$within[-(1:2)] )  # all possible combinations
  measures <- levels(rms$within$measure) # all measures taken 
  rownames(cells) <- NULL  # this is for my sanity
  n.treat <- length(cells[[1]]) # how many treatments

  ### this seems inefficient ###
  
  # loop over subjects
  n.subj <- length( data[[1]] )
  for (i in 1:n.subj ) {
    
    X <- data[ i, rms$between, drop = FALSE ]  # get the between info for this subject
    X <- as.data.frame( lapply(X, rep, n.treat) ) # make one copy per treatment
    X[ measures ] <- NA # initialise the measures
    X[ names(cells) ] <- cells  # attach treatment information
    
    # loop over treatments, adding the relevant information  
    for ( j in 1:n.treat ) {
      tmp <- merge( rms$within, cells[j,,drop=FALSE] ) # matching
      m <- as.character( tmp$measure ) # measure name
      wn <- as.character(tmp$wide.name) # wide variable name
      X[ j,m ] <- data[ i,wn ] # assign
    }
    
    # store it
    if (i == 1) { out <- X }  # initialise output
    else { out <- rbind(out,X) } # or append to it
    
  }
  
  return(out)
  
 }


longToWide <- function(data, rms = NULL, ... ) {
  
  # get the repeated measures info
  if( is.null( rms ) ) {
    rms <- longRM( data , ...)
  }
  
  # data frame containing between-subjects variables
  between <- unique( data[ rms$between ] ) 
  
  # counts
  n.subj <- length( between[[1]] )  # number of unique entities
  n.vars <- length( rms$within[[1]] )  # number of within-variables
  n.obs <- length( data[[1]] )  # number of observations in long form
  n.measures <- nlevels( rms$within$measure ) # number of measure vars
  
  # initialise the within-subjects variables
  within <- matrix( nrow = n.subj, ncol = n.vars )  # empty matrix
  within <- as.data.frame( within )  # convert to data frame
  names( within ) <- rms$within$wide.name # name the variables
  
  
  # define the "inject" function, which matches each row of the "full",
  # data frame against a row of the "unique" data frame. I'm not happy with
  # this function at all, and when I get time to go searching I'll replace
  # this part with a better method.
  inject <- function(full, unique) {
    n.obs <- length( full[[1]] )  # rows in the full data frame
    n.unique <- length( unique[[1]] ) # rows in the unique frame
    n.vars <- length( full )  # columns in both
    out <- vector( length = n.obs ) # initialise the output
    for( u in 1:n.unique) {        # loop over all unique cases
      matches <- rep(TRUE, n.obs)  # initially, everything matches
      for( v in 1:n.vars ) {       # loop over all variables
        matches <- matches & ( full[,v]==unique[u,v] ) # does next variable match?
      }
      out[which(matches)] <- u  # update the output
    }
    return(out)  
  }
 
  # now use the inject function to find indices governing the mapping from 
  # the long form rows to the corresponding wide form rows
  ind <- inject( data[rms$between], between )
  
  M <- levels( rms$within$measure )  # measure names
  Tr <- names( rms$within[-(1:2)] )  # treatment names
  
  
  for ( o in 1:n.obs ) {  # loop over rows in the long form matrix

    # logicals indicating the rms$within  cases that match obs o on the treatments
    ind2 <- inject( rms$within[Tr], data[o,Tr,drop=FALSE] ) == 1
    
    for (m in M) {  # loop over the measure variables
      
      tmp <- ind2 & rms$within$measure == m  # matches to ind2 & to the measure name
      wn <- rms$within$wide.name[ tmp ] # grab corresponding wide name
      wn <- as.character( wn )  # ensure this is character     
      within[ ind[o], wn ] <- data[o,m] # assignment
    }
  }
  
  # output
  return( cbind(between, within) )
  
}  



wideToMV <- function( data, rms = NULL, ...) {

  # convert from wide form to multivariate form
  # don't change the wide form names, even though they seem redundant:
  # useful for specifying RM-ANOVA
  
  # get the repeated measures info
  if( is.null( rms ) ) {
    rms <- wideRM( data , ...)
  }
  
  mv <- levels( rms$within$measure )  # the repeated measures 
  out <- list()  # initialise the outputs
  
  for( i in  1:length(mv) ) { # loop over measures
    
    ln <- with( rms$within, wide.name[measure == mv[i]] )  # long names  
    ln <- as.character( ln )  # convert to character 
    out[[i]] <- as.matrix( data[, ln] )  # construct matrix and store it

  } 
  
  names(out) <- mv  # give the matrices the right names
  out[ rms$between ] <- data[ rms$between ]  # add the between subjects info
    
  return(out)
  
}



wideRM <- function(data, treatments = NULL, sep = "_") {
  
    # wideRM() is a function that takes a data frame (data) as input, and constructs 
    # a repeated measures structure as output. The data frame is assumed to be in 
    # wide form, with the variable naming following a logical scheme. Specifically
    # any repeated measure variable should be named: 
    # 
    #      Measurement_Treatment1Level_Treatment2Level 
    #
    # etc. The separator character doesn't have to be "_" though; but must not be
    # present in any purely between subject variables. The "treatments" input is a
    # character vector naming each of the treatments (since the variable names only
    # specify the levels). If unspecified, assumed to be c("treatment.1", 
    # "treatment.2", etc)
    #
    # The output of the function is a list with two elements, "within" and "between".
    # The "between" part is simply a character vector naming the between-subject 
    # variables. The "within" part is a data frame with several variables. The first
    # is the "wide.name" (name of the variable in wide form), the second is "measure"
    # (the name of the measurement taken) and the others contain the various 
    # treatments.
    
  
    # identify the within and between subject variables, and their names
    var.names <- strsplit( names(data), sep ) # a list, one element per var
    n.names <- unlist(lapply( var.names, length)) # a vector indicating how many names
    between.vars <- n.names == 1  # between subj vars have no separators
    within.vars <- !between.vars  # for convenience
  
    # check that the number of separators is consistent:
    if ( length(unique(n.names[within.vars])) != 1) {
      stop("within-subject variables must have the same number of separators")
    }
  
    # store the betweens by name:
    between <- unlist(var.names[ between.vars ])
  
    # construct a data frame with the within data:
    within <- t(as.data.frame( var.names[ within.vars ]))  # matrix
    within <- cbind( names(data)[within.vars], within)  # append widename
    rownames(within) <- NULL
    n.rm <- n.names[within.vars][1] -1 # number of rm factors
    if (is.null(treatments)) { treatments <- paste("treatment",1:(n.rm),sep=".") }
    else { if(length(treatments) != n.rm) { stop("incorrect number of 'treatments' listed") } }
    colnames(within) <- c('wide.name','measure',treatments)
    within <- as.data.frame( within )
  
    # return
    return(list(within = within, between = between))
    
}



longRM <- function( data, treatments, measures, between, sep = "_") {
  
    # longRM() is a function that takes a data frame (data) as input, and constructs 
    # a repeated measures structure as output. The data frame is assumed to be in 
    # long form: the between subject variables (e.g., subject id) are all replicated
    # several times. Each treatment is a separate variable, as is each measurement
    # taken. The user needs to specify which variables correspond to the treatments,
    # which ones correspond to the repeated measures, and which are the between subject
    # variables.
    #
    # The output of the function is a list with two elements, "within" and "between".
    # The "between" part is simply a character vector naming the between-subject 
    # variables. The "within" part is a data frame with several variables. The first
    # is the "wide.name" (name of the variable in wide form), the second is "measure"
    # (the name of the measurement taken) and the others contain the various 
    # treatments.  
  
    # create W, a matrix containing all unique treatment combinations that
    # exist anywhere in the data set.
    W <- unique( data[treatments] )  # grab all observed treatments
    W$wide.name <- NA # initialise wide.name variable
    W$measure <- NA   # initialise measure
    nm <- length( measures )  # number of repeated measures
    W <- W[ c( nm+1,nm+2,1:nm ) ] # permute the order of columns
    
    # populate W with specific info for each measure, and append
    for ( m in 1:nm ) { # loop over measures
      W$measure <- measures[ m ] # record the measure
      if  (m==1) { within <- W } # either initialise "within"...
      else { within <- rbind( within, W ) } # ... or append W to it.
    }
    
    # collapse the measures & treatments to contruct the wide form names
    within$wide.name <- apply( within[-1], 1, paste, collapse = "_" )
    
    # convert to factors
    within$measure <- as.factor( within$measure )
    
    # return both components in a list
    return( list( within = within, between = between ) )
  
}




who <- function(expand = FALSE) {
     
  # extract a list of objects in the parent environment
  envir <- parent.frame()
  varnames <- objects( envir )
  
  # check to see if it's empty
  if ( length(varnames) == 0 ) {
    obj <- character(0)
    class(obj) <- "whoList"
    return(obj)
  }

  # define getWhoInfo as a subfunction
  getWhoInfo <- function(varnames, envir, is.df) {
    if( is.df ) { # data frame?
      df <- varnames
      x <- eval( as.name( df ), envir = envir ) # get the data frame
      varnames <- names(x) # get names of variables
    }
    n <- length(varnames)           # how many objects
    classes <- vector(length = n)   # all output lengths include class info
    info <- vector(length = n)      # output lengths 3+ use info 
    for (i in 1:n) {    
      if( is.df ) {
         var <- varnames[i]
         c <- call("$", as.name(df), as.name(var)) # construct a call 
         v <- eval(c, envir) # evaluate call (i.e., extract variable)      
      } 
      else  { v <- eval( as.name( varnames[i] ), envir ) }
      classes[i] <- class(v)[1]  # class
      if( mode(v) %in% c("logical","numeric","complex","character","list") ) { 
        if( is.null(dim(v)) ) { info[i] <- length(v) } # size = length 
        else { info[i] <- paste( dim(v), collapse = " x ") }  # size = dimensions
      } else { info[i] <- "" } # size = nothing 
    }
    if ( is.df ) { varnames <- paste( "$",varnames,sep="") }  # expand names?
    obj <- data.frame(varnames, classes, info, stringsAsFactors = FALSE)
    names(obj) <- c("Name","Class","Size")
    return(obj)  
  }

  # get info
  obj <-  getWhoInfo(varnames[1], envir, 0) #hack!!!
  for (v in varnames) {
    obj2 <- getWhoInfo(v, envir, 0)
    obj <- rbind(obj, obj2 )
    if (expand) { 
      if( inherits(eval( as.name(v), envir = envir), "data.frame")) {
      	obj2 <- getWhoInfo(v, envir, 1)   # get info 
      	obj <- rbind(obj,obj2)  # add the expanded info to the output
  	  }
  	}
  }
  obj <- obj[-1,] #hack!!!    
  class(obj) <- "whoList"  
  return(obj)
  
}
  
print.whoList <- function(x,...) {

    if( length(x) ==0 ) { 
      cat("No variables found\n") 
    }     
    else{
      obj <- x # copy x
      class(obj) <- "data.frame"  # this is okay
      ind <- grep('\\$', obj$Name)  # variables inside a data frame 
      obj$Name[ind] <- paste(" ", obj$Name[ind], sep = "")   # add indent
      m <- as.matrix(format.data.frame(obj, digits = NULL, na.encode = FALSE)) # matrix
      dimnames(m)[[2]] <- paste( '--', dimnames(m)[[2]], '--', sep = " " ) # tweak header
      row.names(m) <- rep.int("",dim(m)[1]) # hide row names
      print(m, quote = FALSE, right = FALSE, print.gap = 3) # now print   
    }
    return( invisible(x) )
    
}
  
  


