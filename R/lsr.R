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
# who
# permuteLevels
# colCopy
# rowCopy
# expandFactors
# xfun

# not for export...
# print.whoList
# pooledSD


wideToLong <- function( data, within="within", sep="_", split=TRUE) {
	
	ind <- grep(sep,names(data)) # indices of variables that are repeated
	idvar <- names(data)[-ind] # names of id varibles
	v.names <- unique(gsub( paste(sep,".*$",sep=""), "", names(data)[ind])) # measure var names
	times <- unique(gsub( paste("^(",paste(v.names,collapse="|"),")",sep,sep=""), "", names(data)[ind])) # measure 'time' names
	varying <- list()
	for( i in seq_along(v.names) ) varying[[i]] <- grep( paste("^",v.names[i],sep,sep=""), names(data)[ind], value=TRUE)
	
	tmp <-make.unique(c(names(data),"withintmp"))
	within.tmp <- tmp[length(tmp)]
	x<-reshape( data, idvar=idvar, varying=varying, direction="long", times=times, v.names=v.names, timevar=within.tmp )
	
	if( split==TRUE & length( grep(sep,times))>0 ) { # split multiple treatments into two factors?
		split.treatments <- tFrame(as.data.frame(strsplit(x[,within.tmp],sep)))
		if( length(within)==1) { 
			names(split.treatments) <- paste(within,1:length(split.treatments),sep="") 
		} else {
			if( length(within) == length(split.treatments)) {
				names(split.treatments) <- within 
			} else { stop( "length of 'within' is incorrect" )}
		}
		x <- x[,setdiff(names(x),within.tmp)] # delete collapsed treatment
		x <- cbind(x,split.treatments) # append split treatment
	} else { 
		x[,within.tmp]<- factor(x[,within.tmp])
		names(x)[grep(within.tmp,names(x))] <- within 
	}
	rownames(x) <- NULL
	names(x) <- make.unique(names(x))
	return(x)
	
}

longToWide <- function( data, formula, sep="_") {
	
	within <- all.vars(formula[-2])
	v.names <- all.vars(formula[-3])
	idvar <- setdiff(names(data),c(within,v.names)) 
	
	if( length(within)>1 ) { 
		collapsed.treatments <- apply(as.matrix(data[,within]),1,paste,collapse="_")
		data <- data[,setdiff(names(data),within)] # delete split treatments
		data$within <- collapsed.treatments # append collapsed treatment
		within <- "within"
	}
	times <- unique( data[,within]) # measure 'time' names
	varying <- list()
	for( i in seq_along(v.names) ) varying[[i]] <- paste(v.names[i],times, sep=sep)
	
	x<-reshape( data, idvar=idvar, varying=varying, direction="wide", times=times, v.names=v.names, timevar=within)
	rownames(x) <- NULL
	return(x)
	
}

expandFactors <- function( data, ... ) {
 
  df <- model.matrix( as.formula( paste("~",names(data),collapse="+")), data, ... )
  df <- df[,-1]
  attr(df,"contrasts") <- NULL
  attr(df,"assign") <- NULL

  return(df)
}

xfun <- function( formula, data = NULL, fun, subset = NULL, split = TRUE, ... ) {
  
  # check that formula is a two-sided formula
  if( class(formula) != "formula" ) stop( "formula input must be a formula") 
  if( length(formula) < 3) stop( "formula must be two-sided" )
  
  # get the list of variables included in the formula
  vars <- all.vars( formula )  
  
  # check that the formula has only one LHS variable
  if( length(vars) != length(all.vars(formula[-2])) +1 ) stop( "only one outcome variable may be specified" )
  
  # if variables are not in data frame, check that they are appropriate and then construct one
  if( is.null(data) ) {
    
    # check that the variables exist
    okay <- sapply(vars,exists)
    if ( !all(okay)) stop( paste("unable to find variable(s):",paste(names(okay)[!okay],collapse=", ")) )
    
    # place all the variables in a list
    data<-lapply(vars,get) 
    
    # check that the variables are all vectors
    if(!all(sapply(data,is.vector))) stop( "variables must be vectors" )
    
    # check that the vectors are of the same length
    var.lengths<-sapply(data,length)
    if(all(var.lengths!=var.lengths[1])) stop( "variables are not of the same length" )
    
    # construct data frame
    data<-as.data.frame(data)
    names(data)<-vars
  } else {
    
    # check that the variables exist within the data frame
    
  }
  
  # use tapply to apply fun separately to outcome by group
  if( is.null(subset) ) { 
    x <- tapply( data[,vars[1]], data[,vars[-1]], fun, ... )
  } else {
    x <- tapply( data[subset,vars[1]], data[subset,vars[-1]], fun, ... ) 
  }
  
  # if the output is mode list & the user has requested a split output, split it
  if( mode(x) == "list" & split == TRUE )  {
    
    # check to see if the contents of x are always vectors
    is.vector <- apply(x,1:length(dim(x)),function(l){is.vector(l[[1]])})
    
    if( all(is.vector) ){ 
      max.len <- max(apply(x,1:length(dim(x)),function(l){length(l[[1]])}))
      y <- list()
      for( i in 1:max.len) y[[i]] <- apply(x,1:length(dim(x)),function(l,ind){l[[1]][ind]}, i )
      x <- y
    } else {
      warning( paste("output of function is not a vector: returning array of lists, not list of arrays" ))
    }
    
  }
  
  return(x)
  
}

colCopy <- function(x,times, dimnames=NULL ) {
  if( is.null(dimnames) ) dimnames<-list(names(x),character(0)) 
  matrix( x, length(x), times, byrow=FALSE, dimnames )
}

rowCopy <- function(x,times, dimnames=NULL ) {
  if( is.null(dimnames) ) dimnames<-list(character(0),names(x))
  matrix( x, times, length(x), byrow=TRUE, dimnames )
}

permuteLevels <- function(x,perm,ordered = is.ordered(x),invert=FALSE) {
  if(invert){ perm <- order(perm) }
  y <- factor( x = order(perm)[as.numeric(x)],
                 levels = seq_along(levels(x)),
                 labels = levels(x)[perm],
                 ordered = ordered ) 
  return(y)
}


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


etaSquared<- function( x, type = 2, anova = FALSE ) {
  
  if( type == 1) {
      
    ss <- anova(x)[,"Sum Sq",drop=FALSE]  # Type 1 SS 
    ss.res <- ss[dim(ss)[1],]  # Full model RSS
    ss.tot <- sum( ss )  # Total SS 
    ss <- ss[-dim(ss)[1],,drop=FALSE]
    
  } else { if (type == 2) {
      
    ss.tot <- sum(( x$model[,1] - mean(x$model[,1]) )^2)
    ss.res <- sum(( x$residuals)^2)
    t <- attr(x$terms,"term.labels")
    ord <- attr(x$terms,"order")
    max.ord <- max(ord)
    ss <- matrix(NA,length(t),1)
    rownames(ss) <- t
    for( i in seq_along(ss) ) {
      
      # find the i-th term plus its higher order terms
      if( length( grep(":",t[i]) )>0 ) {
        ti <- seq_along(t) 
        vi <- unlist(strsplit(t[i],":"))
        for( v in vi ) ti <- intersect(ti, grep(v,t)) 
      } else {
        ti <- grep(t[i],t)
      }
    
      # now calculate the corresponding ss value
      m0 <- lm( x$terms[-ti], x$model )  # remove all of these
      if( ord[i] < max.ord) {
        m1 <- lm( x$terms[-setdiff(ti,i)], x$model ) # remove all except i-th term
        ss[i] <- anova(m0,m1)$`Sum of Sq`[2] # get the ss value
      } else {
        ss[i] <- anova(m0,x)$`Sum of Sq`[2]
      }
    }
    
    
  } else { if (type == 3) { 
    
    ss <- drop1(x,scope=x$terms)[-1,"Sum of Sq",drop=FALSE] # Type 3 SS
    ss.res <- ss[dim(ss)[1],]  # Full model RSS
    ss.tot <- sum(( x$model[,1] - mean(x$model[,1]) )^2)
    
  } else { 
    stop("type must be equal to 1,2 or 3") 
  }}} 
  
  
  if( anova == FALSE) {
    eta2 <- ss / ss.tot
    eta2p <- ss / (ss + ss.res)
    E <- cbind(eta2, eta2p)
    rownames(E) <- rownames(ss)
    colnames(E) <- c("eta.sq","eta.sq.part")
  
  } else {
    ss <- rbind( ss, ss.res )
    eta2 <- ss / ss.tot
    eta2p <- ss / (ss + ss.res)  
    k <- length(ss)
    eta2p[k] <- NA
    df <- anova(x)[,"Df"] # lazy!!!
    ms <- ss/df
    Fval <- ms / ms[k]
    p <- 1-pf( Fval, df, rep.int(df[k],k) ) 
    E <- cbind(eta2, eta2p, ss, df, ms, Fval, p)
    E[k,6:7] <- NA
    colnames(E) <- c("eta.sq","eta.sq.part","SS","df","MS","F","p")
    rownames(E) <- rownames(ss) 
    rownames(E)[k] <- "Residuals"
  }
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
  if(class(x)=="factor") modal.values <- as.character(modal.values)
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
  
  


