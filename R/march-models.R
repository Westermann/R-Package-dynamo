
# Bivariate DCC
bidcc.filter <- function(y,param){
  
  T      <- nrow(y)
  
  if( any(!is.finite(param)) ){ 
    filter = list( rho=rep(0,T) , loglik=-Inf )    
    return( filter ) 
  }
  
  result <- .C( 'bidcc_filter', 
                status = as.integer(0), 
                rho    = as.double(rep(0,T)), 
                eps    = as.double(rep(0,T*2)),
                loglik = as.double(0),
                as.double(param),
                as.double(y),
                as.integer(T), 
                PACKAGE="dynamo")

  filter = list( rho=result$rho , loglik=result$loglik )
  
  return(filter)
}

bidcc.fit <- function(y,opts){
  
  # INPUT 
  if( is.null(opts$param) ){  
    param.init <- c( 0.025 , 0.950 ) 
  }
  else {
    param.init <- opts$param.init 
  }
  if( is.null(opts$fit) ){ 
    fit <- TRUE
  }
  else { 
    fit <- as.logical( opts$fit ) 
  }
  
  # MAIN
  obj  <- function(x){ return( -bidcc.filter(y,x)$loglik ) }  
  
  if( fit==TRUE ){ 
    res <- nlminb( param.init, obj, lower=c(1e-5,0), upper=c(1,1) )
    param.est <- res$par
  }
  else {
    param.est <- param.init 
  }
  
  filter <- bidcc.filter(y,param.est)
  vcv    <- vcv.mle( param.est , obj , 0.0001 * param.est )
  rho    <- data.frame( rho=filter$rho )
  eps    <- data.frame( eps=filter$eps )
  
  param.est           <- as.array(param.est)
  dimnames(param.est) <- list( c('alpha','beta') )
  
  list( param=param.est , 
        fit=rho, 
        resid=eps,
        vcv=vcv ,
        obj=filter$loglik )
}

# BEKK
bekk.filter <- function(y,param=list()){
  
  T <- nrow(y)
  N <- ncol(y)

  result <- .C( 'bekk_filter', 
                S=as.double(rep(0,T*N*N)), 
                C=as.double(rep(0,T*T)), 
                eps=as.double(rep(0,T*N)),
                loglik=as.double(0),
                param=as.double(param),
                y=as.double(as.matrix(y)),
                T=as.integer(T), 
                N=as.integer(N), 
                PACKAGE="dynamo")
  
  filter = list( Sig=result$S, eps=result$eps, loglik=result$loglik, Consts=result$C )
  
  return(filter)
}

bekk.fit <- function(y,opts){

  # Input validation and default setting

  T <- nrow(y)
  N <- ncol(y)
  nparams <- T^2 + 2

  if( is.null(opts$param) ){  
    param.init <- rep(0, nparams) 
  } else {
    param.init <- opts$param 
  }
  if( length(param.init) != nparams ){
    stop('Parameter vector has incorrect length')
  }
  if( is.null(opts$fit) ){
    fit <- TRUE
  } else { 
    fit <- as.logical( opts$fit ) 
  }

  # Optimization
  obj  <- function(x){ return( -bekk.filter(y)$loglik ) }  
  
  if( fit==TRUE ){ 
    res <- nlminb( param.init, obj, lower=c(1e-5), upper=c(1-1e-5) )
    param.est <- res$par
  } else {
    param.est <- param.init 
  }
  
  filter <- bekk.filter(y,param.init)
  
  Sig.C  <- array( filter$S , dim=c(nrow(y),ncol(y),ncol(y)) )
  Sig    <- array( 0 , dim=c(nrow(y),ncol(y),ncol(y)) )
  for( t in 1:nrow(y) ) Sig[t,,] = Sig.C[t,,] %*% t(Sig.C[t,,])
  eps    <- matrix( filter$eps , nrow(y) , ncol(y) );

  # vcv    <- vcv.mle( param.est , obj , 0.0001 * param.est )
  Sig    <- list( Sig=Sig )
  eps    <- data.frame( eps=eps )
  param.est           <- as.array(param.est)
  # dimnames(param.est) <- list( c('a','b',r) )

  list( param=param.est , 
        fit=Sig, 
        C=filter$C,
        resid=eps,
        # vcv=vcv ,
        obj=filter$loglik )
}

# MEWMA
mewma.filter <- function(y,param){
  
  T <- nrow(y)
  N <- ncol(y)
  
  result <- .C( 'mewma_filter', 
                status = as.integer(0), 
                S    = as.double(rep(0,T*N*N)), 
                eps    = as.double(rep(0,T*N)),
                loglik = as.double(0),
                as.double(param),
                as.double(y),
                T=as.integer(T), 
                N=as.integer(N), 
                PACKAGE="dynamo")
  
  filter = list( Sig=result$S , eps=result$eps , loglik=result$loglik )
  
  return(filter)
}

mewma.fit <- function(y,opts){

  # INPUT 
  if( is.null(opts$param) ){  
    param.init <- c( 0.94 ) 
  }
  else {
    param.init <- opts$param.init 
  }
  if( is.null(opts$fit) ){ 
    fit <- TRUE
  }
  else { 
    fit <- as.logical( opts$fit ) 
  }
  
  # MAIN
  obj  <- function(x){ return( -mewma.filter(y,x)$loglik ) }  
  
  if( fit==TRUE ){ 
    res <- nlminb( param.init, obj, lower=c(1e-5), upper=c(1-1e-5) )
    param.est <- res$par
  }
  else {
    param.est <- param.init 
  }
  
  filter <- mewma.filter(y,param.est)
  
  Sig.C  <- array( filter$S , dim=c(nrow(y),ncol(y),ncol(y)) )
  Sig    <- array( 0 , dim=c(nrow(y),ncol(y),ncol(y)) )
  for( t in 1:nrow(y) ) Sig[t,,] = Sig.C[t,,] %*% t(Sig.C[t,,])
  eps    <- matrix( filter$eps , nrow(y) , ncol(y) );
  
  vcv    <- vcv.mle( param.est , obj , 0.0001 * param.est )
  Sig    <- list( Sig=Sig )
  eps    <- data.frame( eps=eps )
  
  param.est           <- as.array(param.est)
  dimnames(param.est) <- list( c('lambda') )

  list( param=param.est , 
        fit=Sig, 
        resid=eps,
        vcv=vcv ,
        obj=filter$loglik )
}
