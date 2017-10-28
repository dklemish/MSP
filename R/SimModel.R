simSP <- function(n = 100, dim=1,  
                  model="Brown-Resnick", corrFn = "Matern", 
                  plotFn = TRUE, keepPsi=FALSE,
                  Cstar=15, 
                  sigma2 = 1, nu=2, gam=2, range=0.1, alpha=0.5,
                  seed){
  if(!missing(seed)){
    set.seed(seed)  
  } else{
    seed <- NULL
  }
  
  phi <- 1/range
  
  ## Set up dimensions, grid points, and check nugget
  if(!(dim %in% c(1,2))){
    stop("Only 1 or 2 dimensional simulations supported")
  }else if(dim==1){
    x <- matrix(seq(0, 1, length=n), ncol=1)
  }else if(dim==2){
    x <- as.matrix(expand.grid(seq(0, 1, length=n), seq(0, 1, length=n)))
  }
  
  ## Check model types
  if(model %in% c("Brown-Resnick", "extremeT", "extremeGauss")){
    if(!(corrFn %in% c("Matern", "PowExp", "Exp", "Gaussian"))){
      stop("Correlation function not supported")
    }
    covPar <- list(sigma2, phi, nu, gam)
    
    S <- .Call('_MSP_Sim_SP', PACKAGE='MSP', x, Cstar, model, corrFn, keepPsi, covPar, alpha)
  }else{
    stop("Model not supported")
  }
  
  set.seed(NULL)
  # 
  # oo$x        <- x
  # oo$y        <- y
  # oo$xnames   <- x.names
  # oo$ynames   <- y.names
  # oo$tau.g    <- tau.g
  # oo$dim      <- dimpars
  # 
  # oo$A        <- A
  # oo$R        <- R
  # oo$log.det  <- log.det
  # oo$lp.grid  <- lp.grid
  # 
  # oo$prox     <- prox.grid
  # oo$reg.ix   <- reg.ix
  # oo$hyper    <- hyperPar
  # oo$imcmcpar <- imcmc.par
  # oo$dmcmcpar <- dmcmc.par
  # oo$runtime  <- tm.cpp[3]
  # 
  # samples  <- as.data.frame(oo$parsamp)
  # parnames <- colnames(samples)
  # colnames(samples) <- parnames
  # oo$parsamp <- samples
  # 
  # class(oo) <- "corrQR"
  return(S)
}

simReich <- function(n = 100, dim=1,  
                     plotFn = TRUE, keepPsi=FALSE,
                     alpha=0.5, nknots = 10, bw=1/(nknots+1), rho=0, 
                     seed){
  if(!missing(seed)){
    set.seed(seed)  
  } else{
    seed <- NULL
  }
  
  ## Set up dimensions, grid points, and check nugget
  if(!(dim %in% c(1,2))){
    stop("Only 1 or 2 dimensional simulations supported")
  }else if(dim==1){
    x <- matrix(seq(0, 1, length=n), ncol=1)
  }else if(dim==2){
    x <- as.matrix(expand.grid(seq(0, 1, length=n), seq(0, 1, length=n)))
  }
  
  if(dim==1){
    knots <- matrix(seq(1/(nknots+1), 1-1/(nknots+1), length=nknots), ncol=1)
    Sigma <- matrix(1, nrow=1)
  }else{
    knots <- as.matrix(expand.grid(seq(1/(nknots+1), 1-1/(nknots+1), length=nknots), seq(1/(nknots+1), 1-1/(nknots+1), length=nknots)))
    Sigma <- matrix(c(1, rho, rho, 1), nrow=2)
  }
  S <- .Call('_MSP_Sim_Reich', PACKAGE='MSP', x, knots, Sigma, keepPsi, bw, alpha)
  
  set.seed(NULL)
  
  return(S)
}

simSmith <- function(n = 100, dim=1,  
                     plotFn = TRUE, keepPsi=FALSE,
                     Cstar=10, 
                     sigma2 = 0.1, rho=0, radius = 1,
                     seed){
  if(!missing(seed)){
    set.seed(seed)  
  } else{
    seed <- NULL
  }
  
  ## Set up dimensions, grid points, and check nugget
  if(!(dim %in% c(1,2))){
    stop("Only 1 or 2 dimensional simulations supported")
  }else if(dim==1){
    x <- matrix(seq(0, 1, length=n), ncol=1)
  }else if(dim==2){
    x <- as.matrix(expand.grid(seq(0, 1, length=n), seq(0, 1, length=n)))
  }
  
  if(dim==1){
    Sigma <- matrix(sigma2, nrow=1)
  }else{
    Sigma <- sigma2 * matrix(c(1, rho, rho, 1), nrow=2)
  }
  
  S <- .Call('_MSP_Sim_Smith', PACKAGE='MSP', x, Cstar, keepPsi, Sigma, radius)
  
  set.seed(NULL)
  
  return(S)
}