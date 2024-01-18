# Mix type data (normal and Poisson)
gmf_mix = function(Y,
                   X=NULL,
                   family=poisson(),
                   mixp=NULL,
                   p=2,
                   penaltyU=1,
                   penaltyV=0,
                   penaltyBeta=0,
                   intercept=FALSE,
                   maxIter=NULL,
                   gamma=NULL,
                   tol=1e-3,
                   damping=0,
                   init="random",
                   normalize_uv=TRUE,
                   verbose=TRUE,
                   method="airwls",
                   airwls.internal.steps = 1,
                   parallel=1){
  # Adjust parameters to defaults for given methods
  # AIRWLS requires less steps and takes longer steps, but is slower
  # quasi-Newton is fast but since gradients are less accurate it may need shorters steps
  if (is.null(maxIter)){
    if (method == "quasi"){
      maxIter = 1000
    }
    if (method == "airwls"){
      maxIter = 100
    }
  }
  if (is.null(gamma)){
    if (method == "quasi"){
      gamma = 0.01
    }
    if (method == "airwls"){
      gamma = 0.1
    }
  }
  
  # Derive dimensions
  n = dim(Y)[1]
  m = dim(Y)[2]
  
  # Add the intercept as a column of ones in X
  if (intercept)
    X = cbind(matrix(1,n,1),X)
  # Y.norm <- Y[,1:mixp]
  # Y.pois <- Y[,-(1:mixp)]
  
  # X still may be NULL if no intercept and empty input X. Then d = 0
  d = 0
  if (!is.null(X))
    d = dim(X)[2]
  
  # remember where is missing data so that we can keep replacing it
  # with more and more accurate models
  isna = is.na(Y)
  
  if(sum(is.na(Y))>0){
    isna = is.na(Y)
    fun3 <- function(y,fm1){
      q3 <- glm(y~X-1,family = fm1)
      # q3 <- glm(X.train[,1]~V.train,family = poisson())
      if(sum(is.na(y))>0){
        regpart <- X[is.na(y),] %*% q3$coefficients
        if(fm1$family=="binomial"){
          return(sapply(q3[["family"]]$linkinv(regpart),rbinom,n=1,size=1))
        }else if(fm1$family=="poisson"){
          return(sapply(q3[["family"]]$linkinv(regpart),rpois,n=1))
        }else if(fm1$family=="gaussian"){
          return(sapply(q3[["family"]]$linkinv(regpart),rnorm,n=1,sd=1))
        }
      }
      
    }
    Y[isna] = c(unlist(apply(Y[,1:mixp], 2, fun3,fm1=gaussian())),unlist(apply(Y[,-(1:mixp)], 2, fun3,fm1=poisson())))
  }
  # print(Y)
  # Y[isna] = mean(Y,na.rm=TRUE)
  
  # print(Y)
  
  # Initialize U, V and beta using the selected method
  if (init == "svd"){
    initialization = gmf.initial.mix(Y,X,p,mixp,family)
    U = initialization$u #[,1:p,drop=FALSE]
    V = initialization$v #[,1:p,drop=FALSE]
    beta = rbind(initialization$beta)
  }
  if (init=="random"){
    sd = 1e-1/2
    udim = c(n,p)
    vdim = c(m,p)
    betadim = c(d,m)
    U = array(rnorm(prod(udim))/prod(dim(udim))*sd, udim)
    V = array(rnorm(prod(vdim))/prod(dim(vdim))*sd, vdim)
    V.norm = as.matrix(V[1:mixp,])
    V.pois = as.matrix(V[-(1:mixp),])
    beta = array(rnorm(prod(betadim))/prod(betadim)*sd, betadim)
  }
  
  # Remember the last deviance
  llast = Inf
  
  
  # Since regression and latent predictions are used multiple times
  # throughout the computation, we will store them
  regpart = 0
  if (d)
    regpart = X %*% beta
  latentpart = U %*% t(V)
  
  for (i in 1:maxIter){
    # After the first iteration, replace NAs with model values
    if (i > 1){
      Y[isna] = cbind(gaussian()$linkinv(eta.norm),poisson()$linkinv(eta.pois))[isna]
    }
    
    if (method == "airwls"){
      # Perform airwls.internal.steps internal steps of
      # the regularized IRWLS
      for (j in 1:airwls.internal.steps){
        # Get column coefficients
        coefs = slice.glm(cbind(X,U), Y, m, mixp=mixp, coefs = rbind(beta,t(V)),
                          offset = NULL, penalized=penaltyV, parallel = parallel,
                          family = family, method = method, stepsize = gamma)
        if (d)
          beta = coefs[1:d,,drop=FALSE]
        V = t(coefs[(d+1):(d+p),,drop=FALSE])
        # V.norm = as.matrix(V[1:mixp,])
        # V.pois = as.matrix(V[-(1:mixp),])
        
      }
      for (j in 1:airwls.internal.steps){
        # Get row coefficients, correct for regression using the offset
        offset = NULL
        if (d){
          regpart = X %*% beta
          offset = t(regpart)
          offset.norm=offset[1:mixp,]
          offset.pois=offset[-(1:mixp),]
          
        }
        
        funz <- function(ii){
          # ii=1
          eps=1e-3;diff=1e5;iter=1;U.new = U
          while((sqrt(sum(diff^2))> eps*(sqrt(sum(U.new[ii,]^2)))) && iter <= 100) {
            for (l in 1:p) {
              # l=1
              func1 <- function(x,l){
                U.new[ii,l] <- x
                Eta <- regpart + U.new %*% t(V)
                Eta[,-(1:mixp)]=exp(Eta[,-(1:mixp)])
                return(x+sum(V[,l]*(Eta[ii,]-Y[ii,])))
              }
              
              U.new[ii,l] <- uniroot(func1,l=l,c(-10,10))$root
              
              # tt1 <- try(uniroot(func1,l=l,c(-10,10))$root,silent = T)
              # # print(tt1)
              # if(!("try-error" %in% class(tt1))){
              #   U.new[ii,l] <- tt1
              # }
              
            }
            diff=U.new[ii,]-U[ii,]
            U[ii,]=U.new[ii,]
            iter=iter+1
          }
          return(U[ii,])
        }
        
        U <- matrix(sapply(1:n,funz),nrow = n,ncol = p,byrow = T)
        U <- t(t(U)/sqrt(colSums(U^2)))*sqrt(n)
        
        
      }
    }
    else {
      # Directly compute gradients and Hessians
      eta = latentpart + regpart
      eta.norm = eta[,1:mixp]
      eta.pois = eta[,-(1:mixp)]
      
      M.norm = gaussian()$linkinv(eta.norm)
      M.pois = poisson()$linkinv(eta.pois)
      
      M = cbind(M.norm,M.pois)
      
      # Set up helper matrices for computing differentials 
      mueta.norm = matrix(gaussian()$mu.eta(eta.norm),nrow(eta.norm),mixp)
      mueta.pois = poisson()$mu.eta(eta.pois)
      mueta = cbind(mueta.norm,mueta.pois)
      var.norm = matrix(gaussian()$variance(M.norm),nrow(M.norm),mixp)
      var.pois = poisson()$variance(M.pois)
      var.all = cbind(var.norm,var.pois)
      dratio = mueta/var.all
      ddiff = ((Y - M) * dratio)
      ddratio = dratio * mueta
      
      V = update.params(V,U,penaltyV,t(ddiff),t(ddratio),gamma,damping)
      # Update beta if exists
      if (d){
        beta = t(update.params(t(beta),X,penaltyBeta,t(ddiff),t(ddratio),gamma,damping))
        regpart = X %*% beta
      }
      
      funz <- function(ii){
        # ii=1
        eps=1e-3;diff=1e5;iter=1;U.new = U
        while((sqrt(sum(diff^2))> eps*(sqrt(sum(U.new[ii,]^2)))) && iter <= 500) {
          for (l in 1:p) {
            # l=1
            func1 <- function(x,l){
              U.new[ii,l] <- x
              Eta <- regpart + U.new %*% t(V)
              Eta[,-(1:mixp)]=exp(Eta[,-(1:mixp)])
              return(x+sum(V[,l]*(Eta[ii,]-Y[ii,])))
            }
            tt1 <- try(uniroot(func1,l=l,c(-100,100))$root,silent = T)
            # print(tt1)
            if(!("try-error" %in% class(tt1))){
              U.new[ii,l] <- tt1
            }
          }
          diff=U.new[ii,]-U[ii,]
          U[ii,]=U.new[ii,]
          iter=iter+1
        }
        return(U[ii,])
      }
      # newU <- matrix(sapply(1:n,funz),nrow = n,ncol = p,byrow = T)
      U <- matrix(sapply(1:n,funz),nrow = n,ncol = p,byrow = T)
      # V = newV
      # U = newU
      # U <- t(t(U)/sqrt(colSums(U^2)))*sqrt(n)
      
      
      
      
      # U <- t(t(U)/sqrt(colSums(U^2)))*sqrt(n)
      
      
    }
    
    # Update linear predictors
    if (d)
      regpart = X %*% beta
    latentpart = U %*% t(V)
    eta = latentpart + regpart
    eta.norm = eta[,1:mixp]
    eta.pois = eta[,-(1:mixp)]
    
    # Get predicted means
    predicted.norm = gaussian()$linkinv(eta.norm)
    predicted.pois = poisson()$linkinv(eta.pois)
    
    # Compute penalties and mean deviance
    penalties = norm(U,"F")*penaltyU + norm(V,"F")*penaltyV + norm(beta,"F")*penaltyBeta
    # print(sum(regpart==Inf))
    l = (matrix.deviance(predicted.norm, Y[,1:mixp], gaussian()) + matrix.deviance(predicted.pois, Y[,-(1:mixp)], poisson()) + penalties)/prod(dim(Y))
    # l1 = 
    # Diagnostics output: iteration -> deviance    
    if ( verbose){#i %% 10 == 0 &&
      cat("iteration ",i,", PQL = ",l,"\n",sep = "")
    }
    
    # Check the stopping criteria
    if (l==Inf){
      llast = l
    } else{
      if (abs(llast - l) / l < tol){
        break
      }
      else{
        llast = l
      }
    }
  }
  
  # Rotate U and V so that cov(U) is the identity matrix, and V is upper diagonal
  if (normalize_uv){
    decomp = correct.uv(U, V)
    U = decomp$u
    V = decomp$v
    #    assert_that(norm(latentpart - U%*%t(V) ) < 1e-2)
  }
  
  print(paste("Stopped after",i,"iterations with deviance",l))
  
  eta = latentpart + regpart
  eta.norm = eta[,1:mixp]
  eta.pois = eta[,-(1:mixp)]
  
  list(beta = beta,
       u = U,
       v = V,
       fit = cbind(gaussian()$linkinv(eta.norm),poisson()$linkinv(eta.pois)),
       # family = family,
       deviance = l)
}

slice.glm = function(X, Y, slices, mixp=NULL, coefs, offset=NULL, penalized=0, parallel=1, family = poisson(), method="step", stepsize = 1){
  res = c()
  steps = 10
  if (method == "step")
    steps = 1
  
  res = mclapply(1:slices,function(i){
    if (is.null(offset)){
      offset.slice = 0
    } else {
      offset.slice = offset[,i,drop=FALSE]
    }
    ok = !is.na(Y[,i])
    if(is.null(mixp)){
      family1=family
    } else{
      if(i<=mixp){
        family1=gaussian()
      }else{
        family1=poisson()
      }
    }
    
    glm.basic(X[ok,,drop=FALSE],
              Y[ok,i,drop=FALSE],
              beta=coefs[,i],
              family=family1,
              offset = offset.slice,
              stepsize = stepsize,
              penalized = penalized,
              steps=steps)
  },mc.cores = parallel)
  arr = simplify2array(res)
  if (is.vector(arr))
    arr = rbind(arr, NULL)
  arr
}

norm_vec <- function(x) sqrt(sum(x^2))

glm.step = function(X,Y,beta,family,offset = 0, stepsize = 1, penalized = 0)
{
  eta = offset + X %*% beta
  mu = family$linkinv(eta) 
  
  mu.eta = family$mu.eta( eta )
  var.mu = family$variance(mu)
  
  Winv = var.mu / mu.eta**2
  W = mu.eta**2 / var.mu
  
  Z = (eta - offset) + (Y - mu) * Winv
  
  # Use only rows that give reasonable values  
  thresh = 1e20
  keep.rows = !is.nan(c(W)) & !is.infinite(c(W)) & !is.nan(c(Z)) & !is.infinite(c(Z)) & (c(W)>1/thresh) & (c(W)<thresh) & (c(abs(Z))>1/thresh) & (c(abs(Z))<thresh)
  
  if (sum(keep.rows) < ncol(X)+1)
  {
    stop("Too many rows with infinite values")
  }
  Xfull = X[keep.rows,,drop=FALSE]
  Yfull = Z[keep.rows,,drop=FALSE]
  Wfull = c(W)[keep.rows]
  
  if (penalized){
    Yfull = rbind(Yfull,matrix(0, ncol(Xfull), 1))
    Xfull = rbind(Xfull,diag(1, ncol(Xfull)))
    Wfull = c(Wfull,rep(penalized, ncol(Xfull)))
  }
  
  fit = (lsfit(Xfull,Yfull,Wfull,intercept = FALSE))$coefficients#suppressWarnings
  # print(sum(fit==Inf))
  fit*stepsize + (1-stepsize)*beta
}

glm.basic = function(X,Y,beta=NULL,family=gaussian(),tol=1e-5,offset=0, stepsize = 1, penalized = 0, steps=1){
  if (is.null(beta))
    beta = matrix(0,ncol(X),1)
  for (i in 1:steps){
    betaold = beta
    beta = glm.step(X,Y,beta,family,offset=offset,stepsize=stepsize,penalized=penalized)
    if (norm_vec(betaold - beta)/norm_vec(beta) < tol)
      break
  }
  beta
}

qrrange = function(X,q=c(0.05,0.95)){
  range = quantile(c(X),q)
  X[X>range[2]] = range[2]
  X[X<range[1]] = range[1]
  X
}

gmf.initial.mix = function(Y,X,d=min(dim(Y)),mixp,family = poisson()){
  cf = c()
  res = c()
  for (c in 1:mixp){
    yy = Y[,c]
    m = glm(yy ~ X-1, family = gaussian())
    cf = cbind(cf, m$coefficients)
    res = cbind(res, m$residuals)
  }
  for (c in (mixp+1):ncol(Y)){
    yy = Y[,c]
    m = glm(yy ~ X-1, family = poisson())
    cf = cbind(cf, m$coefficients)
    res = cbind(res, m$residuals)
  }
  
  s = svd(res, nu = d, nv = d)
  if(d>1){u = s$u %*% diag(sqrt(s$d[1:d]))
  v = s$v %*% diag(sqrt(s$d[1:d]))}else{
    u = s$u * (sqrt(s$d[1:d]))
    v = s$v * (sqrt(s$d[1:d]))
  }
  
  list(beta = cf,
       u = u,
       v = v
  )
}

matrix.deviance = function(pred, obs, family){
  if (length(pred) == 1){
    pred.matrix = obs
    pred.matrix[] = pred
    pred = pred.matrix
  }
  isna = is.na(obs)
  pred = pred[!isna]
  obs = obs[!isna]
  
  mean(family$dev.resids(obs, pred, 1),na.rm = TRUE)
}




# Rotates U and V such that the covariance of U is diagonal, and V is upper-triangular with a positive diagonal
#' @importFrom whitening whiteningMatrix
correct.uv = function(U, V){
  S = cov(U)
  
  if (ncol(U)==1){
    return(list(u=U/sqrt(c(S)),v=V*sqrt(c(S))))
  }
  
  # Make cov of U identity
  W = whiteningMatrix(S)
  U = U %*% W
  V = V %*% t(solve(W))
  
  # Make V lower triangular
  V.qr = qr(t(V))
  U = U %*% qr.Q(V.qr)
  V = t(qr.R(V.qr))
  
  # Positive diagonal of V
  d = diag(V)
  V = t(sign(d)*t(V))
  U = t(sign(d)*t(U))
  
  list(u=U,v=V)
}

# compute gradients for the parameters. For stability, in estimation remove values
# that are too large
update.params = function(U,V,penalty,ddiff,ddratio,gamma,damping){
  #  stop()
  th = 1e2
  
  # filter out rows with at least one large value of V 
  keep.rows = rowSums(abs(V)>th) < 0.5
  # print(keep.rows)
  nfull = sum(keep.rows)
  
  # correct rows with large ddiff (remove those and reweight rows)
  corrections = nfull / rowSums(abs(ddiff)<th)
  ddiff[abs(ddiff)>th] = 0
  dU = - sweep(ddiff[,keep.rows] %*% V[keep.rows,], MARGIN=1, corrections, `*`) + U*penalty
  
  # correct rows with large ddratio (remove those and reweight rows)
  corrections = nfull / rowSums(abs(ddratio)<th)
  ddratio[abs(ddratio)>th] = 0
  ddU = sweep(ddratio[,keep.rows] %*% V[keep.rows,]**2, MARGIN=1, corrections, `*`) + penalty + damping
  
  gradU = dU / ddU
  U = matrix(U - gamma * gradU,nrow(U),ncol(U)) 
}