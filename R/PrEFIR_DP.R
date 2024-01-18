gmf_DP = function(Y,X=NULL,
                  family=poisson(),
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
                  trace=TRUE,
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
    if (method == "airwls"){
      gamma = 0.1
    }
    if (method == "quasi"){
      gamma = 0.01
    }
  }
  
  # Derive dimensions
  n = dim(Y)[1]
  m = dim(Y)[2]
  
  # Add the intercept as a column of ones in X
  if (intercept)
    X = cbind(matrix(1,n,1),X)
  
  # X still may be NULL if no intercept and empty input X. Then d = 0
  d = 0
  if (!is.null(X))
    d = dim(X)[2]
  
  # remember where is missing data so that we can keep replacing it
  # with more and more accurate models
  isna = is.na(Y)
  # Y[isna] = mean(Y,na.rm=TRUE)
  fun3 <- function(y){
    q3 <- glm(y~X-1,family = family)
    # q3 <- glm(X.train[,1]~V.train,family = poisson())
    if(sum(is.na(y))>0){
      regpart <- X[is.na(y),] %*% q3$coefficients
      return(sapply(q3[["family"]]$linkinv(regpart),rpois,n=1))
    }
    
  }
  Y[isna] = unlist(apply(Y, 2, fun3))
  
  # Initialize U, V and beta using the selected method
  if (init == "svd"){
    initialization = gmf.initial_DP(Y,X,p,family)
    U = initialization$u #[,1:p,drop=FALSE]
    V = initialization$v #[,1:p,drop=FALSE]
    beta = rbind(initialization$beta)
  }
  if (init=="random"){
    sd = 1e-1
    udim = c(n,p)
    vdim = c(m,p)
    betadim = c(d,m)
    U = array(rnorm(prod(udim))/prod(dim(udim))*sd, udim)
    V = array(rnorm(prod(vdim))/prod(dim(vdim))*sd, vdim)
    beta = array(rnorm(prod(betadim))/prod(betadim)*sd, betadim)
  }
  dispersion <- rep(0.999, m)
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
      Y[isna] = family$linkinv(latentpart + regpart)[isna]
    }
    
    if (method == "airwls"){
      # Perform airwls.internal.steps internal steps of
      # the regularized IRWLS
      Theta=matrix(dispersion, n, m, byrow = TRUE)
      
      for (j in 1:airwls.internal.steps){
        # Get column coefficients
        coefs = slice.glm_DP(cbind(X,U), Y, m, coefs = rbind(beta,t(V)),dispersion=Theta,
                             offset = NULL, penalized=penaltyV, parallel = parallel,
                             family = family, method = method, stepsize = gamma)
        if (d)
          beta = coefs[1:d,,drop=FALSE]
        V = t(coefs[(d+1):(d+p),,drop=FALSE])
      }
      # Get row coefficients, correct for regression using the offset
      offset = NULL
      if (d){
        # regpart = X %*% beta
        offset = t(regpart)
      }
      for (j in 1:airwls.internal.steps){
        U = t(slice.glm_DP(V, t(Y), n, coefs = t(U),dispersion=t(Theta), offset=offset,
                           penalized=penaltyU, parallel = parallel, family = family,
                           method = method, stepsize = gamma))
      }
      
      # Update linear predictors
      if (d)
        regpart = X %*% beta
      latentpart = U %*% t(V)
      MU=exp(latentpart + regpart)
      # print(MU)
      
      ###optim dispersion (C_MLE)
      dispersion.cur.logfunc <- obj_dispersion(X=Y,MU=MU, dispersion=dispersion)
      # print(dispersion.cur.logfunc)
      
      # q <- try(constrOptim(dispersion, ui=rbind(diag(m),-diag(m)),ci=c(rep(1e-10,m),rep(-1,m)), X=as.matrix(X.train),MU=MU, method = "Nelder-Mead", f = obj_dispersion, control = list(trace = 0,  fnscale = -1, maxit = 500)), silent = TRUE)
      q <- try(constrOptim(dispersion, ui=diag(m),ci=rep(1e-10,m), X=Y,MU=MU, method = "Nelder-Mead", f = obj_dispersion, control = list(trace = 0,  fnscale = -1, maxit = 200)), silent = TRUE)
      # q <- try(constrOptim(dispersion, ui=diag(m),ci=rep(1e-10,m), X=Y,MU=MU, method = "BFGS", f = obj_dispersion, grad = gr_dispersion, control = list(trace = 0,  fnscale = -1, maxit = 500)), silent = TRUE)
      # q <- try(constrOptim(dispersion, ui=rbind(diag(m),-diag(m)),ci=c(rep(1e-10,m),rep(-1,m)), X=Y,MU=MU, method = "BFGS", f = obj_dispersion, grad = gr_dispersion, control = list(trace = 0,  fnscale = -1, maxit = 500)), silent = TRUE)
      # q <- try(constrOptim(new.dispersion, ui=rbind(diag(p),-diag(p)),ci=c(rep(0,p),rep(-1,p)), model.coefs=new.model.coefs, va.mu=new.va.mu, va.sigma = new.va.sigma, method = "BFGS", f = obj_dispersion, grad = gr_dispersion, control = list(trace = 0,  fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ print("error")
      }else{
        if(i > 1 && dispersion.cur.logfunc > q$value){if(trace)
          cat("Optimization of dispersion did not improve on iteration step ",i,"\n")
        }else{
          if(trace) cat("Model parameters dispersion updated","\n")
          dispersion <- q$par;
          if(q$convergence != 0) { if(trace) cat("Optimization of dispersion did not converge on iteration step ", i,"\n") }
        }
      }
      
    }    else {
      # Directly compute gradients and Hessians
      eta = latentpart + regpart
      MU = family$linkinv(eta)
      Theta=matrix(dispersion, n, m, byrow = TRUE)
      
      # Set up helper matrices for computing differentials 
      dratio = family$mu.eta(eta) / family$variance(MU) * Theta
      # print(dratio)
      ddiff = ((Y - MU) * dratio)
      ddratio = dratio * family$mu.eta(eta)
      
      # gaussion()$mu.eta returns a vector instead of a matrix
      # potentially, other families too, so we add the following umbrella conversion
      if (!is.matrix(ddratio)){ 
        ddratio = matrix(ddratio,nrow(eta),ncol(eta)) 
      }
      
      # Update model parameters, updating U and V at once
      newU = update.params(U,V,penaltyU,ddiff,ddratio,gamma,damping)
      V = update.params(V,U,penaltyV,t(ddiff),t(ddratio),gamma,damping)
      U = newU
      
      # Update beta if exists
      if (d){
        beta = t(update.params(t(beta),X,penaltyBeta,t(ddiff),t(ddratio),gamma,damping))
      }
      
      # Update linear predictors
      if (d)
        regpart = X %*% beta
      latentpart = U %*% t(V)
      MU=exp(latentpart + regpart)
      
      ###optim dispersion (C_MLE)
      dispersion.cur.logfunc <- obj_dispersion(X=Y,MU=MU, dispersion=dispersion)
      # print(dispersion.cur.logfunc)
      
      # q <- try(constrOptim(dispersion, ui=rbind(diag(m),-diag(m)),ci=c(rep(1e-10,m),rep(-1,m)), X=as.matrix(X.train),MU=MU, method = "Nelder-Mead", f = obj_dispersion, control = list(trace = 0,  fnscale = -1, maxit = 500)), silent = TRUE)
      q <- try(constrOptim(dispersion, ui=diag(m),ci=rep(1e-10,m), X=Y,MU=MU, method = "Nelder-Mead", f = obj_dispersion, control = list(trace = 0,  fnscale = -1, maxit = 200)), silent = TRUE)
      # q <- try(constrOptim(dispersion, ui=diag(m),ci=rep(1e-10,m), X=Y,MU=MU, method = "BFGS", f = obj_dispersion, grad = gr_dispersion, control = list(trace = 0,  fnscale = -1, maxit = 500)), silent = TRUE)
      # q <- try(constrOptim(dispersion, ui=rbind(diag(m),-diag(m)),ci=c(rep(1e-10,m),rep(-1,m)), X=Y,MU=MU, method = "BFGS", f = obj_dispersion, grad = gr_dispersion, control = list(trace = 0,  fnscale = -1, maxit = 500)), silent = TRUE)
      # q <- try(constrOptim(new.dispersion, ui=rbind(diag(p),-diag(p)),ci=c(rep(0,p),rep(-1,p)), model.coefs=new.model.coefs, va.mu=new.va.mu, va.sigma = new.va.sigma, method = "BFGS", f = obj_dispersion, grad = gr_dispersion, control = list(trace = 0,  fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ print("error")
      }else{
        if(i > 1 && dispersion.cur.logfunc > q$value){if(trace)
          cat("Optimization of dispersion did not improve on iteration step ",i,"\n")
        }else{
          if(trace) cat("Model parameters dispersion updated","\n")
          dispersion <- q$par;
          if(q$convergence != 0) { if(trace) cat("Optimization of dispersion did not converge on iteration step ", i,"\n") }
        }
      }
      
    }
    
    
    # Get logfunc
    logfunc = obj_dispersion(X=Y,MU=MU, dispersion=dispersion)
    
    # Compute penalties and mean deviance
    penalties = norm(U,"F")*penaltyU + norm(V,"F")*penaltyV + norm(beta,"F")*penaltyBeta
    predicted = family$linkinv(latentpart + regpart)
    # l = (matrix.deviance(predicted, Y, family) + penalties)/prod(dim(Y))
    l = (-logfunc + penalties)/prod(dim(Y))
    
    # Diagnostics output: iteration -> deviance    
    if (i %% 10 == 0 && verbose){
      cat("iteration ",i,", PQL = ",l,"\n",sep = "")
    }
    
    # Check the stopping criteria
    if (abs(llast - l) / abs(l) < tol){
      break
    }
    else{
      llast = l
    }
  }
  
  # Rotate U and V so that cov(U) is the identity matrix, and V is upper diagonal
  # if (normalize_uv){
  #   decomp = correct.uv(U, V)
  #   U = decomp$u
  #   V = decomp$v
  #   #    assert_that(norm(latentpart - U%*%t(V) ) < 1e-2)
  # }
  
  print(paste("Stopped after",i,"iterations with deviance",l))
  
  list(beta = beta,
       u = U,
       v = V,
       dispersion = dispersion,
       fit = family$linkinv(latentpart + regpart),
       family = family,
       deviance = l)
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


update.params = function(U,V,penalty,ddiff,ddratio,gamma,damping){
  #  stop()
  th = 1e2
  
  # filter out rows with at least one large value of V 
  keep.rows = rowSums(abs(V)>th) < 0.5
  
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

norm_vec <- function(x) sqrt(sum(x^2))

glm.step_DP = function(X,Y,beta,dispersion=1,family,offset = 0, stepsize = 1, penalized = 0)
{
  eta = offset + X %*% beta
  mu = family$linkinv(eta) 
  
  mu.eta = family$mu.eta( eta )
  var.mu = family$variance(mu)
  
  Winv = var.mu / mu.eta**2
  W = mu.eta**2 / var.mu * dispersion
  
  Z = (eta - offset) + (Y - mu) * Winv
  
  # Use only rows that give reasonable values  
  thresh = 1e20
  keep.rows = !is.nan(c(W)) & !is.infinite(c(W)) & !is.nan(c(Z)) & !is.infinite(c(Z)) & (c(W)>1/thresh) & (c(W)<thresh) & (c(abs(Z))>1/thresh) & (c(abs(Z))<thresh)
  # print(keep.rows)
  
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
  
  fit = suppressWarnings(lsfit(Xfull,Yfull,Wfull,intercept = FALSE))$coefficients
  fit*stepsize + (1-stepsize)*beta
}

glm.basic_DP = function(X,Y,beta=NULL,dispersion,family=gaussian(),tol=1e-5,offset=0, stepsize = 1, penalized = 0, steps=1){
  if (is.null(beta))
    beta = matrix(0,ncol(X),1)
  for (i in 1:steps){
    betaold = beta
    beta = glm.step_DP(X,Y,beta,dispersion=dispersion,family,offset=offset,stepsize=stepsize,penalized=penalized)
    if (norm_vec(betaold - beta)/norm_vec(beta) < tol)
      break
  }
  beta
}

slice.glm_DP = function(X, Y, slices, coefs,dispersion, offset=NULL, penalized=0, parallel=1, family = poisson(), method="step", stepsize = 1){
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
    glm.basic_DP(X[ok,,drop=FALSE],
                 Y[ok,i,drop=FALSE],
                 beta=coefs[,i],
                 dispersion=dispersion[,i],
                 family=family,
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

gmf.initial_DP = function(Y,X,d=min(dim(Y)),family = poisson()){
  cf = c()
  res = c()
  for (c in 1:ncol(Y)){
    yy = Y[,c]
    if (family$family=="binomial"){
      yy[is.na(yy)] = rbinom(sum(is.na(yy)), 1, mean(yy,na.rm=TRUE))
      yy = as.factor(yy)
    }
    
    m = glm(yy ~ X-1, family = family)
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

obj_dispersion <- function(X,MU, dispersion=NULL){
  n = dim(X)[1]
  m = dim(X)[2]
  theta = dispersion
  Theta = matrix(theta, n, m, byrow = TRUE)
  X_nz <- X
  X_nz[which(X==0)] <- 1
  
  y <-  0.5*sum(log(Theta)) + sum(Theta*X*log(MU)) + sum(Theta*X) - sum(Theta*X*log(X_nz)) - sum(Theta*MU)
  
  return(y)
}

gr_dispersion <- function(X,MU, dispersion=NULL){
  n = dim(X)[1]
  m = dim(X)[2]
  theta = dispersion
  Theta = matrix(theta, n, m, byrow = TRUE)
  X_nz <- X
  X_nz[which(X==0)] <- 1
  y <-  0.5/dispersion + colSums(X*log(MU)) + colSums(X) - colSums(X*log(X_nz)) - colSums(MU)
  
  return(c(y))
}



