PrEFIR = function (Y, X = NULL, family = poisson(), p = 2, penaltyU = 1, 
                  penaltyV = 0, penaltyBeta = 0, intercept = FALSE, maxIter = NULL, 
                  gamma = NULL, tol = 0.001, damping = 0, init = "random", 
                  normalize_uv = TRUE, verbose = TRUE, method = "airwls", 
                  airwls.internal.steps = 1, parallel = 1) 
{
  if (is.null(maxIter)) {
    if (method == "quasi") {
      maxIter = 1000
    }
    if (method == "airwls") {
      maxIter = 100
    }
  }
  if (is.null(gamma)) {
    if (method == "quasi") {
      gamma = 0.01
    }
    if (method == "airwls") {
      gamma = 0.1
    }
  }
  n = dim(Y)[1]
  m = dim(Y)[2]
  if (intercept) 
    X = cbind(matrix(1, n, 1), X)
  d = 0
  if (!is.null(X)) 
    d = dim(X)[2]
  if (init == "svd") {
    initialization = gmf.initial(Y, X, p, family)
    U = initialization$u
    V = initialization$v
    beta = rbind(initialization$beta)
  }
  if (init == "random") {
    sd = 0.1/2
    udim = c(n, p)
    vdim = c(m, p)
    betadim = c(d, m)
    U = array(rnorm(prod(udim))/prod(dim(udim)) * sd, udim)
    V = array(rnorm(prod(vdim))/prod(dim(vdim)) * sd, vdim)
    beta = array(rnorm(prod(betadim))/prod(betadim) * sd, 
                 betadim)
  }
  llast = Inf
  isna = is.na(Y)
  Y[isna] = mean(Y, na.rm = TRUE)
  regpart = 0
  if (d) 
    regpart = X %*% beta
  latentpart = U %*% t(V)
  for (i in 1:maxIter) {
    if (i > 1) {
      Y[isna] = family$linkinv(latentpart + regpart)[isna]
    }
    if (method == "airwls") {
      for (j in 1:airwls.internal.steps) {
        coefs = slice.glm(cbind(X, U), Y, m, coefs = rbind(beta, 
                                                           t(V)), offset = NULL, penalized = penaltyV, 
                          parallel = parallel, family = family, method = method, 
                          stepsize = gamma)
        if (d) 
          beta = coefs[1:d, , drop = FALSE]
        V = t(coefs[(d + 1):(d + p), , drop = FALSE])
      }
      for (j in 1:airwls.internal.steps) {
        offset = NULL
        if (d) {
          regpart = X %*% beta
          offset = t(regpart)
        }
        U = t(slice.glm(V, t(Y), n, coefs = t(U), offset = offset, 
                        penalized = penaltyU, parallel = parallel, 
                        family = family, method = method, stepsize = gamma))
        U <- t(t(U)/sqrt(colSums(U^2))) * sqrt(n)
      }
    }
    else {
      eta = latentpart + regpart
      M = family$linkinv(eta)
      dratio = family$mu.eta(eta)/family$variance(M)
      ddiff = ((Y - M) * dratio)
      ddratio = dratio * family$mu.eta(eta)
      if (!is.matrix(ddratio)) {
        ddratio = matrix(ddratio, nrow(eta), ncol(eta))
      }
      newU = update.params(U, V, penaltyU, ddiff, ddratio, 
                           gamma, damping)
      V = update.params(V, U, penaltyV, t(ddiff), t(ddratio), 
                        gamma, damping)
      U = newU
      if (d) {
        beta = t(update.params(t(beta), X, penaltyBeta, 
                               t(ddiff), t(ddratio), gamma, damping))
      }
    }
    if (d) 
      regpart = X %*% beta
    latentpart = U %*% t(V)
    predicted = family$linkinv(latentpart + regpart)
    penalties = norm(U, "F") * penaltyU + norm(V, "F") * 
      penaltyV + norm(beta, "F") * penaltyBeta
    l = (matrix.deviance(predicted, Y, family) + penalties)/prod(dim(Y))
    if (verbose) {
      cat("iteration ", i, ", PQL = ", l, "\n", sep = "")
    }
    if (l == Inf) {
      llast = l
    }
    else {
      if (abs(llast - l)/l < tol) {
        break
      }
      else {
        llast = l
      }
    }
  }
  if (normalize_uv) {
    decomp = correct.uv(U, V)
    U = decomp$u
    V = decomp$v
  }
  print(paste("Stopped after", i, "iterations with deviance", 
              l))
  list(beta = beta, u = U, v = V, fit = family$linkinv(latentpart + 
                                                         regpart), family = family, deviance = l)
}

slice.glm = function(X, Y, slices, coefs, offset=NULL, penalized=0, parallel=1, family = poisson(), method="step", stepsize = 1){
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
    glm.basic(X[ok,,drop=FALSE],
              Y[ok,i,drop=FALSE],
              beta=coefs[,i],
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

gmf.initial = function(Y,X,d=min(dim(Y)),family = poisson()){
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



