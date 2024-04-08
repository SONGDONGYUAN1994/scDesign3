## Part from Chenxin Jiang, thresholding for sparse matrix
### Thresholding Operator
thresh_op <- function(z, operator, delta, n){
  if(operator == 'hard'){
    s_method <- s_hard
  }else if(operator == 'soft'){
    s_method <- s_soft
  }else if(operator == 'scad'){
    s_method <- s_scad
  }else if(operator == 'al'){
    s_method <- s_al
  }else{
    stop('Please specify a valid thresholding operator.')
  }
  s_method(z, delta, n)
}

# Operator 1: Hard Thresholding
s_hard <- function(z, delta, n){
  
  p<-dim(z)[2]
  
  # if(is.null(theta)){
  #   theta=diag(p)
  # }
  # lambda = sqrt(theta*log(p)/n)*delta
  
  lambda <- sqrt(log(p)/n)*delta
  output <- (z>lambda)*z
  diag(output) <- diag(z)
  return(output)
}

# Operator 2: Soft Thresholding
s_soft <- function(z, delta, n){
  p<-dim(z)[1]
  lambda <- sqrt(log(p)/n)*delta
  z0 <- abs(z)-lambda
  output <- sign(z)*(z0>0)*z0
  diag(output) <- diag(z)
  return(output)
}

# Operator 3: SCAD (smoothly clipped absolute deviation)
s_scad <- function(z, delta, n, a=3.7){
  p<-dim(z)[1]
  lambda <- sqrt(log(p)/n)*delta
  
  output <- matrix(NA, dim(z)[1], dim(z)[2])
  
  index1 <- which(abs(z)<=2*lambda)
  z0 <- abs(z)-lambda
  output[index1] <- sign(z[index1])*(z0[index1]>0)*z0[index1]
  
  index2 <- which(abs(z)>2*lambda & abs(z)<=a*lambda)
  output[index2] <- ((a-1)*z[index2]-sign(z[index2])*a*lambda)/(a-2)
  
  index3 <- which(abs(z)>a*lambda)
  output[index3] <- z[index3]
  
  diag(output) <- diag(z)
  return(output)
}


# Operator 4: Adaptive lasso
s_al <- function(z, delta, n){
  p<-dim(z)[1]
  lambda <- sqrt(log(p)/n)*delta
  
  eta <- 3
  z0 <- abs(z) - lambda^(eta+1)*abs(z)^(-eta)
  output <- sign(z)*(z0>0)*z0
  diag(output) <- diag(z)
  return(output)
}

# Delta
### The function to select the optimal thresholding level
delta.est = function(data, 
                     method=c('cv', 'qiu'),
                     operator=c('hard', 'soft', 'scad', 'al')){
  n <- dim(data)[1]
  if((method=='qiu') && (operator=='hard')){
    s = covariance(data) *(n-1)/n
    delta = qiu.select(data, s)
  }else if(method=='cv'){
    delta = cv.min(data, operator)
  }else{
    stop('Please specify a valid thresholding method and an operator function.')
  }
  return(delta)
}


### Qiu function to tune delta
qiu.select = function(data, s=NULL){
  n = dim(data)[1]
  p = dim(data)[2]
  
  if(is.null(s)){
    s = covariance(data) *(n-1)/n
  }
  
  # standardized covariance of Sigma
  #theta = theta_est(data)
  #eta = sqrt(n/log(p))*s*theta^(-1/2)
  eta = sqrt(n/log(p))*s 
  # select lower triangular
  ltrig = abs(eta[lower.tri(eta, diag = FALSE)])
  
  # parameter to select optimal thr
  a = min(sqrt(2+log(n)/log(p)), 2)
  a1 = 2 - a
  a0 = (log(log(p)))^(-1/2)
  M = sum(((a1+a0)<ltrig) & (ltrig<=2))
  q = (p-1)*p/2
  plog = (log(p))^(1/2)
  V = 2*q*(pnorm(2*plog)-pnorm((a1+a0)*plog))
  N2 = max(M-V, plog)
  
  #N2 = sum((a1<=ltrig) & (ltrig<=2))
  delta = sqrt(2*(2-log(N2*1/plog)/log(p)))
  
  return(delta)
}

#### CV to tune delta
# The tuning parameter lambda was selected by minimizing the Frobenius norm of the 
# difference between s(sample_cov(train)) and sample_cov(test) by cross validation

# cv.select: For implementation, we choose the optimal delta from a delta.list 
# instead of really solving the optimization problem.
# lambda.list ranges in [0,lambda.max], with a user specified length.

# cv.min: solve the optimization problem directly

cv.select <- function(data, 
                      operator, 
                      fold=5, 
                      delta.length=10, 
                      delta.max=2,
                      random.select=FALSE,
                      random.size=500,
                      seed=123){
  
  n = dim(data)[1]
  p = dim(data)[2]
  
  if(p<100){
    random.select = FALSE
  }
  
  if(random.select){
    set.seed(seed)
    gene_idx = sample(1:p, size=random.size)
    data = data[,gene_idx]
  }
  
  
  # # Randomly shuffle the data
  # data0 <- data[sample(nrow(data)),]
  # sample.cov =  covariance(data, large = TRUE)
  
  # Create delta list
  delta.list = seq(0, delta.max, length.out=delta.length)
  
  # Create equal size folds
  folds <- cut(seq(1, nrow(data)), breaks=fold, labels=FALSE)
  
  fold_losses <- rep(0, fold) # Record fold loss 
  cv_errors <- rep(0, delta.length) # Record average estimated loss
  
  
  ## Perform k fold cross validation
  
  # iterate over the values of lambda
  for(j in seq_len(delta.length)){
    
    # iterate over CV folds
    for(i in seq_len(fold)){
      # Segment the data by fold 
      testIndexes <- which(folds==i, arr.ind=TRUE) # i-th fold as the testing data
      #trainIndexes <- which(folds!=i, arr.ind=TRUE) # i-th fold as the testing data
      testData <- data[testIndexes, ]
      trainData <- data[-testIndexes, ]
      # Get the covariance matrix estimator s based on the training data
      sample.cov.train = covariance(trainData, center=TRUE, large = TRUE)
      theta = theta_est(trainData)
      s = thresh_op(sample.cov.train, operator=operator, delta = delta.list[j], n=n)
      # Compute Frobenius risk = norm(distance of s and covariance of testing data)
      sample.cov.test = covariance(testData, center=TRUE, large = TRUE)
      fold_losses[i] = norm(s - sample.cov.test, type='F')
    }
    # Get the average estimated loss of CV
    cv_errors[j] <-  mean(fold_losses)
  }
  
  # Find the lambda that minimize cv_errors
  delta = delta.list[which.min(cv_errors)]
  
  delta
}

cv.min <- function(data, 
                   operator, 
                   fold=5,
                   delta.max=2,
                   random.select=FALSE,
                   random.size=100){
  n = dim(data)[1]
  p = dim(data)[2]
  
  if(p<100){
    random.select = FALSE
  }
  
  if(random.select){
    #set.seed(seed)
    gene_idx = sample(1:p, size=random.size)
    data = data[,gene_idx]
  }
  
  ## Perform k fold cross validation
  cv.loss = function(delta){
    # Create equal size folds
    folds <- cut(seq(1, nrow(data)), breaks=fold, labels=FALSE)
    fold_losses <- rep(0, fold) # Record fold loss 
    
    # iterate over CV folds
    for(i in seq_len(fold)){
      # Segment the data by fold 
      testIndexes <- which(folds==i, arr.ind=TRUE) # i-th fold as the testing data
      #trainIndexes <- which(folds!=i, arr.ind=TRUE) # i-th fold as the testing data
      testData <- data[testIndexes, ]
      trainData <- data[-testIndexes, ]
      # Get the covariance matrix estimator s based on the training data
      sample.cov.train = covariance(trainData, center=TRUE, large = TRUE)
      #theta = theta_est(trainData)
      #s = s_method(sample.cov.train, delta = delta, n=n)
      s = thresh_op(sample.cov.train, operator=operator, delta = delta, n=n)
      # Compute Frobenius risk = norm(distance of s and covariance of testing data)
      sample.cov.test = covariance(testData, center=TRUE, large = TRUE)
      fold_losses[i] = norm(s - sample.cov.test, type='F')
    }
    # Get the average estimated loss of CV
    cv_errors <-  mean(fold_losses)
    cv_errors
  }
  
  # Solve the optimization problem
  res = optimize(cv.loss, c(0, delta.max))
  delta = res$minimum
  delta
}

# variation of sample variance
theta_est <- function(data){
  n <- dim(data)[1]
  p <- dim(data)[2]
  
  theta <- matrix(0, p, p)
  for(i in 1:n){
    v <- data[i,] - colMeans(data)
    theta <- (Matrix::tcrossprod(v,v) - s)^2 + theta
  }
  theta <- theta/n
  theta
}

### Thresholding
est_threshold <- function(data, 
                          method = c('cv', 'qiu'),
                          operator = c('hard', 'soft', 'scad', 'al'),
                          corr = TRUE){
  p <- dim(data)[1]
  n <- dim(data)[2]
  
  # sample covariance
  z <- covariance(data) *(n-1)/n
  
  # select the optimal thresholding level
  delta.cv <- delta.est(data, method=method, operator=operator)
  s <- thresh_op(z, operator=operator, delta=delta.cv, n=n)
  
  # Modify s to make is psd
  tol <- 1e-6
  ev <- eigen(s, symmetric=TRUE, only.values = TRUE)$values
  s1 <- s + (tol-min(ev))*diag(dim(s)[1]) 
  
  if(corr){
    # make corr
    s1_corr <- stats::cov2cor(s1)
    output <- s1_corr
  }else{
    output <- s1
  }
  
  return(output)
}