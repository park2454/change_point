lowrank = function(T, p, k,
                 max_iter = 1000, 
                 verbose = TRUE, 
                 SVD = TRUE, 
                 signal_mag = 3,
                 rand_init = FALSE){
  
  # data generation
  Theta = matrix(rnorm(T*k),T,k)
  A = matrix(rnorm(k*p),p,k)
  
  m = rnorm(p,signal_mag)
  M = outer(rep(1,T),m)
  
  Y = Theta %*% t(A) +
    M +
    matrix(rnorm(T*p), T, p)
  
  # initialize
  if(rand_init){
    Theta_hat = matrix(0,T,k)
    A_hat = matrix(0,p,k)
    m_hat = rep(0,p)
  } else {
    Theta_hat = Theta
    A_hat = A
    m_hat = m
  }
  
  # loop
  for(iter in 1:max_iter){
    # adjusted mean
    adj_Y = Y - Theta_hat %*% t(A_hat)
    
    # CUSUM
    m_new = apply(adj_Y,2,mean)
    
    # adjusted mean
    adj_Y = Y - outer(rep(1,T), m_new)
    
    if(SVD){
      # SVD
      svd_hat = svd(adj_Y, nu=k, nv=k)
      Theta_new = svd_hat$d[1:k] * t(svd_hat$u)
      Theta_new = t(Theta_new)
      A_new = svd_hat$v
    } else {
      # lm for theta, A
      mlr = lm(t(adj_Y)~A_hat+0)
      Theta_new = t(mlr$coefficients)
      
      mlr = lm(adj_Y~Theta_new+0)
      A_new = t(mlr$coefficients)
    }
    
    # check
    mse = function(x) mean(x^2, na.rm=TRUE)
    dt = mse(m_new - m_hat) + mse(Theta_new%*%t(A_new) - Theta_hat%*%t(A_hat))
    
    if (dt < 1e-6){
      if (verbose) cat("converged at iteration",iter,"\n")
      break
    } else {
      Theta_hat = Theta_new
      A_hat = A_new
      m_hat = m_new
    }
  }
  if(iter==max_iter){
    if (verbose) cat("reached reached max iteration",max_iter,"\n")
  } 
  
  out = list(
    true = list(
      Theta = Theta,
      A = A,
      factor = Theta %*% t(A),
      M = outer(rep(1,T),m),
      m = m,
      Y = Theta %*% t(A) + outer(rep(1,T),m)
    ),
    estimate = list(
      Theta = Theta_new,
      A = A_new,
      factor = Theta_new %*% t(A_new),
      M = outer(rep(1,T),m_new),
      m = m_new,
      Y = Theta_new %*% t(A_new) + outer(rep(1,T),m_new)
    ),
    dim = list(
      T = T,
      p = p,
      k = k
    )
  )
  return(out)
}

lowrank0 = function(T, p, k){
  
  # data generation
  Theta = matrix(rnorm(T*k),T,k)
  A = matrix(rnorm(k*p),p,k)
  
  Y = Theta %*% t(A) + matrix(rnorm(T*p), T, p)
  
  svd_hat = svd(Y, nu=k, nv=k)
  Theta_new = svd_hat$d[1:k] * t(svd_hat$u)
  Theta_new = t(Theta_new)
  A_new = svd_hat$v
  
  out = list(
    true = list(
      Theta = Theta,
      A = A,
      Y = Theta %*% t(A)
    ),
    estimate = list(
      Theta = Theta_new,
      A = A_new,
      Y = Theta_new %*% t(A_new)
    ),
    dim = list(
      T = T,
      p = p,
      k = k
    )
  )
  return(out)
}

mse_factor = array(NA,dim = c(5,5,4,10))
mse_M = array(NA,dim = c(5,5,4,10))
mse_Y = array(NA,dim = c(5,5,4,10))
mse_Y2 = array(NA,dim = c(5,5,4,10))
prior = array(NA,dim = c(5,5,4,10))

distr=TRUE
param=TRUE
rand=TRUE
mag=3

mse = function(x) mean(x^2, na.rm=TRUE)

K = c(2,5,10,50)
for(rep in 1:10){
  for(k in 1:3){
    for(T in 1:3){
      for(p in 1:3){
        out = lowrank(T=T*200,p=p*200,k=K[k], rand_init=TRUE)
        mse_factor[T,p,k,rep] = mse(out$true$factor - out$est$factor)
        mse_M[T,p,k,rep] = mse(out$true$M - out$est$M)
        mse_Y[T,p,k,rep] = mse(out$true$Y - out$est$Y)
        out = lowrank0(T=T*200,p=p*200,k=K[k])
        mse_Y2[T,p,k,rep] = mse(out$true$Y - out$est$Y)
        save(mse_factor, mse_M, mse_Y, mse_Y2,
             file = "baseline.rda")
      }
    }
  }
}
