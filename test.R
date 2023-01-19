source("change_point.R")

set.seed(1)
N=1000
K=3
T=500

res1 = matrix(NA,100,3)
colnames(res1) = c("overall", "factor", "M")
res2 = res1

for(iter in 1:100){
  # data generation
  A <- matrix(runif(N*K, 0.5, 1.5), N,K)
  F <- matrix(rnorm(T*K), T,K)
  F[(T/2+1):T,] = F[(T/2+1):T,] +1
  mu0 = runif(N, 0, 1) 
  mu1 = runif(N, 3, 4)
  tau.vec <- sample((T/10+1):(T*0.9), N, replace = TRUE)
  # tau.vec <- rbinom(N, T-1, prob = 0.5) + 1
  T.matr <- matrix(0, N, T)
  for(i in 1:N){
    T.matr[i, tau.vec[i]:T] = 1
  }
  data = A%*% t(F) + (1-T.matr) * mu0 + T.matr * mu1 + rnorm(N*T)
  
  fit = cpfm(t(data), K, identifiability = FALSE)
  
  # normalize M
  factor_true = t(A%*% t(F) )
  M_true = t((1-T.matr) * mu0 + T.matr * mu1)
  factor_true = factor_true + matrix(M_true[1,],T,N,byrow=TRUE)
  M_true = M_true - matrix(M_true[1,],T,N,byrow=TRUE)
  
  factor_est = fit$Theta %*% t(fit$A)
  M_est = fit$M
  factor_est = factor_est + matrix(M_est[1,],T,N,byrow=TRUE)
  M_est = M_est - matrix(M_est[1,],T,N,byrow=TRUE)
  
  res1[iter,1] = norm(factor_true + M_true - factor_est - M_est, type="F")^2
  res1[iter,2] = norm(factor_true - factor_est, type="F")^2
  res1[iter,3] = norm(M_true - M_est, type="F")^2
  
  # normalize Theta
  factor_true = t(A%*% t(F) )
  M_true = t((1-T.matr) * mu0 + T.matr * mu1)
  M_true = M_true + matrix(apply(F,2,mean),T,K,byrow=TRUE) %*% t(A)
  factor_true = (F - matrix(apply(F,2,mean),T,K,byrow=TRUE) ) %*% t(A) 
  
  factor_est = fit$Theta %*% t(fit$A)
  M_est = fit$M
  M_est = M_est + matrix(apply(fit$Theta,2,mean),T,K,byrow=TRUE) %*% t(fit$A)
  factor_est =  (fit$Theta - matrix(apply(fit$Theta,2,mean),T,K,byrow=TRUE) ) %*% t(fit$A) 
  
  res2[iter,1] = norm(factor_true + M_true - factor_est - M_est, type="F")^2
  res2[iter,2] = norm(factor_true - factor_est, type="F")^2
  res2[iter,3] = norm(M_true - M_est, type="F")^2
}

apply(res1, 2, function(x) median(x, na.rm=TRUE))
apply(res2, 2, function(x) median(x, na.rm=TRUE))
