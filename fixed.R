fixed = function(T, p, k,
                 max_iter = 1000, 
                 verbose = TRUE, 
                 SVD = TRUE, 
                 change_dist = FALSE,
                 prob = 0.5,
                 signal_mag = 3,
                 rand_init = FALSE,
                 identifiability = TRUE){
  
  # data generation
  Theta = matrix(rnorm(T*k),T,k)
  A = matrix(rnorm(k*p),p,k)
  
  adjust = function(x){
    out = c(rep(x[1], x[3]), rep(x[2], T-x[3]))
    return(out) 
  }
  
  pre_mean = rnorm(p)
  post_mean = rnorm(p, 2*signal_mag*rbinom(p,1,0.5)-signal_mag)
  
  if(identifiability){
    # identifiability
    th_bar = apply(Theta, 2, mean)
    Theta = Theta - outer( rep(1,T), th_bar)
    pre_mean = pre_mean + A %*% th_bar
    post_mean = post_mean + A %*% th_bar
  }
  
  if(change_dist){
    tau = 1+rbinom(p,T-1,prob) #binomial values 0 ~ (T-1)
  } else {
    tau = sample(1:T, p, replace=TRUE)
  }
  
  change = cbind(pre_mean,post_mean,tau)
  
  Y = Theta %*% t(A) +
    apply(change,1,adjust) +
    matrix(rnorm(T*p), T, p)
  
  # CUSUM function
  # input: a vector
  # output: pre-, post-mean, change point
  cusum = function(x){
    T = length(x)
    res = matrix(0,T-1,3)
    for (t in 2:T-1){
      res[t,1] = mean(x[1:t])
      res[t,2] = mean(x[-(1:t)])
      res[t,3] = sqrt(t*(T-t))*abs(res[t,1] - res[t,2])
    }
    sol = which.max(res[,3])
    out = c(res[sol,1:2], sol)
    return(out)
  }
  
  # initialize
  if(rand_init){
    Theta_hat = matrix(0,T,k)
    A_hat = matrix(0,p,k)
    change_hat = apply(Y,2,mean)
  } else {
    Theta_hat = Theta
    A_hat = A
    change_hat = change
  }
  
  # loop
  for(iter in 1:max_iter){
    # adjusted mean
    adj_Y = Y - Theta_hat %*% t(A_hat)
    
    # CUSUM
    change_new = t(apply(adj_Y,2,cusum))
    
    # adjusted mean
    adj_Y = Y - apply(change_new,1,adjust)
    
    if(SVD){
      # SVD
      svd_hat = svd(adj_Y, nu=k, nv=k)
      Theta_new = svd_hat$d[1:k] * t(svd_hat$u)
      Theta_new = t(Theta_new)
      if(identifiability){
        Theta_new = Theta_new - outer( rep(1,T), apply(Theta_new, 2, mean) )
      }
      A_new = svd_hat$v
    } else {
      # lm for theta, A
      mlr = lm(t(adj_Y)~A_hat+0)
      Theta_new = t(mlr$coefficients)
      
      mlr = lm(adj_Y~Theta_new+0)
      A_new = t(mlr$coefficients)
    }
    
    # check
    dt = norm(change_new - change_hat, type="F")^2 +
      norm(Theta_new%*%t(A_new) - Theta_hat%*%t(A_hat), type="F")^2
    
    if (dt < 1e-6){
      if (verbose) cat("converged at iteration",iter,"\n")
      break
    } else {
      Theta_hat = Theta_new
      A_hat = A_new
      change_hat = change_new
    }
  }
  if(iter==max_iter){
    if (verbose) cat("reached reached max iteration",max_iter,"\n")
  } 
  
  out = list(
    true = list(
      Theta = Theta,
      A = A,
      factor = Theta %*% t(A) + outer(rep(1,T), change[,1]),
      M = apply(change,1,adjust),
      change = change,
      Y = Theta %*% t(A) + apply(change,1,adjust)
    ),
    estimate = list(
      Theta = Theta_new,
      A = A_new,
      factor = Theta_new %*% t(A_new) + outer(rep(1,T), change_new[,1]),
      M = apply(change_new,1,adjust),
      change = change_new,
      Y = Theta_new %*% t(A_new) + apply(change_new,1,adjust)
    ),
    dim = list(
      T = T,
      p = p,
      k = k
    )
  )
  return(out)
}

# # estimation comparison
# pdf("histogram.pdf",width=12,height=6)
# par(mfrow=c(2,4))
# for(dist in c(FALSE,TRUE)){
#   for(mag in 0:3){
#     samp = c()
#     for(rep in 1:1){
#       out = fixed(T=100, p=100, k=5, change=dist, signal=mag)
#       err = (out$true$change[,3] - out$estimate$change[,3])/out$dim$T
#       samp = c(samp,err)
#     }
#     hist(samp, breaks=seq(-1,1,1/20), 
#          main = paste0("dist,mag=",dist,mag),
#          xlab = paste0("sample var=",round(var(samp),4)))
#   }
# }
# dev.off()
# 
# mean((out$true$Y - out$estimate$Y)^2)

# Most of the estimation is correct
# There are a few extreme errors but the latent part seems to be consistent
# Errors may only occur in the case where pre- and post-mean are similar.

# Use k = 2,5,10,50
# plot a heatmap of average errors for varying T and p

mse = function(x) mean(x^2, na.rm=TRUE)
tv = function(x) sum(abs(x), na.rm=TRUE)/2
mse_factor = array(NA,dim = c(5,5,4,10))
mse_M = array(NA,dim = c(5,5,4,10))
mse_Y = array(NA,dim = c(5,5,4,10))

distr=TRUE
rand=TRUE
mag=3

K = c(2,5,10,50)
for(rep in 1:10){
  for(k in 1:3){
    for(T in 1:3){
      for(p in 1:3){
        out = fixed(T=T*200, p=p*200, k=K[k],
                    change_dist = distr, 
                    rand_init = rand,
                    signal_mag = mag)
        mse_factor[T,p,k,rep] = mse(out$true$factor - out$est$factor)
        mse_M[T,p,k,rep] = mse(out$true$M - out$est$M)
        mse_Y[T,p,k,rep] = mse(out$true$Y - out$est$Y)
        save(mse_factor, mse_M, mse_Y,
             file = paste0("fixed_dist=",distr,"_mag=",mag,"_rand=",rand,".rda"))
      }
    }
  }
}
