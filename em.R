em = function(T, p, k,
              max_iter = 1e3, 
              verbose = TRUE, 
              SVD = TRUE, 
              change_dist = FALSE, #FALSE: uniform, TRUE: binomial,
              parametric = FALSE,
              prob = 0.5,
              signal_mag = 3,
              rand_init = FALSE){
  # data generation
  Theta = matrix(rnorm(T*k),T,k)
  A = matrix(rnorm(k*p),p,k)
  
  adjust = function(x){
    out = c(rep(x[1], x[3]), rep(x[2], T-x[3]))
    return(out) 
  }
  
  pre = rnorm(p)
  post = rnorm(p, 2*signal_mag*rbinom(p,1,0.5)-signal_mag)
  if(change_dist){
    tau = 1+rbinom(p,T-1,prob) #binomial values 0 ~ (T-1)
    prior = dbinom(1:T-1,T-1,prob)
  } else {
    tau = sample(1:T, p, replace=TRUE)
    prior = rep(1/T,T)
  }
  change = cbind(pre,post,tau)
  
  Y = Theta %*% t(A) +
    apply(change,1,adjust) +
    matrix(rnorm(T*p), T, p)
  
  # initialize
  if(rand_init){
    pre_hat = apply(Y,2,mean)
    post_hat = pre_hat
    prior_hat = rep(1/T,T)
    svd_hat = svd(Y - outer(rep(1,T), post_hat), nu=k, nv=k)
    Theta_hat = svd_hat$d[1:k] * t(svd_hat$u)
    Theta_hat = t(Theta_hat)
    A_hat = svd_hat$v
  } else {
    Theta_hat = Theta
    A_hat = A
    pre_hat = pre
    post_hat = post
    prior_hat = prior
  }
  Q_hat = -Inf
  
  for(iter in 1:max_iter){
    # E-step
    rss_pre = ( Y - Theta_hat %*% t(A_hat) - outer( rep(1,T), pre_hat) )^2
    rss_pre = apply(rss_pre,2,cumsum)
    rss_post = ( Y - Theta_hat %*% t(A_hat) - outer( rep(1,T), post_hat) )^2
    rss_post = outer( rep(1,T), apply(rss_post,2,sum) ) - apply(rss_post,2,cumsum)
    w_adj = -(rss_pre+rss_post)
    w_adj = w_adj - outer(rep(1,T), apply(w_adj,2,max))
    posterior = prior_hat * exp(w_adj)
    posterior = posterior / outer( rep(1,T),  apply(posterior,2,sum) )
    
    # M-step
    
    # Theta, A
    w_adj = apply( rbind(0,posterior[-T,]), 2, cumsum)
    Y_adj = Y-outer(rep(1,T), post_hat)*w_adj-outer(rep(1,T), pre_hat)*(1-w_adj)
    if(SVD){
      # SVD
      svd_hat = svd(Y_adj, nu=k, nv=k)
      Theta_new = svd_hat$d[1:k] * t(svd_hat$u)
      Theta_new = t(Theta_new)
      A_new = svd_hat$v
    } else {
      # multiple linear regression
      mlr = lm(t(Y_adj)~A_hat+0)
      Theta_new = t(mlr$coefficients)
      colnames(Theta_new) = NULL
      
      mlr = lm(Y_adj~Theta_new+0)
      A_new = t(mlr$coefficients)
      colnames(A_new) = NULL
    }
    
    # mu
    Y_adj = Y - Theta_new %*% t(A_new)
    pre_new = apply( apply(Y_adj,2,cumsum) * posterior, 2, sum )
    pre_new = pre_new / apply( (1:T)*posterior, 2, sum )
    post_new = outer( rep(1,T), apply(Y_adj,2,sum) ) - apply(Y_adj,2,cumsum)
    post_new = apply(post_new*posterior,2,sum)/apply( (T:1-1)*posterior,2,sum )
    
    # pi
    if(parametric){
      alpha = sum((1:T-1)*posterior)/sum(posterior)/(T-1)
      prior_new = dbinom(1:T-1,T-1,alpha)
    } else {
      prior_new = apply(posterior,1,mean)
    }
    
    # check
    Q_new = sum((-(rss_pre + rss_post) + log(prior)) * posterior)
    
    if (Q_new - Q_hat < 1e-6){
      if (verbose) cat("converged at iteration",iter,"\n")
      break
    } else {
      Theta_hat = Theta_new
      A_hat = A_new
      pre_hat = pre_new
      post_hat = post_new
      prior_hat = prior_new
      Q_hat = Q_new
    }
  }
  if(iter==max_iter){
    if (verbose){
      cat(Q_new-Q_hat,"\n")
      cat("reached reached max iteration",max_iter,"\n")
    } 
  } 
  
  change_new = cbind(pre_new, post_new, apply(posterior,2,which.max))
  
  out = list(
    true = list(
      Theta = Theta,
      A = A,
      factor = Theta %*% t(A),
      M = apply(change,1,adjust),
      change = change,
      prior = prior,
      Y = Theta %*% t(A) + apply(change,1,adjust)
    ),
    estimate = list(
      Theta = Theta_new,
      A = A_new,
      factor = Theta_new %*% t(A_new),
      M = apply(change_new,1,adjust),
      change = change_new,
      prior = prior_new,
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

mse = function(x) mean(x^2, na.rm=TRUE)
tv = function(x) sum(abs(x), na.rm=TRUE)/2
mse_factor = array(NA,dim = c(5,5,4,10))
mse_M = array(NA,dim = c(5,5,4,10))
mse_Y = array(NA,dim = c(5,5,4,10))
prior = array(NA,dim = c(5,5,4,10))

distr=TRUE
param=TRUE
rand=TRUE
mag=3

K = c(2,5,10,50)
for(rep in 1:10){
  for(k in 1:3){
    for(T in 1:3){
      for(p in 1:3){
        out = em(T=T*200, p=p*200, k=K[k],
                 change_dist = distr, 
                 parametric = param,
                 rand_init = rand,
                 signal_mag = mag)
        mse_factor[T,p,k,rep] = mse(out$true$factor - out$est$factor)
        mse_M[T,p,k,rep] = mse(out$true$M - out$est$M)
        mse_Y[T,p,k,rep] = mse(out$true$Y - out$est$Y)
        prior[T,p,k,rep] = tv(out$true$prior - out$est$prior)
        save(mse_factor, mse_M, mse_Y, prior,
             file = paste0("em_dist=",distr,"_mag=",mag,"_rand=",rand,"_param=",param,".rda"))
      }
    }
  }
}