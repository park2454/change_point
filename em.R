em = function(T, p, k,
              max_iter = 1e3, 
              verbose = TRUE, 
              SVD = TRUE, 
              change_dist = FALSE, #FALSE: uniform, TRUE: binomial,
              structured = FALSE,
              prob = 0.5,
              signal_mag = 3){
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
  Theta_hat = Theta
  A_hat = A
  pre_hat = pre
  post_hat = post
  prior_hat = prior
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
    if(structured){
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
  
  out = list(
    true = list(
      Theta = Theta,
      A = A,
      prior = prior
    ),
    estimate = list(
      Theta = Theta_new,
      A = A_new,
      prior = prior_new
    ),
    dim = list(
      T = T,
      p = p,
      k = k
    )
  )
  return(out)
}

# par(mfrow=c(2,1))
# out = em(T=1000,p=1000,k=50, structured = FALSE)
# plot(out$estimate$prior, type="l",col="red")
# lines(1:T,dbinom(1:T-1,T-1,0.5))
# out = em(T=1000,p=1000,k=50, structured = TRUE)
# plot(out$estimate$prior, type="l",col="red")
# lines(1:T,dbinom(1:T-1,T-1,0.5))

heat_factor = array(NA,dim = c(5,5,4,30))
heat_prior = array(NA,dim = c(5,5,4,30))
K = c(2,5,10,50)
for(rep in 1:30){
  for(k in 1:4){
    for(T in 1:5){
      for(p in 1:5){
        out = em(T=T*200, p=p*200, k=K[k])
        err = (out$true$Theta)%*%t(out$true$A)-(out$estimate$Theta)%*%t(out$estimate$A)
        heat_factor[T,p,k,rep] = sqrt(mean(err^2))
        err = sum(abs(out$true$prior - out$estimate$prior))/2
        heat_prior[T,p,k,rep] = err
        save(heat_factor, heat_prior, file = "em.RData")
      }
    }
  }
}
