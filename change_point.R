cpfm = function(Y, k,
                method = c("fixed", "em"), #default fixed
                parametric = TRUE,
                stop_criterion = ifelse(method=="fixed", 1e-6, 1e-3),
                max_iter = 1e3, 
                verbose = TRUE,
                identifiability = TRUE
                ){
  start = Sys.time()
  
  cusum = function(x){
    n = length(x)
    mu0 = cumsum(x)[-n]
    mu1 = sum(x) - mu0
    mu0 = mu0 / 1:(n-1)
    mu1 = mu1 / (n-1):1
    CUSUM = sqrt(1:(n-1) * (n-1):1)*abs(mu0-mu1)
    ind = which.max(CUSUM)
    out = c(mu0[ind], mu1[ind], ind)
    return(out)
  }
  
  create_M = function(change,n){
    p = nrow(change)
    change[,3] = floor(change[,3])
    M = matrix(change[,2], n, p, byrow = TRUE) # matrix of mu1's
    for(j in 1:p){
      M[1:change[j,3],j] = change[j,1] # replace mu0's according to tau's
    }
    return(M)
  }
  
  fixed_change = function(Y, k, max_iter, stop_criterion, verbose){
    
    n = nrow(Y)
    p = ncol(Y)
    
    # initialize
    Theta = matrix(Inf,n,k)
    A = matrix(Inf,p,k)
    change = matrix(0,p,3)
    M = matrix(0,n,p)
    
    for(iter in 1:max_iter){
      # svd
      svd_out = svd(Y-M, nu=k, nv=k)
      Theta_new = svd_out$d[1:k] * t(svd_out$u)
      Theta_new = t(Theta_new)
      # Theta_new = Theta_new - matrix(apply(Theta_new, 2, mean),n,k,byrow=TRUE)
      A_new = svd_out$v
      
      # cusum
      change_new = t( apply( Y- Theta_new %*% t(A_new), 2, cusum) )
      M_new = create_M(change_new,n)
      
      dt = norm(Theta-Theta_new, type="F") +
        norm(A-A_new, type="F") +
        norm(M-M_new, type="F")
      if(iter > 1 & dt < stop_criterion){
        if(verbose) cat("converged at iteration ", iter, "\n")
        break
      }
      Theta = Theta_new
      A = A_new
      change = change_new
      M = M_new
    }
    
    if(verbose & iter==max_iter) cat("reached max iteration",max_iter,"\n")
    
    out = list(
      Theta = Theta_new,
      A = A_new,
      change = change_new,
      M = M_new
    )
    return(out)
  }
  
  random_change = function(Y, k, parametric, max_iter, stop_criterion, verbose){
    
    n = nrow(Y)
    p = ncol(Y)
    
    # initialize
    mu0 = mu1 = apply(Y,2,mean)
    prior = rep(1/n,n)
    svd_out = svd(Y - matrix(mu1,n,p,byrow=TRUE), nu=k, nv=k)
    Theta = svd_out$d[1:k] * t(svd_out$u)
    Theta = t(Theta)
    A = svd_out$v
    Q = -Inf
    
    Q_trace = rep(NA,max_iter)
    for(iter in 1:max_iter){
      # E-step
      rss0 = ( Y - Theta %*% t(A) - matrix(mu0,n,p,byrow=TRUE) )^2
      rss0 = apply(rss0,2,cumsum)
      rss1 = ( Y - Theta %*% t(A) - matrix(mu1,n,p,byrow=TRUE) )^2
      rss1 = matrix(apply(rss1,2,sum),n,p,byrow=TRUE) - apply(rss1,2,cumsum)
      w_adj = -(rss0+rss1)
      w_adj = w_adj - matrix(apply(w_adj,2,max),n,p,byrow=TRUE)
      posterior = prior * exp(w_adj)
      posterior = posterior / matrix(apply(posterior,2,sum),n,p,byrow=TRUE)
      
      # M-step
      
      # Theta, A
      w_adj = apply( rbind(0,posterior[-n,]), 2, cumsum)
      Y_adj = Y- matrix(mu1,n,p,byrow=TRUE)*w_adj-matrix(mu0,n,p,byrow=TRUE)*(1-w_adj)
      svd_out = svd(Y_adj, nu=k, nv=k)
      Theta = svd_out$d[1:k] * t(svd_out$u)
      Theta = t(Theta)
      # Theta = Theta - matrix(apply(Theta, 2, mean),n,k,byrow=TRUE)
      A = svd_out$v
      
      # mu
      Y_adj = Y - Theta %*% t(A)
      mu0 = apply( apply(Y_adj,2,cumsum) * posterior, 2, sum )
      mu0 = mu0 / apply( (1:n)*posterior, 2, sum )
      mu1 = matrix(apply(Y_adj,2,sum),n,p,byrow=TRUE) - apply(Y_adj,2,cumsum)
      mu1 = apply(mu1*posterior,2,sum) / apply( (n:1-1)*posterior,2,sum )
      
      # pi
      if(parametric){
        alpha = sum((1:n-1)*posterior)/sum(posterior)/(n-1)
        prior = dbinom(1:n-1,n-1,alpha)
      } else {
        prior = apply(posterior,1,mean)
      }
      
      # check
      Q_new = sum((-(rss0 + rss1) + log(prior)) * posterior)
      Q_trace[iter] = Q_new
      
      if (Q_new - Q < stop_criterion){
        if (verbose) cat("converged at iteration",iter,"\n")
        break
      } else {
        Q = Q_new
      }
    }
    
    if(verbose & iter==max_iter) cat("reached max iteration",max_iter,"\n")
    
    change = cbind( mu0, mu1, apply(posterior*(1:n),2,sum) )
    M = create_M(change,n)
    out = list(
      Theta = Theta,
      A = A,
      change = change,
      M = M,
      prior = prior,
      Q_trace = Q_trace
    )
    return(out)
  }
  
  method = method[1]
  
  if(method == "fixed"){
    out = fixed_change(Y, k, max_iter, stop_criterion, verbose)
  } else if(method == "em"){
    out = random_change(Y, k, parametric, max_iter, stop_criterion, verbose)
  } else {
    stop("method requires either \"fixed\" or \"em\"")
  }
  
  if(identifiability){
    n = nrow(out$Theta)
    k = ncol(out$Theta)
    Tbar = apply(out$Theta, 2, mean)
    out$Theta = out$Theta - matrix(Tbar,n,k,byrow=TRUE)
    out$M = out$M + matrix(Tbar,n,k,byrow=TRUE) %*% t(out$A)
    out$change[,1:2] = out$change[,1:2] + out$A %*% cbind(Tbar, Tbar)
  }
  
  if(verbose) print(Sys.time() - start)
  return(out)
}
