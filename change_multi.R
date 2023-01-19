cpfm_multi = function(Y, k,
                method = c("fixed", "em"), #default fixed
                parametric = TRUE,
                stop_criterion = ifelse(method=="fixed", 1e-6, 1e-3),
                gamma = 0.1*log(nrow(Y)),
                delta = round(nrow(Y)*0.1),
                max_iter = 1e3, 
                verbose = TRUE,
                identifiability = TRUE
){
  start = Sys.time()
  
  fixed_change = function(Y, k, max_iter, stop_criterion, gamma_set, delta, verbose){
    
    n = nrow(Y)
    p = ncol(Y)
    
    # initialize
    Theta = matrix(Inf,n,k)
    A = matrix(Inf,p,k)
    M = matrix(0,n,p)
    M_new = matrix(0,n,p)
    # change_point = list()
    # change_mean = list()
    
    for(iter in 1:max_iter){
      # svd
      svd_out = svd(Y-M, nu=k, nv=k)
      Theta_new = svd_out$d[1:k] * t(svd_out$u)
      Theta_new = t(Theta_new)
      A_new = svd_out$v
      
      # change points
      Y_adj = Y - Theta_new %*% t(A_new)
      
      for(j in 1:p){
        DP_result = DP.univar(Y_adj[,j],gamma,delta)
        cpt = DP_result$cpt |> floor()
        g_mean = c(0,cumsum(Y_adj[,j])[c(cpt,n)]) |> diff()
        g_mean = g_mean/diff(c(0,cpt,n))
        M_new[,j] =  rep(g_mean, diff(c(0,cpt,n)) )
        # change_point$j = cpt
        # change_mean$j = g_mean
      }
      
      dt = norm(Theta-Theta_new, type="F") +
        norm(A-A_new, type="F") +
        norm(M-M_new, type="F")
      if(iter > 1 & dt < stop_criterion){
        if(verbose) cat("converged at iteration ", iter, "\n")
        break
      }
      Theta = Theta_new
      A = A_new
      M = M_new
      cat("iteration ", iter,"\n")
    }
    
    if(verbose & iter==max_iter) cat("reached max iteration",max_iter,"\n")
    
    out = list(
      Theta = Theta_new,
      A = A_new,
      M = M_new
      # change_point = change_point,
      # change_mean = change_mean
    )
    return(out)
  }
  
  method = method[1]
  
  if(method == "fixed"){
    out = fixed_change(Y, k, max_iter, stop_criterion, gamma_set, delta, verbose)
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
    # for(j in 1:p){
    #   out$change_mean = out$change_mean + A[j,] %*% Tbar
    # }
  }
  
  if(verbose) print(Sys.time() - start)
  return(out)
}
