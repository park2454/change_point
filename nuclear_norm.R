low_rank1 = array(NA,c(5,5,4,2,2,30))
nuclear_norm1 = array(NA,c(5,5,4,2,2,30))
low_rank2 = array(NA,c(5,5,4,2,2,30))
nuclear_norm2 = array(NA,c(5,5,4,2,2,30))
low_rank3 = array(NA,c(5,5,4,2,2,30))
nuclear_norm3 = array(NA,c(5,5,4,2,2,30))

adjust = function(x){
  out = c(rep(x[1], x[3]), rep(x[2], T-x[3]))
  return(out) 
}

K = c(2,5,10,50)
for(rep in 1:30){
  for(i in 1:2){
    change_dist = (i==1)
    for(j in 1:2){
      signal_mag = ifelse(j==1,3,0)
      for(l in 1:4){
        k = K[l]
        for(m in 1:5){
          T = 200*m
          for(n in 1:5){
            p = 200*n
            Theta = matrix(rnorm(T*k),T,k)
            A = matrix(rnorm(k*p),p,k)
            pre_mean = rnorm(p)
            post_mean = rnorm(p, 2*signal_mag*rbinom(p,1,0.5)-signal_mag)
            if(change_dist){
              tau = 1+rbinom(p,T-1,0.5) #binomial values 0 ~ (T-1)
            } else {
              tau = sample(1:T, p, replace=TRUE)
            }
            change = cbind(pre_mean,post_mean,tau)
            Y = Theta %*% t(A)
            SVD = svd(Y, nu=0, nv=0)
            low_rank1[m,n,l,i,j,rep] = sum(SVD$d>1e-6)
            nuclear_norm1[m,n,l,i,j,rep] = sum(SVD$d)
            Y = Y + apply(change,1,adjust)
            SVD = svd(Y, nu=0, nv=0)
            low_rank2[m,n,l,i,j,rep] = sum(SVD$d>1e-6)
            nuclear_norm2[m,n,l,i,j,rep] = sum(SVD$d)
            Y = apply(change,1,adjust)
            SVD = svd(Y, nu=0, nv=0)
            low_rank3[m,n,l,i,j,rep] = sum(SVD$d>1e-6)
            nuclear_norm3[m,n,l,i,j,rep] = sum(SVD$d)
            cat(n,m,l,j,i,rep,"\n")
          }
        }
      }
    }
  }
}

setwd("~/research/change_point")
save(low_rank1,low_rank2,low_rank3,
     nuclear_norm1,nuclear_norm2,nuclear_norm3, 
     file = "nuclear_norm.RData")
load("nuclear_norm.RData")

# summation of factor and mean-change
apply(nuclear_norm2[,,,2,2,],c(1,2,3),mean)
# factor
apply(nuclear_norm1[,,,2,2,],c(1,2,3),mean)
# mean-change
apply(nuclear_norm3[,,,2,2,],c(1,2,3),mean)


# summation of factor and mean-change
apply(nuclear_norm2[,,,2,1,],c(1,2,3),mean)
# factor
apply(nuclear_norm1[,,,2,1,],c(1,2,3),mean)
# mean-change
apply(nuclear_norm3[,,,2,1,],c(1,2,3),mean)


M = apply(nuclear_norm3[,,,2,1,],c(1,2),mean)
M / M[1,1]
