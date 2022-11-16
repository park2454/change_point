unique_change2 = array(NA,c(6,6,2,30))

for(rep in 1:30){
  for(i in 1:2){
  change_dist = (i==1)
    for(m in 1:6){
      T = 10^m
      for(n in 1:6){
        p = 10^n
        if(change_dist){
          tau = rbinom(p,T-1,0.5) #binomial values 0 ~ (T-1)
        } else {
          tau = sample(1:T, p, replace=TRUE)
        }
        unique_change2[m,n,i,rep] = log10(length(unique(tau)))
        cat(n,m,l,j,i,rep,"\n")
      }
    }
  }
}

for(k in 1:4){
  apply(low_rank[,,k,1,1,],c(1,2),mean) |> print()
}

for(k in 1:4){
  apply(low_rank[,,k,2,1,],c(1,2),mean) |> print()
}


apply(low_rank[,,1,1,1,],c(1,2),mean)
apply(unique_change[,,1,],c(1,2),mean)
apply(unique_change2[,,1,],c(1,2),mean)

apply(low_rank[,,1,2,1,],c(1,2),mean)
apply(unique_change[,,2,],c(1,2),mean)
apply(unique_change2[,,2,],c(1,2),mean)

save(low_rank, nuclear_norm, unique_change, unique_change2, file="rank.RData")
