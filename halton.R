
primer=function(v){
 return(regexpr("^1$|^(11+?)\\1+$",unlist(lapply(v,function(z){paste(rep("1",z),sep='',collapse='')})),perl=TRUE)
== -1)
}


halton <- function (nodraws,dim) 
{ 

  primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
              103,107,109,113,127,131,137,139, 149,151, 157, 163, 167, 173,179, 181,191, 193,197,199, 211, 223,227,229)
   
  discard <- 2*primes[dim]
  draws <- rep(0,nodraws)
  for (k in 1:dim)
  {
  s <- 0; n <- 0; i <- 1
  
  while(length(n)<(nodraws+discard+1))
    {
      n <- s
      for(j in 1:(primes[k]-1))
        {
        add <- s + j/primes[k]^i
        n <- c(n,add)
        }
      s <- n
      i <- i+1      
    }
  draws <- rbind(draws,s[(1+discard):(nodraws+discard)])
  
   }
  draws <- draws[-1,]
  return(qnorm(draws))
}


halton.s <- function (nodraws,dim) 
{ 

  primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
              103,107,109,113,127,131,137,139, 149,151, 157, 163, 167, 173,179, 181,191, 193,197,199, 211, 223,227,229)
   
  discard <- 2*primes[dim]
  draws <- rep(0,nodraws)
  for (k in 1:dim)
  {
  s <- 0; n <- 0; i <- 1
  oseq <- 1:(primes[k]-1)
  oseq <- sample(oseq,length(oseq))
  while(length(n)<(nodraws+discard+1))
    {
      n <- s
      for(j in 1:length(oseq))
        {
        add <- s + oseq[j]/primes[k]^i
        n <- c(n,add)
        }
      s <- n
      i <- i+1      
    }
  draws <- rbind(draws,s[(1+discard):(nodraws+discard)])
  
   }
  draws <- draws[-1,]
  return(list(draws=qnorm(draws),halton=draws))
}

halton.u <- function (nodraws) 
{ 
   out1 <- halton.s(nodraws,1)$halton + runif(1)
   out1[out1>1] <- out1[out1>1]-1
   return(qnorm(out1))
}
