## STANDARD FUNCTIONS
## Function to draw n columns of lenght 'size' of independent N(0,1)
inorms <- function(n,size)
  {
    one.matrix <- matrix(rnorm(n*size),ncol=n,nrow=size)
    return(one.matrix)
  }

## Function to Draw from a Matric Variate Normal
## MATVN(Mu, S (kron) L)
## L is a k x k matrix, k=the number of rows of Mu
## S is a n x n matrix, n=the number of columns of Mu
matnorm <- function (S,L,Mu) 
{ 
 Mu <- as.matrix(Mu)
 S.chol <- chol(S)
 L.chol <- chol(L)
 k <- nrow(L)
 j <- nrow(S)
 
 U <- inorms(j,k)
 #Z <- L.chol %*% U %*% t(S.chol) + Mu
 Z <- t(L.chol) %*% U %*% S.chol + Mu  
 Z <- as.matrix(Z)
 colnames(Z) <- NULL
 return(Z)
}

## Function to Draw from a I Wishart
## IW(p,V), from Richard Paap lecture week 7
iwishart <- function (v,S) 
{ 
  v <- round(v)
  U <- inorms(v,ncol(S))  # inorms() is a function defined below
  S.chol <- chol(S)
  M <- solve(U %*% t(U))
  Z.IW <- t(S.chol) %*% M %*% S.chol
  return(Z.IW)
}



f.bayes.lr <- function (Y,X,B_bar,B_COV_P,sigma2,sig_df,sig_scale) 
{
 # Y = X * B + E  where E ~ N(0,sigma2*I)
 # B_COV_P   = prior covariance matrix for B
 # B_bar     = prior mean of B
 # sigma2    = variance of E
 # sig_df    = prior degrees of freedom for sigma2
 # sig_scale = prior scale for sigma2
 
 ## this function draws the betas and sigma using standard conjugate piors
 ## on the linear regression model

 ## draw for beta
 xxmat <- t(X) %*% X + qr.solve(B_COV_P)
 beta.part1 <- qr.solve(xxmat) %*% ( t(X) %*% yvec + qr.solve(B_COV_P) %*% B_bar ) 
 beta.part2 <- chol(sigma2 * qr.solve(xxmat)) %*% rnorm(ncol(X))
 B_tilde <- beta.part1 + beta.part2
 rownames(B_tilde) <- NULL

 ## draw for sigma
 resid.part1 <- t(yvec - X %*% B_tilde) %*% (yvec - X %*% B_tilde)
 #resid.part2 <- t(B_tilde - B_bar) %*% solve(B_COV_P) %*% (B_tilde - B_bar)
 resid.part3 <- sig_df * sig_scale
 resid2 <- resid.part1 + resid.part3
 nvec <- rnorm(length(yvec)+sig_df)
 chi2 <- t(nvec) %*% nvec 
 sigma2_tilde <- resid2 / chi2 

 return(list(beta=B_tilde,sigma2=sigma2_tilde))
}


rdirichlet <- function (alpha) 
{ 
dim = length(alpha)
y = rep(0, dim)
for (i in 1:dim) y[i] = rgamma(1, alpha[i])
return(y/sum(y))
}
