
#sink(file="logNewR.txt" ,split="TRUE")
rm(list=ls())
load("inits.Rdata")
seed1 <- as.numeric(Sys.time())
set.seed(seed1)
##########
## DATA SECTION
load("datahazard.Rdata")
source("halton.R")                ## loads halton functions for random number generation
source("std_fun.R")


## data formating min(t,T) in columns
data.all[,1][is.na(data.all$t_price)] <- data.all$cens_t[is.na(data.all$t_price)]
data.all[,1][data.all$t_price>data.all$cens_t] <- data.all$cens_t[data.all$t_price>data.all$cens_t]
data.all[,3][is.na(data.all[,3])] <- data.all$cens_t[is.na(data.all[,3])]
data.all[,4][is.na(data.all[,4])] <- 0
data.all[,4] <- data.all$g_saddle *  data.all$g_recov  * data.all$d_saddle + data.all$g_saddle * ( 1 - data.all$g_recov ) * ( data.all$cens_t - data.all$t_saddle )

## three cases corrected
data.all[1158,4] <- 0; data.all[1158,3] <- 14; data.all[1158,5] <- 0
data.all[1197,4] <- 0; data.all[1197,3] <- 11; data.all[1197,5] <- 0
data.all[682,2] <- 0

## Data Matrix with t^tilde
t_tilde <- data.all[,c(1,3,4)]

## Regressors Data Matrix
Z <- data.z[,-1]                                  ## this matrix contains data for 48 regressors
no_z <- ncol(Z)

X <- Z #[,c(1,(no_z-1),no_z)]                     ## select 2 regressors for first estimation
no_x <- ncol(X)
invXX <- solve(crossprod(X))                  ## inverse to use later
invXXprior <- solve(crossprod(X)+diag(rep(1,no_x)))


## indicator functions for price, saddle, and end of saddle (1=it occurred, 0=it did not occur)
gamma_index <- data.all[,c(2,5,6)]




##########
## FIXED CONSTANTS
no_x <- ncol(X)       ## number of columns in X (regressors)
no_r <- 7

H <- 200              ## number of draws to compute integrals

no.treat <- 4         ## <- number of treatment effects
k.months <- 2         ## length of treatment effect intervals

N <- nrow(data.z)     ## number of products


k_theta <- rep(0,4)

treat1 <- data.all$t_saddle - data.all$t_price
k_theta[1] <- max(treat1[treat1>0])

treat2 <- (data.all$d_saddle+data.all$t_saddle)-data.all$t_price
k_theta[2] <- max(treat2[treat2>0])

treat3 <- data.all$t_price - data.all$t_saddle
k_theta[3] <- max(treat3[treat3>0])

treat4 <- data.all$t_price - (data.all$d_saddle+data.all$t_saddle)
k_theta[4] <- max(treat4[treat4>0])

Kt <- round(min(k_theta/k.months))
Kt <- 7

##########
#### INITIAL VALUES
kappa <- matrix(nrow=3,ncol=no_x,runif(no_x*3))           ## unform(0,1) random numbers for initial values of kappa coefficients)

no_v <- 3
P <- matrix(ncol=no_v,nrow=no_v,rnorm(no_v*no_v))                   ## initial values for covariance matrix of random effects
Om <- cov2cor(t(P) %*% P)
S <- diag(rep(0.4,no_v))
Omega <- S %*% Om %*% S

fix.draws <- halton.s(H,no_v)$draws                          ## fixed random normal draws

vs <- exp(t(chol(Omega)) %*% fix.draws)                   ## initial values for random effects

alpha <- c(0.08,-0.17,-0.06)                              ## initial values for alpha_s, alpha_p, alpha_d

theta <- matrix(ncol=Kt,nrow=no.treat,0) ## matrix with initial values for treatment effects
rownames(theta) <- c("theta_p_on_s","theta_p_on_d","theta_s_on_p","theta_d_on_p")

theta[1:2,] <- 0
theta[2,] <- 0
theta[3,] <- 0
theta[4,] <- 0

phi_x <- matrix(nrow=N,ncol=no_r,rnorm(N*no_r,mean=0,sd=0.5))

#B_tilde <- matrix(rnorm(no_x*no_r,sd=0.1),ncol=no_r,no_x)
B_tilde <- matrix(0,ncol=no_r,no_x)



##########
#### FUNCTIONS

## regressor function (works both as vector or as individual function)
f.phi<- function (x,kappa)
{
  temp1 <- exp(x %*% kappa)
  return(temp1)
}

## regressor function in treatment effect (likelihood 2)
f.phi_ds<- function (t,kappa)
{
  temp1 <- exp(t %*% kappa)
  return(temp1)
}

## hazard function
f.hazard <- function (t,alpha)
{
  temp1 <- exp(alpha * t)
  return(temp1)
}

## int. hazard function
f.int_hazard <- function (t,alpha)
{
  if(t<0){return(0)}
  temp1 <- 1/alpha * (exp(alpha*t)-1)
  return(temp1)
}


## function gets the treatment effect coefficient
f.sel_treat <- function (t,t_treat,theta,Kt)
{
  #browser()
  # t = time of event
  # t_treat = time of treatment
  # theta = treatment coefficients
  # Kt = number of intervals in treatment

  t_dif <- t-t_treat
  if(t_dif<0){return(0)}

  ## check positive t_dif
  ## if(t_dif<0){return("negative, treatment must be < t ")}
  if(t_dif<=0){return(0)}

  limit_left <- 0:(Kt-1)*k.months
  limit_right <-1:(Kt)*k.months

  if(t_dif>(Kt*k.months))
    {
     limit <- Kt
    }else{
          ind1 <- t_dif > limit_left & t_dif < limit_right
          limit <- which(ind1==TRUE)
    }

  sum1 <- exp(theta[limit])
  return(sum1)
}

f.int_treat <- function (t,t_treat,alpha,theta,Kt,phi_x_t)
{
  #browser()
  t <- as.numeric(t)
  t_treat <- as.numeric(t_treat)
  #browser()
  if(t<t_treat | t<0 | t_treat<0){return(0)} ## if treatment occurs after t, then returns 0

  ##alpha = hazard coefficient, theta=vector with treatment effects
  limit_vec <- sort(c((0:Kt*k.months)+t_treat,t))
  limit_vec <- limit_vec[limit_vec<=t]
  nlimit <- length(limit_vec)

  #int_hazard_vec <- apply(limit_vec,1,function(t){f.int_hazard(t,alpha=alpha)})
  int_hazard_vec <- sapply(limit_vec,function(t){f.int_hazard(t,alpha)})

  int_left <- int_hazard_vec[1:(nlimit-1)]
  int_right <- int_hazard_vec[2:nlimit]

  int_vec <- int_right - int_left

    if(nlimit==(Kt+2))
    {
     int_t1 <- int_vec[1:Kt] * exp(theta[1:Kt])
     int_t2 <- int_vec[Kt+1] * exp(theta[Kt])
     int_treat <- c(int_t1,int_t2)
    }else{
    int_treat <- int_vec[1:(nlimit-1)] * exp(theta[1:(nlimit-1)])
    }

  int_c <- sum(int_treat)*phi_x_t #*vs
  int_d <- mean(int_c)
  return(int_d)

}


f.int_treat_two <- function (t,t_treat_a,t_treat_b,alpha,theta_a,theta_b,Kt,phi_x_t)
{
  #browser()
  t <- as.numeric(t)
  t_treat_a <- as.numeric(t_treat_a)
  t_treat_b <- as.numeric(t_treat_b)
  #browser()
  if(t<t_treat_b | t<0 | t_treat_b<0){return(0)} ## if treatment occurs after t, then returns 0

  ##alpha = hazard coefficient, theta=vector with treatment effects
  int_a <- c((0:Kt*k.months))+t_treat_a ; Ka <- length(int_a)
  int_b <- c((0:Kt*k.months))+t_treat_b ; Kb <- length(int_b)

  int_a_left <- int_a[1:(Ka-1)]
  int_a_right <- int_a[2:Ka]
  all_a <- cbind(int_a_left,int_a_right,theta_a)

  int_b_left <- int_b[1:(Kb-1)]
  int_b_right <- int_b[2:Kb]
  all_b <- cbind(int_b_left,int_b_right,theta_b)
             
  limit_vec <- sort(c(int_a,int_b,t))
  limit_vec <- limit_vec[limit_vec<=t & limit_vec >= t_treat_b]             
  nlimit <- length(limit_vec)

  limits_col <- cbind(limit_vec[1:(nlimit-1)],limit_vec[2:nlimit],rep(0,(nlimit-1)),rep(0,(nlimit-1)))

  limits_col;
  
  for (i in 1:(nlimit-1))
  {
      if(limits_col[i,2]<=max(all_a[,2]))
        {
        limits_col[i,3] <- all_a[,3][limits_col[i,1] >= all_a[,1] & limits_col[i,2] <= all_a[,2]]
        }else{
        limits_col[i,3] <- all_a[Ka-1,3]
        }
    
      if(limits_col[i,2]<=max(all_b[,2]))
        {
        limits_col[i,4] <- all_b[,3][limits_col[i,1] >= all_b[,1] & limits_col[i,2] <= all_b[,2]]
        }else{
        limits_col[i,4] <- all_b[Kb-1,3]
        }
  }             
             
  ## hazard integral
  int_hazard_vec <- sapply(limit_vec,function(t){f.int_hazard(t,alpha)})
             
  int_left <- int_hazard_vec[1:(nlimit-1)]
  int_right <- int_hazard_vec[2:nlimit]
             
  int_vec <- int_right - int_left

  if(nlimit==(Kt+2))
  {
   int_t1 <- int_vec[1:Kt] * exp(limits_col[1:Kt,3]) * exp(limits_col[1:Kt,4])
   int_t2 <- int_vec[Kt+1] * exp(limits_col[Kt,3]) * exp(limits_col[Kt,4]) 
   int_treat <- c(int_t1,int_t2)
  }else{
    int_treat <- int_vec * exp(limits_col[,3]) * exp(limits_col[,4])
  }
             
  int_c <- sum(int_treat)*phi_x_t #*vs
  int_d <- mean(int_c)
  return(int_d)
}




f.group.like <- vector("list",3)

## likelihood for t_s
f.group.like[[1]] <- f.like.ts <- function (t_tilde,gamma_index,theta,phi_x,alpha,vs)
{
   #browser()
   t_tilde <- as.numeric(t_tilde)
   t_s <- t_tilde[2]       # saddle
   t_p <- t_tilde[1]       # price
   theta <- theta[1,]                  # price on saddle

   gamma_index <- as.numeric(gamma_index)
   gamma_p <- gamma_index[1]
   gamma_s <- gamma_index[2]
   gamma_r <- gamma_index[3]

   phi_xm <- exp(phi_x[1]) #f.phi(x,kappa)
   phi_xt <- exp(phi_x[4])

   vs_ts <- vs[1,]
   #vs_tpts <- vs[4,]

   lambda_s <- f.hazard(t_s,alpha)

   line1_index <- gamma_s * as.numeric (t_s < t_p)
   if (line1_index==1)
   {
   line1_a <- line1_index * phi_xm * lambda_s * vs_ts
   line1_b <- - phi_xm * f.int_hazard(t_s,alpha) * vs_ts
   line1 <-  line1_a * exp(line1_b)
   lineA <- line1
   }

   line2_index <- gamma_p * as.numeric(t_s>t_p)
   if (line2_index==1)
   {
   line2_a <- {phi_xm * lambda_s * f.sel_treat(t_s,t_p,theta,Kt) * vs_ts }^gamma_s
   line2_b <- -phi_xm * { {f.int_hazard(t_p,alpha) + f.int_treat(t_s,t_p,alpha,theta,Kt,phi_xt)} * vs_ts }
   line2 <- line2_index * line2_a * exp(line2_b)
   lineA <- line2
   }

   line3_index <- (1-gamma_s) * (1-gamma_p)
   if (line3_index==1)
   {
   line3 <- line3_index * exp(- phi_xm * f.int_hazard(t_s,alpha) * vs_ts )
   lineA <- line3
   }

   #like.int <- as.numeric(lineA)
   like.int <- lineA
   #if(like.int==0){like.int <- exp(-745)}
   return(like.int)
}

## likelihood for t_p
f.group.like[[2]] <- f.like.tp <- function (t_tilde,gamma_index,theta,phi_x,alpha,vs)
{
  #browser()
  gamma_index <- as.numeric(gamma_index)
  gamma_p <- gamma_index[1]
  gamma_s <- gamma_index[2]
  gamma_r <- gamma_index[3]

  t_tilde <- as.numeric(t_tilde)
  t_p   <-  t_tilde[1] ## price
  t_s   <-  t_tilde[2] ## saddle
  t_end <-  t_tilde[2]+t_tilde[3] ## saddle end

  theta_ts <- theta[3,] # ts on tp
  theta_td <- theta[4,] # td on tp

  phi_xm <- exp(phi_x[2]) #f.phi(x,kappa)
  vs_tp <- vs[2,]
  #vs_ts <- vs[5,]
  #vs_te <- vs[6,]

  lambda_p <- f.hazard(t_p,alpha)
  Lambda_p <- f.int_hazard(t_p,alpha)
  Lambda_s <- f.int_hazard(t_s,alpha)

  phi_x_ts <- exp(phi_x[5])
  phi_x_tend <- exp(phi_x[6])

  line1_index <- gamma_p * as.numeric(t_p < t_s)
  if (line1_index==1)
  {
  line1_a <- phi_xm * lambda_p * vs_tp
  line1_b <- -phi_xm * Lambda_p * vs_tp
  line1 <- line1_index * line1_a * exp(line1_b)
  lineA <- line1
  }

  line2_index <- gamma_s * as.numeric ( t_s < t_p & t_p <= t_end )
  if (line2_index==1)
  {
  line2_a <- {phi_xm * lambda_p * f.sel_treat(t_p,t_s,theta_ts,Kt) * vs_tp  }^gamma_p
  line2_b <- -phi_xm * {Lambda_p + f.int_treat(t_p,t_s,alpha,theta_ts,Kt,phi_x_ts)} * vs_tp
  line2 <- line2_index * line2_a * exp(line2_b)
  lineA <- line2
  }

  line3_index <- gamma_s  * gamma_r * as.numeric(t_p > t_end)
  if (line3_index==1)
  {
  line3_a <- {phi_xm * lambda_p * f.sel_treat(t_p,t_s,theta_ts,Kt) * f.sel_treat(t_p,t_end,theta_td,Kt) * vs_tp }^gamma_p
  line3_b <- -phi_xm * { Lambda_s + f.int_treat(t_p,t_s,alpha,theta_ts,Kt,phi_x_ts) + f.int_treat_two(t_p,t_s,t_end,alpha,theta_ts,theta_td,Kt,phi_x_tend) } * vs_tp
  line3 <- line3_index * line3_a * exp(line3_b)
  lineA <- line3
  }

  line4_index <- (1-gamma_p) * (1-gamma_s)
  if (line4_index==1)
  {
  line4 <- line4_index * exp(-phi_xm * Lambda_p * vs_tp )
  lineA <- line4
  }

  #like.int <-  as.numeric(lineA)
  like.int <- lineA
  return(like.int)
}


f.group.like[[3]] <- f.like.td <- function (t_tilde,gamma_index,theta,phi_x,alpha,vs)
{
  #browser()
  gamma_index <- as.numeric(gamma_index)
  gamma_p <- gamma_index[1]
  gamma_s <- gamma_index[2]
  gamma_d <- gamma_index[3]

  t_tilde <- as.numeric(t_tilde)
  t_d <- t_tilde[3]  # duration
  t_s <- t_tilde[2]
  t_e <- t_d + t_s
  t_p <- t_tilde[1]  #price

  theta_p <- theta[2,] ## treatment p on d

  phi_xm    <- exp(phi_x[3]) #f.phi(x,kappa)
  phi_xt <- exp(phi_x[7])
  vs_td <- vs[3,]
  #vs_tp <- vs[7,]

  lambda_d <- f.hazard(t_d,alpha)
  Lambda_d <- f.int_hazard(t_d,alpha)

  h_d <- 1 ## pending to add h() function here

  line0 <- 1 - gamma_s
  if (line0==1)
  {
  lineA <- rep(line0,200)
  }


  line1_index <- gamma_s * gamma_p * as.numeric(t_p < t_s)
  if (line1_index==1)
  {
  line1_a <- exp(-phi_xm * Lambda_d * h_d * vs_td)
  line1_b <-  {phi_xm * lambda_d * h_d * vs_td }^gamma_d
  line1 <- line1_index * line1_a * line1_b
  lineA <- line1
  }


  line2_index <- gamma_s * gamma_p * as.numeric(t_p > t_s & t_p <= t_e )
  if (line2_index==1)
  {
  line2_a <- - phi_xm * h_d * { f.int_hazard(t_p-t_s,alpha) + f.int_treat(t_d,(t_p-t_s),alpha,theta_p,Kt,phi_xt) } * vs_td
  line2_b <- {phi_xm * h_d * lambda_d * f.sel_treat(t_d+t_s,t_p,theta_p,Kt) * vs_td }^gamma_d
  line2 <- line2_index * exp(line2_a) * line2_b
  lineA <- line2
  }

  line3_index <- gamma_s * gamma_p * as.numeric(t_p>t_e)
  if (line3_index==1)
  {
  line3 <- line3_index * exp(-phi_xm * Lambda_d * h_d * vs_td )
  lineA <- line3
  }


  line4_index <- gamma_s * (1-gamma_p)
  if (line4_index==1)
  {
  line4_a <- exp(-phi_xm * Lambda_d * h_d * vs_td )
  line4_b <- {phi_xm * lambda_d * h_d * vs_td }^gamma_d
  line4 <- line4_index * line4_a * line4_b
  lineA <- line4
  }

  like.int <- lineA
  return(like.int)
}


f.like.all<- function (t_tilde,gamma_index,theta,phi_x,alpha,vs)
{
  #browser()
  like_ts <- f.like.ts(t_tilde,gamma_index,theta,phi_x,alpha[1],vs)
  like_tp <- f.like.tp(t_tilde,gamma_index,theta,phi_x,alpha[2],vs)
  like_td <- f.like.td(t_tilde,gamma_index,theta,phi_x,alpha[3],vs)

  like.p <- like_ts * like_tp * like_td
  like.int <- sum(like.p) / length(like.p)
  return(like.int)
}


##########
## PRIORS FUNCTIONS

f.log.prior.phi<- function (phi_r,Omega)
{
  vec1 <- t(chol(solve(Omega))) %*% phi_r
  vecp <- sum(log(dnorm(vec1)))
  return(vecp)
}

f.log.prior.theta<- function (theta,sigma2theta,sigma2one)
{
  nt <- ncol(theta)
  dif <- theta[2:nt] - theta[1:(nt-1)]
  vecp1 <- sum(log(dnorm(theta[1],mean=0,sd=sqrt(sigma2one))))
  vecp0 <- sum(log(dnorm(dif,sd=sqrt(sigma2theta))))
  return(vecp0 + vecp1)
}

f.log.post.phi <- function (t_tilde,gamma_index,theta,phi_x,alpha,phi_bar,Omega)
{
  ## posterior of phi values
  like <- log(f.like.all(t_tilde,gamma_index,theta,phi_x,alpha))
  phi_r <- phi_x - phi_bar
  prior <- f.log.prior.phi(phi_r, Omega)
  post <- like+prior
  return(post)
}

##########
## METROPOLIS HASTINGS FUNCTIONS

counter_met_phi <- rep(0,2)
f.metrop.phi <- function (t_tilde,gamma_index,theta,phi_x,alpha,phi_bar,Omega)
{
  #browser()
  phi_current <- phi_x
  phi_draw <-  phi_x + rnorm(7,sd=0.2)

  den <- post_current <- f.log.post.phi(t_tilde,gamma_index,theta,phi_current,alpha,phi_bar,Omega)
  num <- post_next <- f.log.post.phi(t_tilde,gamma_index,theta,phi_draw,alpha,phi_bar,Omega)

  if (runif(1)<=min(exp(num-den),1))
  {
      phi_next <- phi_draw
      counter_met_phi[1] <<- counter_met_phi[1] + 1
      #print("accept")
  }
  else
  {
  phi_next <- phi_current
  counter_met_phi[2] <<- counter_met_phi[2] + 1
  }

  return(phi_next)
}

counter_alpha <- rep(0,2)
f.metrop.alpha<- function (t_tilde,gamma_index,theta,phi_x,alpha,k,vs)
{
  #browser()
  alpha_current <- alpha[k]
  alpha_draw <- alpha[k] + rnorm(1,sd=0.1/15)
  prior.alpha <- 1

  like_ts_a <- sapply(1:N,function(i)f.group.like[[k]](t_tilde[i,],gamma_index[i,],theta,phi_x[i,],alpha_current,vs))
  like_ts_b <- apply(like_ts_a,2,mean); like_ts_b[like_ts_b==0] <- exp(-745)
  like_ts_current <- sum(log(like_ts_b))

  den <- post_log_current <- like_ts_current +log(prior.alpha)

  like_ts_a <- sapply(1:N,function(i)f.group.like[[k]](t_tilde[i,],gamma_index[i,],theta,phi_x[i,],alpha_draw,vs))
  like_ts_b <- apply(like_ts_a,2,mean); like_ts_b[like_ts_b==0] <- exp(-745)
  like_ts_draw <- sum(log(like_ts_b))

  num <- post_log_next <- like_ts_draw + log(prior.alpha)

  if (runif(1)<=min(exp(num-den),1))
  {
   alpha_next <- alpha_draw
   counter_alpha[1] <<- counter_alpha[1] + 1
  }else{
  alpha_next <- alpha_current
  counter_alpha[2] <<- counter_alpha[2] + 1
  }

  return(alpha_next)
}


counter_theta <- rep(0,2)
f.metrop.theta.vec <- function (t_tilde,gamma_index,theta,phi_x,alpha,k,vs)
{
  #browser()
  if(k==1)f <- 1
  if(k==2)f <- 3
  if(k==3)f <- 2
  if(k==4)f <- 2

  theta_current <- theta
  theta_draw <- theta;
  #theta_draw[k,] <- theta_draw[k,] + t(chol(cov_theta[[k]]*1/10)) %*% rnorm(Kt)
  theta_draw[k,] <- theta_draw[k,] + rnorm(Kt,sd=0.1)

  #theta_draw[k,2:Kt] <- theta[k,2:Kt] + rnorm((Kt-1),sd=0.1)
  #theta_draw[k,1] <- theta[k,1] + rnorm(1,sd=0.1/10)


  like_ts_a <- sapply(1:N,function(i)f.group.like[[f]](t_tilde[i,],gamma_index[i,],theta_current,phi_x[i,],alpha[f],vs))
  like_ts_b <- apply(like_ts_a,2,mean); like_ts_b[like_ts_b==0] <- exp(-745)
  like_ts_current <- sum(log(like_ts_b))
  log_prior_theta_c <- f.log.prior.theta(theta_current,1,1/10)
  den <- like_ts_current + log_prior_theta_c

  like_ts_a <- sapply(1:N,function(i)f.group.like[[f]](t_tilde[i,],gamma_index[i,],theta_draw,phi_x[i,],alpha[f],vs))
  like_ts_b <- apply(like_ts_a,2,mean); like_ts_b[like_ts_b==0] <- exp(-745)
  like_ts_draw <- sum(log(like_ts_b))
  log_prior_theta_d <- f.log.prior.theta(theta_draw,1,1/10)
  num <- like_ts_draw + log_prior_theta_d

  if (runif(1)<=min(exp(num-den),1))
  {
   theta_next <- theta_draw
   counter_theta[1] <<- counter_theta[1] + 1
  }else{
  theta_next <- theta_current
  counter_theta[2] <<- counter_theta[2] + 1
  }

  return(theta_next)
}


counter_beta <- rep(0,2)
f.metrop.beta<- function (t_tilde,gamma_index,theta,b_tilde,alpha,k,vs,X,sigma2met,sigmaBprior)
{
  #browser()
  if(k==1)f <- 1
  if(k==4)f <- 1
  if(k==2)f <- 2
  if(k==5)f <- 2
  if(k==6)f <- 2
  if(k==3)f <- 3
  if(k==7)f <- 3

  beta_draw <- beta_current <- b_tilde
  beta_draw[,k] <- beta_current[,k] + rnorm(no_x,sd=sigma2met)

  phi_x_current <-  X %*% beta_current
  phi_x_draw <- X %*% beta_draw

  like_ts_a <- sapply(1:N,function(i)f.group.like[[f]](t_tilde[i,],gamma_index[i,],theta,phi_x_current[i,],alpha[f],vs))
  like_ts_b <- apply(like_ts_a,2,mean); like_ts_b[like_ts_b==0] <- exp(-745)
  like_ts_current <- sum(log(like_ts_b))
  log_prior_beta_c <- sum(log(dnorm(beta_current[,k],sd=sigmaBprior)))
  den <- post_log_current <- like_ts_current + log_prior_beta_c

  like_ts_a <- sapply(1:N,function(i)f.group.like[[f]](t_tilde[i,],gamma_index[i,],theta,phi_x_draw[i,],alpha[f],vs))
  like_ts_b <- apply(like_ts_a,2,mean); like_ts_b[like_ts_b==0] <- exp(-745)
  like_ts_draw <- sum(log(like_ts_b))
  log_prior_beta_c <- sum(log(dnorm(beta_draw[,k],sd=sigmaBprior)))
  num <- post_log_next <- like_ts_draw +  log_prior_beta_c

  if (runif(1)<=min(exp(num-den),1))
  {
   beta_next <- beta_draw
   counter_beta[1] <<- counter_beta[1] + 1
  }else{
  beta_next <- beta_current
  counter_beta[2] <<- counter_beta[2] + 1
  }

  return(beta_next)
}


counter_omega <- rep(0,2)
f.metrop.omega <- function (t_tilde,gamma_index,theta,phi_x,alpha,vs,P,S)
{
   #browser()
   vs_current <- vs
   P_current <- P
   S_current <- S

   s_diag <- c(rnorm(no_v,sd=0.1/7)) #,rnorm(4,sd=0.1/no_v))
   P_draw <- P + matrix(nrow=no_v,ncol=no_v,rnorm(no_v*no_v,sd=1/7))
   Om <- cov2cor(t(P_draw) %*% P_draw)
   S_draw <- S + diag(s_diag)
   Omega <- S_draw %*% Om %*% S_draw

   vs_draw <- exp(t(chol(Omega)) %*% fix.draws)

   like_a <- sapply(1:N,function(i){f.like.all(t_tilde[i,],gamma_index[i,],theta,phi_x[i,],alpha,vs_current)})
   like_a <- like_a + min(like_a[like_a!=0])
   like_b <- sum(log(like_a))
   log_prior_c <- sum(log(dnorm(S_current,sd=1/20)))
   den <- log_post_current <- like_b + log_prior_c

   like_a <- sapply(1:N,function(i){f.like.all(t_tilde[i,],gamma_index[i,],theta,phi_x[i,],alpha,vs_draw)})
   like_a <- like_a + min(like_a[like_a!=0])
   like_b <- sum(log(like_a))
   log_prior_c <- sum(log(dnorm(S_draw,sd=1/20)))
   num <- log_post_draw <- like_b + log_prior_c

   if (runif(1)<=min(exp(num-den),1))
   {
   p_next <- P_draw
   vs_next <- vs_draw
   s_next <- S_draw
   counter_omega[1] <<- counter_omega[1] + 1
   }else{
   p_next <- P_current
   vs_next <- vs_current
   s_next <- S_current
   counter_omega[2] <<- counter_omega[2] + 1
   }

   return(list(P=p_next,vs=vs_next,S=s_next))
}




##########
## SIMULATION SET UP
nosim <- 6000
nothin <- 1
nosave <- nosim/nothin

noinits <- 100
load_inits <- 1
save_inits <- 1


##########
## STORAGE OF DRAWS
phi_x_store <- array(dim=c(N,7,nosave))
beta_store <- array(dim=c(no_x,7,nosave))
omega_store <- array(dim=c(no_v,no_v,nosave))
alpha_store <- array(dim=c(3,1,nosave))
theta_store <- array(dim=c(4,Kt,nosave))

phi_x_hat <- matrix(nrow=N,ncol=7,0)


##########
### ESTIMATION ALGORITHM
counter_alpha <- rep(0,2)
counter_met_phi <- rep(0,2)
counter_theta <- rep(0,2)

report <- cbind(0,0,0,0,0,0)
colnames(report) <- c("Iter||","Time per Iteration||","Acc. Rate Alpha||", "Acc. Rate Beta||", "Acc. Rate Theta||","Acc. Rate Omega||")


if(load_inits==1)
{
seed1 <- seed1_init
set.seed(seed1)
alpha <- alpha_init
phi_x <- phi_x_init
theta <- theta_init
theta[4,] <- 0
#Omega <- omega_init
B_tilde <- beta_init
#P <- p_init
#vs <- vs_init
#S <- diag(sqrt(diag(omega_init)))

}


time1 <- Sys.time()
for (j in 1:nosim)
{


 ##DRAW OF ALPHA
 for (k in 1:3)
 {
 #print(k)
 alpha[k] <- f.metrop.alpha(t_tilde,gamma_index,theta,phi_x,alpha,k,vs)
 }



 ## DRAW OF THETA
 for (k in 1:4)
 {
 theta <- f.metrop.theta.vec(t_tilde,gamma_index,theta,phi_x,alpha,k,vs)
 }

 ## DRAW OF BETA (#1/30)
 for (k in 1:3)
 {
     B_tilde <- f.metrop.beta(t_tilde,gamma_index,theta,B_tilde,alpha,k,vs,X,1/50,0.1)
 }

 ## DRAW OF PHI
 phi_x <-  X %*% B_tilde


 ## DRAW OF OMEGA
 domega <- f.metrop.omega(t_tilde,gamma_index,theta,phi_x,alpha,vs,P,S)
 P <- domega$P
 vs <- domega$vs
 S <- domega$S
 Omega <- S %*% cov2cor(t(P) %*% P) %*% S


 print(alpha)
 print(Omega)
 print(cov2cor(Omega))
 print(theta)

 ## SAVE DRAWS
 if(j%%nothin==0)
   {
        j_s <- j/nothin
        phi_x_store[,,j_s]  <- phi_x
        beta_store[,,j_s] <- B_tilde
        omega_store[,,j_s] <- Omega
        alpha_store[,,j_s] <- alpha
        theta_store[,,j_s] <- theta
        par(mfcol=c(4,3))
        par(mar=c(2,4,3,3))
        plot(alpha_store[1,,1:j_s],type="l")
        plot(alpha_store[2,,1:j_s],type="l")
        plot(omega_store[1,1,1:j_s],type="l")
        plot(theta_store[1,1,1:j_s],type="l")
        plot(rep(theta[1,],each=2),type="l")
        plot(rep(theta[2,],each=2),type="l")
        plot(rep(theta[3,],each=2),type="l")
        plot(rep(theta[4,],each=2),type="l")
        plot(beta_store[1,1,1:j_s],type="l")
        plot(beta_store[1,2,1:j_s],type="l")
        plot(beta_store[1,3,1:j_s],type="l")
        plot(B_tilde[,1],type="l")
   }


 if(j%%noinits==0)
 {
 if(save_inits==1)
 {
 alpha_init <- alpha
 phi_x_init <- phi_x
 theta_init <- theta
 omega_init <- Omega
 beta_init <- B_tilde
 p_init <- P
 vs_init <- vs
 seed1_init <- seed1
 fix.draws_init <- fix.draws
 save(list=ls(pattern="_init"),file="inits.Rdata")
 }
 }



 ##TIME REPORT
 report[1,1] <- j
 report[1,2] <- round(Sys.time() - time1,2) ; time1 <- Sys.time()
 report[1,3] <- counter_alpha[1]/sum(counter_alpha)
 report[1,4] <- counter_beta[1]/sum(counter_beta)
 report[1,5] <- counter_theta[1]/sum(counter_theta)
 report[1,6] <- counter_omega[1]/sum(counter_omega)
 print(report)

}


save(list=ls(),file="output.Rdata")











## stop()
## ##########
## ## Short Test

## ## calculates time needed for computing likelihood

## phi_x <- cbind(phi_x_ts,phi_x_tp,phi_x_td)

## system.time(
## like.all <- sapply(1:N,function(i){f.like.all(t_tilde[i,],X[i,],gamma_index[i,],theta,kappa,alpha,vs)})
##             )


## ## profile likelihood
## ev <- seq(-0.5,1,by=0.03)
## J <- length(ev)
## like.eval <- rep(0,J)
## for (i in 1:J)
## {
##     alpha[1] <- ev[i]
##     like.all <- sapply(1:N,function(i){f.like.all(t_tilde[i,],X[i,],gamma_index[i,],theta,kappa,alpha,vs)})
##     like.all <- like.all + min(like.all[like.all!=0])
##     like.eval[i] <- sum(log(like.all))
## }


## ## profiles the functions
## Rprof()
## like.all <- sapply(1:N,function(i){f.like.all(t_tilde[i,],X[i,],gamma_index[i,],theta,kappa,alpha,vs)})
## Rprof(NULL)
## summaryRprof()

g1 <- rep(0,N)
for (i in 1:N)
{
    g1a <- f.like.td(t_tilde[i,],gamma_index[i,],theta,phi_x[i,],alpha[1],vs)
    #g1a <- f.like.all(t_tilde[i,],gamma_index[i,],theta,phi_x[i,],alpha,vs)
    g1[i] <- sum(g1a)/200
}

i <- 1
like_ts_current<- sapply(1:N,function(i)f.group.like[[3]](t_tilde[i,],gamma_index[i,],theta,phi_x[i,],alpha[3],vs))

sum(log(apply(like_ts_current,2,sum)/200))



