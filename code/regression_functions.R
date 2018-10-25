controller <- function(X,y,B) {
  set.seed(414)
  r <- regress(X, y, B)
  
  # Burn-in done by the regress function
  # Convergence
  rm <- apply(r$coeff,2,function(x) (cumsum(x)/(1:length(x)))-mean(x))
  
  # Summary
  s <- summary_coeff(r$coeff)
  
  # Capture densities
  all_d <- apply(r$coeff,2,density)
  
  return(list(regress=r, converge=rm, summary=s, densities=all_d))
}

norm <- function(x) {
  return(sum(x^2))
}

rinvchisq <- function(n,v,s2) {
  return(rinvgamma(n,v/2,v*s2/2))
}

# Extract data components, center each column
regress <- function(X,y,B) {
  require(mvtnorm)
  require(MCMCpack)
  
  # Determine parameters
  sizes <- dim(X)
  n <- sizes[1]
  k <- sizes[2]+1
  nk <- n-k
  
  # Initialize storage
  beta_mat <- matrix(NA, nrow = B, ncol = k)
  sigma2_vec <- rep(NA,B)
  
  # Prepare data structures
  X <- apply(X, 2, function(x) x - mean(x)) # Center data
  X <- cbind(intercept = rep(1,n),X) # Add intercept column
  colnames(beta_mat) <- colnames(X)
  
  XTy <- t(X)%*%y
  V_beta <- solve(t(X)%*%X)
  beta_hat <- V_beta%*%XTy
  s2 <- norm(y-X%*%beta_hat)/nk
  
  # Sample from Joint
  for (i in 1:B) {
    sigma2_vec[i] <- rinvchisq(1,nk,s2)
    beta_mat[i,] <- rmvnorm(1,beta_hat,V_beta*sigma2_vec[i])
  }
  return(list(dmat=X,coeff=beta_mat[((B/2)+1):B,],
              variance=sigma2_vec[((B/2)+1):B]))
}

summary_coeff <- function(coeff_mat) {
  m1 <- apply(coeff_mat,2,mean)
  m2 <- apply(coeff_mat,2,median)
  ci <- apply(coeff_mat,2,quantile,probs=c(0.025,0.975))
  s <- rbind(mean=m1, median=m2, ci)
  return(s)
}

plot_dist <- function(title,credInt,d) {
  plot(d,col='blue',xlab='Value',main=title)
  idStart <- max(which(d$x < credInt[1])) + 1
  idEnd <- min(which(d$x > credInt[2])) - 1
  gx <- d$x[idStart:idEnd]
  gy <- d$y[idStart:idEnd]
  px <- rep(0, length(gy))
  polygon(c(gx, rev(gx)), c(px, rev(gy)), border = FALSE,
          col = rgb(0, 0, 1, alpha = 0.5))
}

dic_norm <- function(y,results) {
  # Extract needed variables
  all_beta <- results$regress$coeff
  beta_mean <- results$summary[1,]
  dmat <- results$regress$dmat
  sig2 <- results$regress$variance
  y <- matrix(y,nrow=1)
  
  # Computer mean parameters
  I <- diag(length(y))
  mu_mean <- dmat%*%matrix(beta_mean,ncol=1)
  sigma_mean <- sqrt(mean(sig2))*I
  
  lik_mean <- dmvnorm(y,mean=mu_mean,sigma=sigma_mean,log=TRUE)
  
  # Initialize storage
  B <- length(sig2)
  ylik <- rep(NA,B)
  
  # Compute likelihood at each sample
  for (i in 1:B) {
    b <- matrix(all_beta[i,],ncol=1)
    s <- sig2[i]
    ylik[i] <- dmvnorm(y,mean=dmat%*%b,sigma=sqrt(s)*I,log=TRUE)
  }
  pdic <- 2*(lik_mean-(sum(ylik)/B))
  DIC <- -2*lik_mean+2*pdic
  return(DIC)
}