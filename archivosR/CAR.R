CAR<-function(Y,A,iters=25000,burn=5000){

   tick   <- proc.time()[3]

   # Bookkeeping

    n     <- length(Y)  
    m     <- rowSums(A)
    adj   <- apply(A==1,1,which)
    nei1  <- row(A)[A==1]
    nei2  <- col(A)[A==1]

   # Initial values

    sig2e <- var(Y)/2
    sig2s <- var(Y)/2
    rho   <- 0.99
    mu    <- mean(theta)
    theta <- 0.5*(Y-mu)

   # Pre-compute the determinant for all rho of interest

    C        <- diag(1/sqrt(m))%*%A%*%diag(1/sqrt(m))
    lambda   <- eigen(C)$values
    rho.grid <- seq(0.5,0.999,0.001)
    logd     <- rho.grid
    for(j in 1:length(logd)){
      logd[j] <- sum(log(1-rho.grid[j]*lambda))
    }   


   # Keep track of stuff

    keepers <- matrix(0,iters,4)
    colnames(keepers) <- c("sig2e","sig2s","mu","rho")
    theta1 <- theta2 <- 0

   # GO!!!

   for(iter in 1:iters){

      # THETA

       for(j in 1:n){
         neig     <- adj[[j]]
         ybar     <- mean(theta[neig])
         VVV      <- m[j]/sig2s + 1/sig2e
         MMM      <- rho*ybar*m[j]/sig2s + (Y[j]-mu)/sig2e
         theta[j] <- rnorm(1,MMM/VVV,1/sqrt(VVV))
       }

      # VARIANCES

       TAT   <- sum(theta[nei1]*theta[nei2])
       TMT   <- sum(m*theta^2)
       sig2s <- 1/rgamma(1,n/2+.1,(TMT-rho*TAT)/2+.1)
       sig2e <- 1/rgamma(1,n/2+.1,sum((Y-mu-theta)^2)+.1)

      # MEAN

       VVV  <- n/sig2e + 0.001
       MMM  <- sum(Y-theta)/sig2e
       mu   <- rnorm(1,MMM/VVV,1/sqrt(VVV))

      # CAR DEPENDENCE PARAMETER

       R    <- 0.5*logd + 0.5*rho.grid*TAT/sig2s
       rho  <- sample(rho.grid,1,prob=exp(R-max(R)))

      # KEEP TRACK OF STUFF

       keepers[iter,] <- c(sig2e,sig2s,mu,rho)
       if(iter>burn){
         theta1 <- theta1 + theta/(iters-burn)
         theta2 <- theta2 + theta*theta/(iters-burn)
       }

   }
   tock   <- proc.time()[3]
   output <- list(samps=keepers,
                  theta.mn=theta1,theta.var=theta2-theta1^2,
                  minutes=(tock-tick)/60)

 return(output)}