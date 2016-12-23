# Function for computing the likelihood used in optimization.
lll <- function(lambda, a=NULL, A=NULL, n=NULL, m=NULL, nlatent=NULL)
{
  lambda = matrix(lambda, m, nlatent);
  val = 0;
  for (i in 1:n)
  {
    for (j in 1:m)
    {
      # Calculate helper terms.
      eta = lambda[j,]%*%a[i,];
      psi = pnorm(eta);
      # For numerical stability.
      if (psi > 0.5) psi = psi - 1e-10
      else psi = psi + 1e-10
      # Add current terms of the likelihood.
      val = val + y[i,j]*log(psi) + (1-y[i,j])*log(1 - psi);
      val = val - 1/2 * t(lambda[j,]) %*% A[[i]] %*% lambda[j,];
    }
    # Add terms that do not depend on j.
    val = val + 1/2*(log(determinant(A[[i]], logarithm=FALSE)$modulus[1]) - sum(diag(A[[i]])) - sum(a[i,]^2));
  }
  return(val);
}

# Function for computing the likelihood used in optimization.
lla <- function(a, lambda=NULL, A=NULL, n=NULL, m=NULL, nlatent=NULL)
{
  a = matrix(a, n, nlatent);
  val = 0;
  for (i in 1:n)
  {
    for (j in 1:m)
    {
      # Calculate helper terms.
      eta = lambda[j,]%*%a[i,];
      psi = pnorm(eta);
      # For numerical stability.
      if (psi > 0.5) psi = psi - 1e-10
      else psi = psi + 1e-10
      # Add likelihood terms.
      val = val + y[i,j]*log(psi) + (1-y[i,j])*log(1 - psi);
      val = val - 1/2 * t(lambda[j,]) %*% A[[i]] %*% lambda[j,];
    }
    # Add terms that do not depend on j.
    val = val + 1/2*(log(determinant(A[[i]], logarithm=FALSE)$modulus[1]) - sum(diag(A[[i]])) - sum(a[i,]^2));
  }
  return(val);
}

# Calculate the gradient with respect to the variational means, a.
grad.a <- function(a, lambda=NULL, A=NULL, n=NULL, m=NULL, nlatent=NULL)
{
  a = matrix(a, n, nlatent);
  g.a = NULL; 
  
  for (i in 1:n) 
  {
    grad_a = rep(0, nlatent);
    for (j in 1:m)
    {
      # Calculate helper terms.
      eta = lambda[j,]%*%a[i,];
      psi = pnorm(eta)+ 1e-10;
      # Calculate gradients.
      grad_a = grad_a + (y[i,j]- psi) *dnorm(eta)/(psi*(1-psi)+1e-10)*lambda[j,];
    }
    g.a <- c(g.a, grad_a - a[i,]);
  }
  # Make sure the matrix is formatted properly.
  return(matrix(g.a, n, nlatent, byrow=TRUE));
}

# Caclulate the gradient with respect to the model parameters lambda.
grad.lambda <- function(lambda, a=NULL, A=NULL, n=NULL, m=NULL, nlatent=NULL)
{
  lambda <- matrix(lambda, m, nlatent)
  g.l = NULL;
  for (j in 1:m)
  {
    grad_l = 0
    for (i in 1:n)
    {
      # Calculate helper terms.
      eta = a[i,]%*%lambda[j,];
      psi = pnorm(eta)+ 1e-10;
      grad_l = grad_l + (y[i,j]-psi)*dnorm(eta)/(psi*(1-psi)+1e-10)*a[i,] - 0.5 * (A[[i]] + t(A[[i]]))%*%lambda[j,]
    }
    
    # Update parameters.
    g.l <- c(g.l, grad_l);
  }
  # Make sure the gradient is formatted properly.
  return(matrix(g.l, m, nlatent, byrow=TRUE));
}

gllvm.fit <- function(ys, nlatent=2, maxiters=1000)
{
  n <- dim(ys)[c(1)];
  m <- dim(ys)[c(2)];
  
  # Initialize parameters.
  lambda <- matrix(rnorm(m*nlatent), nrow=m, ncol=2);
  a <- matrix(rnorm(n*nlatent), nrow=n, ncol=2);
  A <- list(); for(i in 1:n) { A[[i]] <- diag(nlatent);}
  
  n_it = 0
  while(TRUE)
  {
    n_it = n_it + 1;
    
    # Stage (1): Update a_i's.
    q <- optim(a, method="BFGS", fn=lla, gr=grad.a, lambda=lambda, A=A, n=n, m=m, nlatent=nlatent, control=list(trace=0, fnscale=-1, maxit = 100, reltol=1e-2))
    new.a <- matrix(q$par[1:(nlatent*n)], n, nlatent, byrow=FALSE);
    
    # Stage (2): Update lambda_j's.
    lower.vec <- rep(-30,length(c(lambda))); 
    upper.vec <- rep(30,length(c(lambda)))
    lower.vec[which(upper.tri(lambda))] <- -1e-2; 
    upper.vec[which(upper.tri(lambda))] <- 1e-2;
    q <- optim(lambda, method="L-BFGS-B", fn=lll, gr=grad.lambda, a=a, A=A, n=n, m=m, nlatent=nlatent, lower=lower.vec, upper=upper.vec, control=list(trace=0, fnscale=-1, maxit = 100, factr=1e-2))
    new.lambda <- matrix(q$par[1:(nlatent*m)], m, nlatent, byrow=FALSE)

    # Stage (3): Update A_i's.
    A[[1]] = solve(diag(rep(1,nlatent)) + diag(apply(lambda^2,2,sum),nlatent,nlatent)) 
    for (i in 2:n)
    {
      # Update parameters. Note these are all the same in our case.
      A[[i]] = A[[1]]
    }
    
    # Update the parameters.
    a <- new.a;
    lambda <- new.lambda;
    
    # Check if we are done.
    if (n_it > maxiters) break;
  }
  
  return(list(a=a, A=A, lambda=lambda));
}
