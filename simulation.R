library(mvtnorm);
library(vegan);
source('gllvm.R');
# Simulate data according to the probit norm.
generate <- function(latent.true) 
{
  n <- dim(latent.true$lv)[c(1)];
  m <- dim(latent.true$coefs)[c(1)];
  y <- matrix(0, n, m);
  
  eta <- cbind(1, as.matrix(latent.true$lv))%*%t(latent.true$coefs);
  # Generate normal data with mean eta then check which is greater than 0.
  for (j in 1:m)
  {
    z.j <- rnorm(n, eta[,j], 1);
    y[which(z.j > 0),j] <- 1;
  }
  
  return(y);
}
set.seed(680);
# Number of responses per observation.
m <- 10;
n <- 50;
# Generate data that is 0.5 from bivariate normal with mean c(-2,2), 
#                       0.3 from bivariate normal with mean c(0,-1),
#                       0.2 from bivariate normal with mean c(1,1).

# latent.true$lv is the realization of the latent variable.
# latent.true$label shows which bivariate normal it came from.
latent.true <- list(lv=rbind(rmvnorm(25,c(-2,2)), 
                             rmvnorm(15, c(0, -1)),
                             rmvnorm(10, c(1,1))), 
                    label=rep(1:3,c(25,15,10)));

# latent.true$coefs are the coefficients for each of the m responses plus an intercept.
latent.true$coefs <- cbind(0, seq(-1,1,length=m),seq(-2,2,length=m))
# y is the response matrix we care to model.
y <- generate(latent.true);

fit <- gllvm.fit(y, 2, 500);
# Check the gradients.
grad.l <- grad.lambda(fit$lambda, a=fit$a, A=fit$A, n=50, m=10, nlatent=2)
grad.v <- grad.a(fit$a, lambda=fit$lambda, A=fit$A, n=50, m=10, nlatent=2)
# Check the Procrustes error.
err.lv <- procrustes(fit$a+rnorm(n,0,1e-4), latent.true$lv+rnorm(n,0,1e-4), symmetric=T)$ss
err.lam <- procrustes(fit$lambda+rnorm(m,0,1e-4), latent.true$coefs[,2:3] + rnorm(m,0,1e-4), symmetric = T)$ss 
print(err.lv);
print(err.lam);