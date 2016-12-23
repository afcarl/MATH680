library(mvtnorm);
source('gllvm.R');

generate <- function(latent.true) 
{
  n <- dim(latent.true$lv)[c(1)];
  m <- dim(latent.true$coefs)[c(1)];
  y <- matrix(0, n, m);

  eta <- cbind(1, as.matrix(latent.true$lv))%*%t(latent.true$coefs);
  
  for (j in 1:m)
  {
    z.j <- rnorm(n, eta[,j], 1);
    y[which(z.j > 0),j] <- 1;
  }
  
  return(y);
  
}
# Number of responses per observation.
m <- 10;

# Generate data that is 0.5 from bivariate normal with mean c(-2,2), 
#                       0.3 from bivariate normal with mean c(0,-1),
#                       0.2 from bivariate normal with mean c(1,1).

# latent.true$lv is the realization of the latent variable.
# latent.true$label shows which bivariate normal it came from.
latent.true <- list(lv=rbind(rmvnorm(50,c(-2,2)), 
                              rmvnorm(30, c(0, -1)),
                              rmvnorm(20, c(1,1))), 
                    label=rep(1:3,c(25,15,10)*1));
# latent.true$coefs are the coefficients for each of the m responses plus an intercept.
latent.true$coefs <- cbind(runif(m, -1, 1), seq(-2,2,length=m),seq(1,-1,length=m))

y <- generate(latent.true);
sims <- replicate(1000, generate(latent.true));

fit <- bernoulli.fit(y, 2, 2, 0.01);

