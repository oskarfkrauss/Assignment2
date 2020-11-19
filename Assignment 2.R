## Q2 b)

install.packages('maxLik')

load('dataex2.Rdata')

df <- dataex2

l <- function(param, data){
  x <- df[,1] 
  r <- df[,2]
  mu <- param
  sum(r*dnorm(x, mean = mu, sd = 1.5^2, log = T) + (1-r)*pnorm(x, mean = mu, sd = 1.5^2, log = T))
}

require(maxLik)

mle_maxLik <- maxLik(logLik = l, data = df, start = 0)

mle_optim <- optim(0, l, control = list(fnscale = -1))

## Q4

load('dataex4.Rdata')

df <- dataex4

df <- df[order(df$Y, na.last = T),] 

ind_obs <- which(is.na(df$Y) == F)
ind_mis <- which(is.na(df$Y) == T)

x_obs <- df[ind_obs,1]
y_obs <- df[ind_obs,2]
x_mis <-df[ind_mis, 1]

Mstep<-function(beta_0, start, tol){
  
  diff <- 1
  beta_new <- beta_0
  
   while(diff > tol){
    
    print(beta_new)
    
    beta_old <- beta_new
    
    l2 <- function(param, data){
  
      x <- df[,1]
      x_obs <- df[ind_obs,1]
      y_obs <- df[ind_obs,2]
      x_mis <-df[ind_mis, 1]
      
      beta0 <- param[1]
      beta1 <- param[2]
  
      sum((beta0 + beta1*x_mis)*(exp(beta_old[1] + x_mis*beta_old[2]))/(1 + exp(beta_old[1] + x_mis*beta_old[2]))) - sum(log(1 + exp(beta0 + x*beta1))) + sum(y_obs*(beta0 + x_obs*beta1))
    }

  beta_new <- maxLik(logLik = l2, data = df, start = start)$estimate
  diff <- sum(abs(beta_old - beta_new))
  
  }
  
  return(beta_new)
  
}

Mstep(c(0, 0), c(0,0), 1e-6)

## Q5 b)

load('dataex5.Rdata')

y <- dataex5

EM <- function(y, theta0, tol){
  
  theta <- theta0
  
  p <- theta[1]
  mu <- theta[2]
  sigma2 <- theta[3] 
  lam <- theta[4]
  
  diff <- 1
  
  while(diff > tol){
    
  theta.old <- theta
  
  #E-step
  ptilde1 <- p*dlnorm(y, meanlog = mu, sdlog = sqrt(sigma2))
  ptilde2 <- (1-p)*dexp(y, rate = lam)
  
  ptilde <- ptilde1/(ptilde1 + ptilde2)
  
  #M-step
  p <- mean(ptilde)
  
  mu <- sum(log(y)*ptilde)/(sum(ptilde))
 
  sigma2 <- sum(ptilde*(log(y)-mu)^2)/(sum(ptilde))
  
  lam <- sum(1-ptilde)/sum(y*(1-ptilde))
  
  theta <- c(p, mu, sigma2, lam)
  
  diff <- sum(abs(theta-theta.old))
  }
  return(theta)
  print(theta)
}

theta0 <- c(0.1, 1, 0.5^2,2)

theta <- EM(y, theta0, 1e-6)

pdf(file = 'EM_model.pdf',width = 12,height = 6)
h <- hist(y, breaks = 20, freq = F, 
          col = 'pink', xlim = c(0, max(y)), main = '')

y_plot<-seq(0,max(y),length=10000)
fy_plot <- theta[1]*dlnorm(y_plot, meanlog = theta[2], sdlog = sqrt(theta[3])) + (1 - theta[1])*theta[4]*exp(-theta[4]*y_plot)

lines(y_plot, fy_plot, col="black", lwd=1.5)

dev.off()
