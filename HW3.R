#Statistical Learning_Homework3
#107064522

##### 1. #####
library(MASS)
set.seed(36)
n <- 100; sigma <- 5; beta0 <- c(2, -2, 0.5, 1, -3)
cormat <- diag(1, nrow = 5, ncol = 5); cormat[cormat == 0] <- 0.5
cholmat <- chol(cormat)
x <- matrix(rnorm(5*n, 0, 1), ncol = 5)%*%cholmat
err <- rnorm(n, 0, sigma)
y <- x%*%beta0+err  

##### 2.(a) #####
#beta_hat
beta_hat = ginv(t(x)%*%x)%*%t(x)%*%y
#beta_tilde_a
z = matrix(0, 100, 5)
for(j in 1:5){
  x_bar = mean(x[,j])
  sigma = sqrt(mean((x[,j]-x_bar)^2))
  z[,j] = (x[,j]-x_bar)/sigma 
}
y_tilde = y-mean(y)
beta_tilde = ginv(t(z)%*%z)%*%t(z)%*%y_tilde
#rescale beta_tilde back
beta_tilde_rescale = ginv(t(x)%*%x)%*%t(x)%*%(z%*%beta_tilde+mean(y))
#check
beta_tilde_rescale
beta_hat

##### 2.(b) #####
beta_t = matrix(0, 5, 16)
for(i in 1:16){
  lambda = 2^(i-11)
  beta_t[,i] = ginv(t(z)%*%z+2*n*lambda*diag(1,5,5))%*%t(z)%*%y_tilde
}

plot(2^c(-10:5), beta_t[1,],
     ylim = c(-4.5, 4.5),
     type = 'b',
     xlab = 'lambda',
     ylab = 'beta of lambda', 
     main = 'Solution Path',
     lwd = 2,
     col = 1)
for(i in 2:5){
  lines(2^c(-10:5), beta_t[i,], type = 'b',lwd = 2, col = i)
}
legend("topright", c("beta(1)","beta(2)","beta(3)","beta(4)","beta(5)"),
       lwd = 2,
       col = c(1,2,3,4,5))

#beta_hat(lambda=2)
beta_hat_2 = ginv(t(x)%*%x)%*%t(x)%*%(z%*%beta_t[,12]+mean(y))

##### 2.(c) #####
library(glmnet)
beta_glmnet = matrix(0, 5, 16)
sigma_y = sqrt(mean((y-mean(y))^2))
y_standard = y_tilde/sigma_y
for(i in 1:16){
  beta_ridge = glmnet(z, y_standard,
                      alpha = 0,
                      lambda = 2^(i-11)*2,
                      intercept = F,
                      standardize = F,
                      thresh = 1e-20)
  beta_glmnet[,i] = as.matrix(beta_ridge$beta*sigma_y)
}

#check
beta_glmnet
beta_t

##### 2.(d) #####
## Step 1 ## Rescale beta_tilde form 2(b) to get beta_hat
beta_h = matrix(0, 5, 16)
for(i in 1:16){
  beta_h[,i] = ginv(t(x)%*%x)%*%t(x)%*%(z%*%beta_t[,i]+mean(y))
}

MSE_h = matrix(0, 16, 1)
MSE_h_temp = matrix(0, 16, 1)
MSE_t = matrix(0, 16, 1)
MSE_t_temp = matrix(0, 16, 1)
for(i in 1:100){
  ## Step 2 ## Cteate test data
  set.seed(i)
  n <- 100; sigma <- 5; beta0 <- c(2, -2, 0.5, 1, -3)
  cormat <- diag(1, nrow = 5, ncol = 5); cormat[cormat == 0] <- 0.5
  cholmat <- chol(cormat)
  x_test <- matrix(rnorm(5*n, 0, 1), ncol = 5)%*%cholmat
  err <- rnorm(n, 0, sigma)
  y_test <- x_test%*%beta0+err 
  
  ## Step 3 ## Estimate y_test by beta_hat model
  est_y_h = x_test%*%beta_h
  ## Step 4 ## Calculate MSE of beta_hat
  for(j in 1:16){ 
    MSE_h_temp[j,] = mean((est_y_h[,j]-y_test)^2)
  }
  ## Step 5 ## Record the value of MSE
  MSE_h = MSE_h+MSE_h_temp
  
  ## Step 6 ## Standardization of test data
  z_test = matrix(0, 100, 5)
  for(j in 1:5){
    sigma = sqrt(mean((x_test[,j]-mean(x_test[,j]))^2))
    z_test[,j] = (x_test[,j]-mean(x_test[,j]))/sigma 
  }
  y_test_tilde = y_test-mean(y_test)
  
  ## Step 7 ## Calculate MSE of beta_tilde
  est_y_t = z_test%*%beta_h
  for(j in 1:16){ 
    MSE_t_temp[j,] = mean((est_y_t[,j]-y_test_tilde)^2)
  }
  MSE_t = MSE_t+MSE_t_temp
}

## Step 8 ## Plot MSE
plot(c(-10:5), MSE_h/100, type = 'b', xlab = 'lambda (log scale)', ylab = 'MSE', main = 'MSE', lwd = 2, col= 1)
lines(c(-10:5), MSE_t/100, type = 'b', lwd = 2, col= 2, lty = 2)
legend("topleft", c("MSE of beta_hat","MSE of beta_tilde"),
       lwd = 2,
       col = c(1,2),
       lty = c(1,2))
