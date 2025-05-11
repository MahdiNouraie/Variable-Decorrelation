# This script reproduces Figure 7(a) of the paper
options(warn=-1) # Turn off warnings
if (!requireNamespace("hdi")) {install.packages("hdi")} # install package if not already installed
if (!requireNamespace("glmnet")) {install.packages("glmnet")} # install package if not already installed
if (!requireNamespace("latex2exp")) {install.packages("latex2exp")} # install package if not already installed
if (!requireNamespace("ggplot2")) {install.packages("ggplot2")} # install package if not already installed
if (!requireNamespace("cmna")) {install.packages("cmna")} # install package if not already installed


library(hdi) # Load the hdi package for riboflavin data
library(glmnet) # Load the glmnet package for Lasso
library(latex2exp) # Load the latex2exp package for TeX formatting
library(ggplot2) # Load the ggplot2 package for plotting
library(cmna) # Load the cmna package for newtons method

# sessionInfo() # Display session information
#R version 4.4.2 (2024-10-31)
#Platform: aarch64-apple-darwin20
#Running under: macOS Sequoia 15.4.1

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#time zone: Australia/Sydney
#tzcode source: internal

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] cmna_1.0.5      ggplot2_3.5.2   latex2exp_0.9.6 glmnet_4.1-8    Matrix_1.7-1   
#[6] hdi_0.1-9       scalreg_1.0.1   lars_1.3       

#loaded via a namespace (and not attached):
#  [1] gtable_0.3.6       dplyr_1.1.4        compiler_4.4.2     tidyselect_1.2.1  
#[5] Rcpp_1.0.14        stringr_1.5.1      parallel_4.4.2     splines_4.4.2     
#[9] scales_1.4.0       lattice_0.22-6     R6_2.6.1           linprog_0.9-4     
#[13] generics_0.1.3     shape_1.4.6.1      iterators_1.0.14   MASS_7.3-64       
#[17] tibble_3.2.1       pillar_1.10.2      RColorBrewer_1.1-3 rlang_1.1.6       
#[21] stringi_1.8.4      cli_3.6.5          withr_3.0.2        magrittr_2.0.3    
#[25] foreach_1.5.2      grid_4.4.2         rstudioapi_0.17.1  lifecycle_1.0.4   
#[29] vctrs_0.6.5        lpSolve_5.6.23     glue_1.8.0         farver_2.1.2      
#[33] codetools_0.2-20   survival_3.8-3     tools_4.4.2        pkgconfig_2.0.3  

# AirHOLP (2025) from "https://github.com/Logic314/Air-HOLP"
AirHOLP <- function(X, y, Threshold, r0 = 10, adapt = TRUE,
                    iter = 10, Lambda, U, XU) {
  # Arguments:-
  # X: matrix of features (Matrix)
  # y: response vector (Vector)
  # Threshold: screening threshold (Integer)
  # r0: initial penalties (Vector)
  # adapt: if >= 1 adaptive penalty will be used (Binary)
  # iter: maximum number of iterations for adaptive penalty selection (Integer)
  # Lambda: eigenvalues of XXT, if missing the function will compute it (Vector)
  # U: eigenvectors of XXT, if missing the function will compute it (Matrix)
  # XU: X transpose times U, if missing the function will compute it (Matrix)
  
  # Output:-
  # index_r: ranking of features by Air-HOLP (Matrix)
  # index_r0: ranking of features by Ridge-HOLP (Matrix)
  # Beta_r: regression coefficients of Air-HOLP (Matrix)
  # Beta_r0: regression coefficients of Ridge-HOLP (Matrix)
  # r: selected penalty parameters by Air-HOLP (Vector)
  # iter_last: number of iterations used in Air-HOLP (Vector)
  
  n <- dim(X)[1] # sample size
  p <- dim(X)[2] # number of features
  q <- length(r0) # number of penalty parameters
  iter_temp2 <- 0*(1:q) # used for calculating iter_last
  iter_temp1 <- iter_temp2 - 1 # used for calculating iter_last
  
  # Standardizing X and y:
  X <- X - matrix(rep(colMeans(X),each = n),n,p)
  X <- X/matrix(rep(sqrt(colMeans(X^2)),each = n),n,p)
  y <- (y - mean(y))/sd(y)
  
  if(adapt){
    # Main computations:
    if(missing(Lambda)|missing(U)){
      XXT <- tcrossprod(X)
      eXXT <- eigen(XXT)
      Lambda <- eXXT$values
      U <- eXXT$vectors
    }
    if(missing(XU)){
      XU <- crossprod(X,U)
    }
    Dn <- diag(Lambda)
    UTy <- crossprod(U,y)
    yUD2UTy <- UTy^2*(Lambda^2)
    
    # Penalty selection:
    r_max <- 1000*sqrt(n) # maximum penalty
    max.iter <- 30 # maximum number of iterations for Newtons method
    index_r <- matrix(1:(p*q), nrow = p, ncol = q)
    index_r0 <- index_r
    Beta_r <- index_r
    Beta_r0 <- index_r
    r <- r0
    r_temp <- r0
    for (j in 1:iter) {
      for (i in 1:q) {
        # Initial screening:
        Beta_temp <- XU%*%((Lambda+r[i])^(-1)*UTy)
        index_temp <- match(1:p,rank(-abs(Beta_temp), na.last = NA,
                                     ties.method = c("random")))
        Xs <- X[,index_temp[1:Threshold]] # Screened features
        if(j<2) {
          Beta_r0[,i] <- Beta_temp
          index_r0[,i] <- rank(-abs(Beta_temp), na.last = NA,
                               ties.method = c("random"))
        }
        
        # Estimating the expected response:
        ys <- Xs%*%(solve(crossprod(Xs) +
                            diag(Threshold)*10^-12)%*%crossprod(Xs,y))
        
        # MSE functions:
        ysUDUTy <- t(crossprod(ys,U)*Lambda)*UTy
        Z <- function(lam) { # The function we minimize
          t((Lambda+lam)^-2)%*%yUD2UTy - 2*t((Lambda+lam)^-1)%*%ysUDUTy
        }
        Z1 <- function(lam) { # First derivative
          -2*t((Lambda+lam)^-3)%*%yUD2UTy + 2*t((Lambda+lam)^-2)%*%ysUDUTy
        }
        Z2 <- function(lam) { # Second derivative
          6*t((Lambda+lam)^-4)%*%yUD2UTy - 4*t((Lambda+lam)^-3)%*%ysUDUTy
        }
        
        # MSE minimization:
        sol <- newton(Z1, Z2, 0.0001, tol = 0.001, m = max.iter)
        r[i] <- sol
        if(r[i] > r_max) {r[i] <- r_max}
        if(r[i] < 0.0001) {r[i] <- 0.0001}
        if(Z(r_max) < Z(r[i])) {r[i] <- r_max} # Checking boundaries
        if(Z(0.0001) < Z(r[i])) {r[i] <- 0.0001}
        
        # Feature screening:
        Beta_r[,i] <- XU%*%((Lambda+r[i])^(-1)*UTy)
        index_r[,i] <- rank(-abs(Beta_r[,i]), na.last = NA,
                            ties.method = c("random"))
        
        # Calculations for the number of iterations:
        if(abs(r[i] - r_temp[i]) < 0.01*r[i]){ # Checking relative error
          iter_temp1[i] <- j
          iter_temp2[i] <- iter_temp2[i] + 1
        }
      }
      if(sum(abs(r - r_temp) < 0.01*r) == q){ # Checking relative error
        break
      }
      r_temp <- r
    }
    iter_last <- iter_temp1 - iter_temp2 + 1 # Number of iterations
    AirHOLP <- list(index_r = index_r, index_r0 = index_r0, Beta_r = Beta_r,
                    Beta_r0 = Beta_r0, r = r, iter_last = iter_last)
  } else{
    if(q < 2) {
      # Feature screening:
      if(missing(Lambda)|missing(U)){
        Beta_r0 <- crossprod(X, solve(tcrossprod(X)+r0*diag(n),y))
      } else{
        UTy <- crossprod(U,y)
        Beta_r0 <- XU%*%((Lambda+r0)^(-1)*UTy)
      }
      index_r0 <- rank(-abs(Beta_r0), na.last = NA, ties.method = c("random"))
      AirHOLP <- list(index_r0 = index_r0, Beta_r0 = Beta_r0)
    } else{
      # Main computations:
      if(missing(Lambda)|missing(U)){
        XXT <- tcrossprod(X)
        eXXT <- eigen(XXT)
        Lambda <- eXXT$values
        U <- eXXT$vectors
      }
      if(missing(XU)){
        XU <- crossprod(X,U)
      }
      Dn <- diag(Lambda)
      UTy <- crossprod(U,y)
      
      # Feature screening:
      index_r <- matrix(1:(p*q), nrow = p, ncol = q)
      for (i in 1:q) {
        Beta_r0[,i] <- XU%*%((Lambda+r0[i])^(-1)*UTy)
        index_r0[,i] <- rank(-abs(Beta_r0[,i]), na.last = NA,
                             ties.method = c("random"))
      }
      AirHOLP <- list(index_r0 = index_r0, Beta_r0 = Beta_r0)
    }
  }
}

# Grahm-Schmidt QR decomposition from https://stackoverflow.com/questions/15584221/gram-schmidt-with-r
grahm_schimdtR <- function(A) {
  m <- nrow(A)
  n <- ncol(A)
  Q <- matrix(0, nrow = m, ncol = n)
  R <- matrix(0, nrow = n, ncol = n)
  for (j in 1:n) {
    v <- A[ , j, drop = FALSE]
    if (j > 1) {
      for(i in 1:(j-1)) {
        R[i, j] <- t(Q[,i,drop = FALSE]) %*% A[ , j, drop = FALSE]
        v <- v - R[i, j] * Q[ ,i]
      }
    }
    R[j, j] = norm(v, type = "2")
    Q[ ,j] = v / R[j, j]
  }
  
  list("Q" = Q, "R" = R)
  
}


# Stability measure (2018) from "https://github.com/nogueirs/JMLR2018/blob/master/R/getStability.R"
getStability <- function(X,alpha=0.05) {
  ## the input X is a binary matrix of size M*d where:
  ## M is the number of bootstrap replicates
  ## d is the total number of features
  ## alpha is the level of significance (e.g. if alpha=0.05, we will get 95% confidence intervals)
  ## it's an optional argument and is set to 5% by default
  ### first we compute the stability
  
  M<-nrow(X)
  d<-ncol(X)
  hatPF<-colMeans(X)
  kbar<-sum(hatPF)
  v_rand=(kbar/d)*(1-kbar/d)
  stability<-1-(M/(M-1))*mean(hatPF*(1-hatPF))/v_rand ## this is the stability estimate
  
  ## then we compute the variance of the estimate
  ki<-rowSums(X)
  phi_i<-rep(0,M)
  for(i in 1:M){ 
    phi_i[i]<-(1/v_rand)*((1/d)*sum(X[i,]*hatPF)-(ki[i]*kbar)/d^2-(stability/2)*((2*kbar*ki[i])/d^2-ki[i]/d-kbar/d+1))
  }
  phi_bar=mean(phi_i)
  var_stab=(4/M^2)*sum((phi_i-phi_bar)^2) ## this is the variance of the stability estimate
  
  ## then we calculate lower and upper limits of the confidence intervals
  z<-qnorm(1-alpha/2) # this is the standard normal cumulative inverse at a level 1-alpha/2
  upper<-stability+z*sqrt(var_stab) ## the upper bound of the (1-alpha) confidence interval
  lower<-stability-z*sqrt(var_stab) ## the lower bound of the (1-alpha) confidence interval
  
  return(list("stability"=stability,"variance"=var_stab,"lower"=lower,"upper"=upper))
  
}

set.seed(26) # Set seed for reproducibility

data(riboflavin) # Load the riboflavin data
data <- as.data.frame(cbind(Y=riboflavin$y - 1, X=riboflavin$x)) # Convert the data to a data frame
rm(riboflavin) # Remove the original data to save memory
x <- as.matrix(data[, -1]) # Extract the predictors
y <- data[,1] # Extract the response

p <- ncol(x) # Number of predictors
B <- 200 # Number of subsamples
x <- scale(x)  # Standardize the predictors
y <- scale(y, scale = FALSE) # Center the response

Threshold <- 10 # Screening threshold
AHOLP <- AirHOLP(x, y, Threshold = Threshold, r0 = 10, adapt = TRUE, iter = 10) # Perform AirHOLP
ranked_features <- AHOLP$index_r  # Ranking of features
# Convert ranks to column order
column_order <- order(ranked_features)  # Find the order that sorts ranks from 1 to p

x_ranked <- x[, column_order] # Reorder columns of X
names <- colnames(x_ranked) # Store the names of the predictors
gram <- grahm_schimdtR(x_ranked) # Perform Gram-Schmidt orthogonalization
x <- gram$Q # Extract the orthogonalized matrix
colnames(x) <- names # Set the column names of the orthogonalized matrix
cv_lasso <- cv.glmnet(x, y, nfolds = 10, alpha = 1) # Perform cross-validation for Lasso

candidate_set <- cv_lasso$lambda # Candidate set of lambda values

S_list <- vector("list", length(candidate_set)) # Initialize a list to store the selection matrices
names(S_list) <- paste0("lambda_", seq_along(candidate_set)) # Name the list entries

# Stability Selection for each lambda in candidate_set
for (lambda_idx in seq_along(candidate_set)) {
  
  lambda <- candidate_set[lambda_idx]  # Current lambda value
  S <- matrix(0, nrow = B, ncol = p)  # Initialize selection matrix for the current lambda
  colnames(S) <- colnames(x) # Set the column names of S to the predictors
  
  for (i in 1:B) {
    # Subsample the data
    sample_index <- sample(1:nrow(data), nrow(data) / 2, replace = FALSE) 
    
    # Prepare the response and predictors
    x_sub <- x[sample_index, ] # Predictors
    y_sub <- y[sample_index] # Response
    
    # Fit the Lasso model with the current lambda
    lasso_model <- glmnet(x_sub, y_sub, alpha = 1, lambda = lambda)
    
    # Extract significant predictors (ignoring the intercept, hence [-1])
    significant_predictors <- ifelse(coef(lasso_model) != 0, 1, 0)[-1]
    
    # Store the significant predictors in matrix S
    S[i, ] <- significant_predictors
  }
  
  # Store the matrix S for the current lambda in the corresponding list entry
  S_list[[lambda_idx]] <- S
}

stability_results <- lapply(S_list, function(S) { # Loop through regularization set and compute stability values
  getStability(S) # Compute stability values
})


stab_values <- c() # Initialize an empty vector to store stability values of each lambda
for (j in 1:length(candidate_set)) {
  temp <- stability_results[[paste0("lambda_", j)]]$stability # Extract stability values
  stab_values <- c(stab_values, temp) # Append stability values to the vector
}


stable_values <- which(stab_values > 0.75) # Index of stable lambda values
lambda_stable <- min(candidate_set[stable_values]) # Minimum stable lambda value
index_of_lambda_stable <- which(candidate_set == lambda_stable) # Index of lambda_stable

S_stable <- S_list[[index_of_lambda_stable]] # Extract the selection matrix of lambda_stable
stability <- data.frame() # Initialize an empty data frame to store stability values
for (k in 2:nrow(S_stable)){ # Loop through sub-samples results for lambda_stable
  output <- getStability(S_stable[1:k,]) # Compute stability values
  stability <- rbind(stability, data.frame(k, output$stability, output$variance, output$lower, output$upper)) # Append stability values to the data frame
}
colnames(stability) <- c('Iteration', 'Stability', 'Variance', 'Lower', 'Upper') # Set column names of the data frame


# Plot the stability values with confidence intervals
ggplot(stability, aes(x = Iteration, y = Stability)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = 'blue', alpha = 0.7) + # Add ribbon for confidence interval
  labs(title = TeX('Stability of Stability Selection ($\\lambda = \\lambda_{stable}$)'), 
       x = 'Iteration (sub-sample)', y = TeX('Stability ($\\hat{\\Phi}$)'))+
  theme_bw() +
  theme(
    plot.title = element_text(size = 18),       # Title text size
    axis.title.x = element_text(size = 18),     # X-axis label size
    axis.title.y = element_text(size = 18),     # Y-axis label size
    axis.text.x = element_text(size = 16),      # X-axis tick text size
    axis.text.y = element_text(size = 16)       # Y-axis tick text size
  )

# Calculate selection frequencies
col_means <- colMeans(S_stable)

# Filter columns with selection frequencies greater than 0.6
filtered_cols <- col_means[col_means > 0.6]

# Create a data frame with the column names and their means
stable_vars <- data.frame(Column = names(filtered_cols), Mean = filtered_cols)

# Display the result
print(stable_vars)
