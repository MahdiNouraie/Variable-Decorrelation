# This script reproduces the toy exa mple given in Sction 2
if (!requireNamespace("glmnet")) {install.packages("glmnet")} # install package if not already installed

# Load necessary libraries
library(glmnet)

#sessionInfo() # Check R version and loaded packages
#R version 4.4.2 (2024-10-31)
#Platform: aarch64-apple-darwin20
#Running under: macOS Sequoia 15.3.1

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
#  [1] glmnet_4.1-8 Matrix_1.7-1

#loaded via a namespace (and not attached):
#  [1] compiler_4.4.2    tools_4.4.2       survival_3.8-3    rstudioapi_0.17.1
#[5] Rcpp_1.0.14       splines_4.4.2     codetools_0.2-20  grid_4.4.2       
#[9] iterators_1.0.14  foreach_1.5.2     shape_1.4.6.1     lattice_0.22-6 

# Grahm-Schmidt QR decomposition function
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


# Set seed for reproducibility
set.seed(26)

# Create a correlated design matrix X
n <- 100  # Number of samples
p <- 5    # Number of predictors

# Generate correlated predictors
X1 <- rnorm(n)
X2 <- 0.9 * X1  # X2 is highly correlated with X1
X3 <- rnorm(n)
X4 <- rnorm(n)
X5 <- rnorm(n)

# Combine the predictors into a matrix
X <- cbind(X1, X2, X3, X4, X5)

# Create response variable y with some noise
beta <- c(1, 0, 0, 0, 0)  # Only X1 should be important
y <- X %*% beta + rnorm(n)

# Perform QR decomposition using Gram-Schmidt
QR1 <- grahm_schimdtR(X)
# Extract Q 
Q1 <- QR1$Q
# Perform cross-validation 
cv_fit1 <- cv.glmnet(Q1, y, alpha = 1)

# Get the lambda (lambda.1se)
lambda_opt1 <- cv_fit1$lambda.1se

# Apply Lasso using the lambda
lasso1 <- glmnet(Q1, y, alpha = 1, lambda = lambda_opt1)

# Change the order of columns in X
X_reordered <- X[, c(2, 1, 3, 4, 5)]  # Swap X1 and X2 (highly correlated columns)

# Perform QR decomposition using Gram-Schmidt on the reordered X
QR2 <- grahm_schimdtR(X_reordered)
# Extract Q and R matrices for the reordered X
Q2 <- QR2$Q
# Perform cross-validation to find lambda for reordered X
cv_fit2 <- cv.glmnet(Q2, y, alpha = 1)

# Get the lambda (lambda.1se) for the reordered X
lambda_opt2 <- cv_fit2$lambda.1se

# Apply Lasso using the lambda for reordered X
lasso2 <- glmnet(Q2, y, alpha = 1, lambda = lambda_opt2)

# Display coefficients for both models
cat("Coefficients for the first model (original order):\n")
print(coef(lasso1))

cat("\nCoefficients for the second model (reordered columns):\n")
print(coef(lasso2))
