
[![Build Status](https://travis-ci.org/AnthonyChristidis/nnGarrote.svg?branch=master)](https://travis-ci.com/AnthonyChristidis/nnGarrote) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/nnGarrote)](https://cran.r-project.org/package=nnGarrote) [![Downloads](http://cranlogs.r-pkg.org/badges/nnGarrote)](https://cran.r-project.org/package=nnGarrote)

nnGarrote
=====

This package provides functions to compute the non-negative garrote estimator with (or without) a penalized initial estimator.

------------------------------------------------------------------------

### Installation

You can install the **stable** version on [R CRAN](https://cran.r-project.org/package=nnGarrote).

``` r
install.packages("nnGarrote", dependencies = TRUE)
```

You can install the **development** version from [GitHub](https://github.com/AnthonyChristidis/nnGarrote).

``` r
library(devtools)
devtools::install_github("AnthonyChristidis/nnGarrote")
```

### Usage

Here is some code to compute the non-negative garrote estimator with ridge regression as an initial estimator, and compare it with ridge regression without the additional garrote shrinkage.

``` r
# Setting the parameters
p <- 500
n <- 100
n.test <- 5000
sparsity <- 0.15
rho <- 0.5
SNR <- 3
set.seed(0)
# Generating the coefficient
p.active <- floor(p*sparsity)
a <- 4*log(n)/sqrt(n)
neg.prob <- 0.2
nonzero.betas <- (-1)^(rbinom(p.active, 1, neg.prob))*(a + abs(rnorm(p.active)))
true.beta <- c(nonzero.betas, rep(0, p-p.active))
# Two groups correlation structure
Sigma.rho <- matrix(0, p, p)
Sigma.rho[1:p.active, 1:p.active] <- rho
diag(Sigma.rho) <- 1
sigma.epsilon <- as.numeric(sqrt((t(true.beta) %*% Sigma.rho %*% true.beta)/SNR))

# Simulate some data
x.train <- mvnfast::rmvn(n, mu=rep(0,p), sigma=Sigma.rho)
y.train <- 1 + x.train %*% true.beta + rnorm(n=n, mean=0, sd=sigma.epsilon)
x.test <- mvnfast::rmvn(n.test, mu=rep(0,p), sigma=Sigma.rho)
y.test <- 1 + x.test %*% true.beta + rnorm(n.test, sd=sigma.epsilon)

# Applying the NNG with Ridge as an initial estimator
nng.out <- cv.nnGarrote(x.train, y.train, intercept=TRUE,
                        initial.model=c("LS", "glmnet")[2],
                        lambda.nng=NULL, lambda.initial=NULL, alpha=0,
                        nfolds=5)
nng.predictions <- predict(nng.out, newx=x.test)
mean((nng.predictions-y.test)^2)/sigma.epsilon^2

# Ridge Regression
cv.ridge <- glmnet::cv.glmnet(x.train, y.train, alpha=0)
ridge <- glmnet::glmnet(x.train, y.train, alpha=0, lambda=cv.ridge$lambda.min)
ridge.predictions <- predict(ridge, newx=x.test)
mean((ridge.predictions-y.test)^2)/sigma.epsilon^2

# Comparisons of the coefficients
coef(nng.out)
coef(ridge)
```

Note that the prediction accuracy is nearly identical, but the non-negative garrote output for the coefficient is much closer to the true one than the ridge regression output.

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
