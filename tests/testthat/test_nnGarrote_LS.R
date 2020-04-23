# -----------------------------------------------
# Test Script - Error nnGarrote with LS function
# -----------------------------------------------

# Required libraries
library(nnGarrote)
library(mvnfast)

# Context of test script
context("Ensure output is returned.")

# There should be an error if we want to compute the IF TS, and no returns are provided
test_that("The cross-validation function \"cv.nnGarrote\" is returning an output.", {

  # Setting the parameters
  p <- 50
  n <- 500
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
  x.train <- rmvn(n, mu=rep(0,p), sigma=Sigma.rho)
  y.train <- 1 + x.train %*% true.beta + rnorm(n=n, mean=0, sd=sigma.epsilon)
  x.test <- rmvn(n.test, mu=rep(0,p), sigma=Sigma.rho)
  y.test <- 1 + x.test %*% true.beta + rnorm(n.test, sd=sigma.epsilon)

  expect_output(
    nng.out <- cv.nnGarrote(x.train, y.train, intercept=TRUE,
                            initial.model=c("LS", "glmnet")[1],
                            lambda.nng=NULL, lambda.initial=NULL, alpha=0,
                            nfolds=5)
    )
})

