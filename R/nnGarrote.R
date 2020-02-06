#'
#' @importFrom stats coef predict sd
#'
#' @title Non-negative Garrote Estimator
#'
#' @description \code{nnGarrote} computes the non-negative garrote estimator.
#'
#' @param x Design matrix.
#' @param y Response vector.
#' @param intercept Boolean variable to determine if there is intercept (default is TRUE) or not.
#' @param initial.model Model used for the groups. Must be one of "LS" (default) or "glmnet".
#' @param lambda.nng Shinkage parameter for the non-negative garrote. If NULL(default), it will be computed based on data.
#' @param lambda.initial The shinkrage parameter for the "glmnet" regularization. If NULL (default), optimal value is chosen by cross-validation.
#' @param alpha Elastic net mixing parameter for initial estimate. Should be between 0 (default) and 1.
#'
#' @return An object of class nnGarrote.
#'
#' @export
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @seealso \code{\link{coef.nnGarrote}}, \code{\link{predict.nnGarrote}}
#'
#' @examples
#' # Setting the parameters
#' p <- 500
#' n <- 100
#' n.test <- 5000
#' sparsity <- 0.15
#' rho <- 0.5
#' SNR <- 3
#' set.seed(0)
#' # Generating the coefficient
#' p.active <- floor(p*sparsity)
#' a <- 4*log(n)/sqrt(n)
#' neg.prob <- 0.2
#' nonzero.betas <- (-1)^(rbinom(p.active, 1, neg.prob))*(a + abs(rnorm(p.active)))
#' true.beta <- c(nonzero.betas, rep(0, p-p.active))
#' # Two groups correlation structure
#' Sigma.rho <- matrix(0, p, p)
#' Sigma.rho[1:p.active, 1:p.active] <- rho
#' diag(Sigma.rho) <- 1
#' sigma.epsilon <- as.numeric(sqrt((t(true.beta) %*% Sigma.rho %*% true.beta)/SNR))
#'
#' # Simulate some data
#' # x.train <- mvnfast::rmvn(n, mu=rep(0,p), sigma=Sigma.rho)
#' # y.train <- 1 + x.train %*% true.beta + rnorm(n=n, mean=0, sd=sigma.epsilon)
#' # x.test <- mvnfast::rmvn(n.test, mu=rep(0,p), sigma=Sigma.rho)
#' # y.test <- 1 + x.test %*% true.beta + rnorm(n.test, sd=sigma.epsilon)
#'
#' # Applying the NNG with Ridge as an initial estimator
#' # nng.out <- nnGarrote(x.train, y.train, intercept=TRUE,
#' #                      initial.model=c("LS", "glmnet")[2],
#' #                      lambda.nng=NULL, lambda.initial=NULL, alpha=0)
#' # nng.predictions <- predict(nng.out, newx=x.test)
#'
#' # Ridge Regression
#' # cv.ridge <- glmnet::cv.glmnet(x.train, y.train, alpha=0)
#' # ridge <- glmnet::glmnet(x.train, y.train, alpha=0, lambda=cv.ridge$lambda.min)
#' # ridge.predictions <- predict(ridge, newx=x.test)
#' # mean((ridge.predictions-y.test)^2)/sigma.epsilon^2
#'
#' # Comparisons of the coefficients
#' # coef(nng.out)
#' # coef(ridge)#'
nnGarrote <- function(x, y, intercept = TRUE,
                      initial.model = c("LS", "glmnet")[1],
                      lambda.nng = NULL, lambda.initial = NULL, alpha = 0){

  # Check input data
  if (all(!inherits(x, "matrix"), !inherits(x, "data.frame"))) {
    stop("x should belong to one of the following classes: matrix, data.frame.")
  } else if (all(!inherits(y, "matrix"), all(!inherits(y, "numeric")))) {
    stop("y should belong to one of the following classes: matrix, numeric.")
  } else if (any(anyNA(x), any(is.nan(x)), any(is.infinite(x)))) {
    stop("x should not have missing, infinite or nan values.")
  } else if (any(anyNA(y), any(is.nan(y)), any(is.infinite(y)))) {
    stop("y should not have missing, infinite or nan values.")
  } else {
    if(inherits(y, "matrix")) {
      if (ncol(y)>1){
        stop("y should be a vector.")
      }
      y <- as.numeric(y)
    }
    len_y <- length(y)
    if (len_y != nrow(x)) {
      stop("y and x should have the same number of rows.")
    }
  }
  if(!is.null(alpha)){
    if (!inherits(alpha, "numeric")) {
      stop("alpha should be numeric.")
    } else if (any(alpha < 0, alpha > 1)) {
      stop("alpha should be a numeric value between 0 and 1.")
    }
  }
  if(!is.null(lambda.nng)){
    if (!inherits(lambda.nng, "numeric")) {
      stop("lambda.nng should be numeric.")
    } else if (any(lambda.nng < 0)) {
      stop("lambda.nng should be a numeric non-negative vector.")
    }
  }
  if(!is.null(lambda.initial)){
    if (!inherits(lambda.initial, "numeric")) {
      stop("lambda.initial should be numeric.")
    } else if (any(lambda.initial < 0, length(lambda.initial)!=1)) {
      stop("lambda.initial should be a numeric non-negative vector of length 1.")
    }
  }
  if(!(initial.model %in% c("LS", "glmnet"))){
    stop("initial.model should be one of \"LS\" or \"glmnet\".")
  }

  # Centering and scaling data
  x.s <- scale(x, center=TRUE, scale=TRUE)
  y.s <- scale(y, center=TRUE, scale=TRUE)

  # Case where initial estimator is LS
  if(initial.model=="LS"){

    # Stop algorithm if LS intial estimate and p>n
    if(ncol(x.s)>nrow(x.s))
      stop("Case where p variables greater than n observations. Use \"initial.model=glmnet\" option if needed.")

    # Getting the initial shrinkage parameters
    initial.beta <- solve(t(x.s)%*%x.s)%*%t(x.s)%*%y.s
    z <- sapply(1:ncol(x.s),
                     function(x, x.s, beta) {return(x.s[,x, drop=FALSE]*beta[x])},
                     x.s=x.s, beta=initial.beta)

    # Applying the NNG
    z.fit <- glmnet::glmnet(z, y.s, alpha=alpha, intercept=FALSE, lower.limits=0)
    # Storing the lambda.nng vector
    if(is.null(lambda.nng))
      lambda.nng <- z.fit$lambda

  } else if(initial.model=="glmnet"){

    if(is.null(lambda.initial))
      initial.beta <- coef(glmnet::cv.glmnet(x.s, y.s, alpha=alpha), s="lambda.min") else
        initial.beta <- coef(glmnet::glmnet(x.s, y.s, alpha=alpha, lambda=lambda.initial))

    # Computing the z matrix
    z <- sapply(1:ncol(x.s),
                function(x, x.s, beta) {return(x.s[,x, drop=FALSE]*beta[x])},
                x.s=x.s, beta=initial.beta)

    # Applying the NNG
    z.fit <- glmnet::glmnet(z, y.s, alpha=alpha, intercept=FALSE, lower.limits=0, lambda=lambda.nng)
    # Storing the lambda.nng vector
    if(is.null(lambda.nng))
      lambda.nng <- z.fit$lambda
  }

  # Computing the fitted NNG vector for each value of the NNG
  nng.beta <- list()
  for(lambda.ind in 1:length(lambda.nng)){
    nng.beta[[lambda.ind]] <- sapply(1:ncol(x.s),
                                     function(x, beta, z.beta) {return(z.beta[x]*beta[x])},
                                     beta=initial.beta, z.beta=z.fit$beta[,lambda.ind], simplify = TRUE)
  }
  # Unlisting the betas
  nng.beta <- matrix(unlist(nng.beta), byrow=FALSE, nrow=ncol(x.s))

  # Scaling back to original scale
  if(intercept){

    nng.beta.temp <- matrix(nrow=nrow(nng.beta)+1, ncol=ncol(nng.beta))
    nng.intercept <- numeric(ncol(nng.beta))
    for(lambda.ind in 1:length(lambda.nng)){
      nng.beta[, lambda.ind] <- nng.beta[, lambda.ind] * (sd(y)/apply(x, 2, sd))
      # Computing the intercept
      nng.intercept[lambda.ind] <- as.numeric(mean(y) - apply(x, 2, mean)%*%nng.beta[,lambda.ind])
      nng.beta.temp[,lambda.ind] <- c(nng.intercept[lambda.ind], nng.beta[,lambda.ind])
    }
    nng.beta <- nng.beta.temp
    rm(nng.beta.temp)
  } else{

      for(lambda.ind in 1:length(lambda.ind))
        nng.beta[,lambda.ind] <- nng.beta[, lambda.ind] * (sd(y)/apply(x, 2, sd))
  }

  # Return the output
  nng.out <- list(betas=nng.beta, intercepts=nng.intercept, lambda.nng=lambda.nng)
  nng.out <- construct.nnGarrote(nng.out, fn_call=match.call())
  return(nng.out)
}








