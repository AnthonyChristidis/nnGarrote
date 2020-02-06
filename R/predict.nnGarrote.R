#'
#' @title Predictions for nnGarrote Object
#'
#' @description \code{predict.nnGarrote} returns the prediction for nnGarrote for new data.
#'
#' @param object An object of class nnGarrote
#' @param newx A matrix with the new data.
#' @param ... Additional arguments for compatibility.
#'
#' @return A matrix with the predictions of the \code{nnGarrote} object.
#'
#' @export
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
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
#' # coef(ridge)
#'
#' @seealso \code{\link{nnGarrote}}
#'
predict.nnGarrote <- function(object, newx, ...){

  # Check input data
  if(!any(class(object) %in% "nnGarrote"))
    stop("The object should be of class \"nnGarrote\"")
  # Storing the number of variables
  if(is.null(object$intercepts))
    p <- nrow(object$betas) else
      p <- nrow(object$betas)-1
  if(is.matrix(newx)){
    if(ncol(newx)!=p)
      stop("The dimension of newx is invalid.")
  } else if(length(newx)!=p)
    stop("The number of variables for newx is invalid.")

  # Matrix to store the predictions
  predictions <- matrix(nrow=nrow(newx), ncol=ncol(object$betas))

  # Removing the intercepts
  if(!is.null(object$intercepts))
    object$betas <- object$betas[-1,,drop=FALSE]

  # Computing the predictions
  for(newx.id in 1:nrow(newx)){
    for(beta.id in 1:ncol(object$betas)){
      predictions[newx.id, beta.id] <- newx[newx.id,,drop=FALSE] %*% object$betas[,beta.id, drop=FALSE]
    }
  }

  # Adding the intercepts
  if(!is.null(object$intercepts))
    for(beta.id in 1:ncol(object$betas))
      predictions[,beta.id] <- predictions[,beta.id] + object$intercepts[beta.id]

  # Returning the coefficients
  return(predictions)
}








