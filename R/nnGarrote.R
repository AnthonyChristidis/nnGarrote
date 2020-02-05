#'
#' @title Non-negative Garrote Estimator
#'
#' @description \code{nnGarrote} computes the non-negative garrote estimator.
#'
#' @param x Design matrix.
#' @param y Response vector.
#' @param intercept Boolean variable to determine if there is intercept (default is TRUE) or not.
#' @param initial.model Model used for the groups. Must be one of "LS" (default) or "glmnet".
#' @param lambda.nnGarrote Shinkage parameter for the non-negative garrote. If NULL(default), it will be computed based on data.
#' @param lambda.initial The shinkrage parameters for the "glmnet" regularization. If NULL (default), optimal value is chosen by cross-validation.
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
#'
#'
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
      stop("lambda.nng should be numeric.")
    } else if (any(lambda.initial < 0, length(lambda.initial)!=1)) {
      stop("lambda.initial should be a numeric non-negative value.")
    }
  }
  if(!(initial.model %in% c("LS", "glmnet"))){
    stop("initial.model should be one of \"LS\" or \"glmnet\".")
  }

  # Centering and scaling data
  x.s <- scale(x, center=TRUE, scale=TRUE)
  y.s <- scale(y, center=TRUE, scale=TRUE)

  # Case where initial estimator is LS
  if(initial.model="LS"){

    # Stop algorithm if LS intial estimate and p>n
    if(ncol(x.s)>nrow(x.s))
      stop("Case where p variables greater than n observations. Use \"glmnet\" option if needed.")

    # Getting the initial shrinkage parameters
    full.beta <- solve(t(x.s)%*%x.s)%*%t(x.s)%*%y.s
  }

}








