# -----------------------------------------------------------------------
# Object Construction for nnGarrote object
#
# object: the nnGarrote object
# fn_call: the function call
construct.nnGarrote <- function(object, fn_call){
  class(object) <- append("nnGarrote", class(object))
  object$call <- fn_call
  return(object)
}

# -----------------------------------------------------------------------
# Object Construction for cv.nnGarrote object
#
# object: the cv.nnGarrote object
# fn_call: the function call
construct.cv.nnGarrote <- function(object, fn_call){
  class(object) <- append("cv.nnGarrote", class(object))
  object$call <- fn_call
  return(object)
}
