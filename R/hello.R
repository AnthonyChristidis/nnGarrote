# # Storing the data for the NNG
# z <- matrix(nrow=nrow(x), ncol=ncol(x))
#
# # Computing the LS estimates for the variables
# for(g in G){
#   x.g <- x.s[, variables.split==g, drop=FALSE]
#   beta.g <- solve(t(x.g)%*%x.g)%*%t(x.g)%*%y
#   final.beta[variables.split==g] <- beta.g
#   z.g <- sapply(1:ncol(x.g), function(x, x.g, beta.g) {return(x.g[,x]*beta.g[x])}, x.g=x.g, beta.g=beta.g)
#   z[,variables.split==g] <- z.g
# }
#
# # Compute the NNG estimate
# nng.cv <- glmnet::cv.glmnet(z, y.s, alpha=1, intercept=FALSE, lower.limit=0)
# nng <- glmnet::glmnet(z, y.s, alpha=1, intercept=FALSE, lambda=nng.cv$lambda.min, lower.limit=0)
# nng.coef <- coef(nng)[-1] * final.beta * (sd(y)/apply(x, 2, sd))
# if(intercept){
#   nng.intercept <- as.numeric(mean(y) - apply(x, 2, mean)%*%nng.coef)
#   nng.coef <- c(nng.intercept, nng.coef)
# }
