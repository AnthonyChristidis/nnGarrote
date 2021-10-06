# ------------------------
# Function - Create Folds
# ------------------------

create_folds <- function(n_data, nfolds = 10){

  folds <- rep(list(numeric(0)), nfolds)
  available_samples <- 1:n_data
  fold_size <- floor(n_data/nfolds)

  if(fold_size==n_data/nfolds){

    for(fold_id in 1:nfolds){

      folds[[fold_id]] <- sample(available_samples, fold_size)
      available_samples <- available_samples[!(available_samples %in% folds[[fold_id]])]
    }

  } else{

    for(fold_id in 1:(nfolds-1)){

      folds[[fold_id]] <- sample(available_samples, min(fold_size + rbinom(1, 1, 0.5), length(available_samples)))
      available_samples <- available_samples[!(available_samples %in% folds[[fold_id]])]
    }
    folds[[nfolds]] <- available_samples

  }

  return(folds)

}
