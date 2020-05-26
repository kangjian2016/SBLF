
Rcsblf <- function(outpath) {
  .Call("csblf", as.character(outpath))
}

#' Fit Spatial Bayesian Latent Factor Models
#'
#' Use this function to fit spatial Bayesian latent factor models given
#' the prepared input data sets. Please also comment on the important
#' output the user should be expecting.
#' 
#' The parameters prefixed with "x" are imaging predictor data, and
#' the parameters prefixed with "z" are imaging outcomes.
#' 
#' @param xtrain training data set of imaging predictors. Sample size is `ntrain * (P * image_len)`
#' @param xtest test data set of imaging predictors. sample size is `(ntotal-ntrain) * (P * image_len)`
#' @param ztrain training data set of imaging outcomes. Sample size is `(ntrain * image_len)`
#' @param ztest test data set of imaging outcomes. Sample size is `((ntotal-ntrain) * image_len)`
#' @param voxel_loc matrix of voxel coordinates.
#'
#' @details Specify important details to the user about the use
#' of this function. Please comment.
#' 
#' The function expects all data to be entered as matrices with the dimensions
#' outlined in each of the parameter arguments.
#' 
#' @references Please add your paper reference.
#' 
#' @return A list which contains
#' \itemize{
#'  \item{"pred_train"}{}
#'  \item{"pred_test"}{}
#'  \item{"latent"}{}
#'  \item{"draws"}{}
#' }
#'
#' @examples 
#' \donttest{
#' 
#' }
#'
#' @export

SBLF <- function(xtrain, xtest, ztrain, ztest, voxel_loc) {
  
  temp_path = paste(tempdir(), '/', sep = '')
  
  data_path = paste(temp_path, 'Data/', sep = '')
  result_path = paste(temp_path, 'Result/', sep = '')
  
  if(file.exists(data_path)) {
    do.call(file.remove, list(list.files(data_path, full.names = TRUE)))
    file.remove(data_path)
  }
  if(file.exists(result_path)) {
    do.call(file.remove, list(list.files(result_path, full.names = TRUE)))
    file.remove(result_path)
  }

  dir.create(data_path)
  dir.create(result_path)

  write.table(x = xtrain, file = paste(data_path, 'ROI_dat_mat.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  write.table(x = xtest, file = paste(data_path, 'ROI_dat_mat_test.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  write.table(x = ztrain, file = paste(data_path, 'ROI_task.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  write.table(x = ztest, file = paste(data_path, 'ROI_task_test.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  write.table(x = voxel_loc, file = paste(data_path, 'ROI_axes.txt', sep = '/'), row.names = FALSE, col.names = FALSE)

  image_len <- ncol(ztrain) # number of voxels per image
  write.table(x = image_len, file = paste(data_path, 'ROI_count.txt', sep = '/'), row.names = FALSE, col.names = FALSE)

  n_predictors <- ncol(xtrain) / image_len # number of imaging predictors
  n_train <- nrow(xtrain)
  n_test <- nrow(xtest)
  ROI_sizes <- matrix(data = c(1, n_train, image_len, n_predictors, n_test), ncol = 1)
  write.table(x = ROI_sizes, file = paste(data_path, 'ROI_sizes.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  cat('test')
  #Call c routine, which writes results to tmp/Result/
  creturn = Rcsblf(temp_path)
  
  #read results into R from tmp/Result/
  results = vector(mode = 'list', length = length(list.files(result_path)))
  names(results) <- list.files(result_path)
  X = list.files(result_path)
  Xnames = gsub('.txt', '', X)
  
  results = 
  lapply(X = X,
    FUN = function(X) {
      as.matrix(read.table(paste0(result_path, X, sep = ''), quote="\"", comment.char=""))
    })
  names(results) <- Xnames
  
  pred_train = results[['PostMean_Out_train']]
  pred_test = results[['PostMean_Out_test']]
  latent = results[['PostMean_Latent']]
  draws = results
  
  err_train = pred_train - ztrain
  mse_train = mean(err_train^2)  
  
  err_test <- pred_test - ztest
  mspe_test <- mean(err_test^2)
  
  err <- matrix(c(mse_train, mspe_test), ncol=2)
  colnames(err) <- c("MSE(training)", "MSPE(test)")
  
  # return results
  return(list(
    err = err, 
    pred_train = pred_train, 
    pred_test = pred_test, 
    latent = latent
  ))
  
}
