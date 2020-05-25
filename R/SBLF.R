
Rcsblf <- function(outpath) {
  .Call("csblf", as.character(outpath))
}

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

  #Call c routine, which writes results to tmp/Result/
  creturn = Rcsblf(temp_path)
  
  #read results into R from tmp/Result/
  results = vector(mode = 'list', length = length(list.files(result_path)))
  names(results) <- list.files(result_path)
  mapply(X = list.files(result_path), Y = 1:length(results),
    FUN = function(file, position) {
      results[[position]] <- as.matrix(read.table(paste0(result_path, X, sep = ''), quote="\"", comment.char=""))
    })
  
  
  # pred_train <- as.matrix(read.table(paste0(resultpath_temp, "PostMean_Out_train.txt", sep=""),
  #                                    quote="\"", comment.char=""))
  # pred_test <- as.matrix(read.table(paste0(resultpath_temp, "PostMean_Out_test.txt", sep=""),
  #                                   quote="\"", comment.char=""))
  # latent <- as.matrix(read.table(paste0(resultpath_temp, "PostMean_Latent.txt"), quote="\"", comment.char=""))
  # out <- list(pred_train=pred_train, pred_test=pred_test, latent=latent)

  
  # return results
  return(results)
}
