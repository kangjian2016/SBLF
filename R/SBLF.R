
Rcsblf <- function(outpath) {
  .Call("csblf", as.character(outpath))
}

SBLF <- function(xtrain, xtest, ztrain, ztest, voxel_loc) {
  # write to temp with proper path names
  temp_path = tempdir()
  
  if(file.exists(data_path)) {
    do.call(file.remove, list(list.files(data_path, full.names = TRUE)))
    file.remove(data_path)
  }
  if(file.exists(result_path)) {
    do.call(file.remove, list(list.files(result_path, full.names = TRUE)))
    file.remove(result_path)
  }
  
  data_path = paste(temp_path, 'Data/', sep = '/')
  result_path = paste(temp_path, 'Result/', sep = '/') 
  
  dir.create(data_path, result_path)
  
  write.table(x = xtrain, file = paste(data_path, 'xtrain.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  write.table(x = xtest, file = paste(data_path, 'xtest.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  write.table(x = ztrain, file = paste(data_path, 'ztrain.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  write.table(x = ztest, file = paste(data_path, 'ztest.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  write.table(x = voxel_loc, file = paste(data_path, 'voxel_loc.txt', sep = '/'), row.names = FALSE, col.names = FALSE)
  
  # figure out how to call c code referencing these paths
  result = Rcsblf(temp_path)
  
  # read results into return list
  
  
  # return results
  
}
