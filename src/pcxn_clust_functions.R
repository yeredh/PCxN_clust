# ==== Functions ====
getTopNeighbors = function(path_name,pcxn,top_n=10){
  # get the top_n neighbors of path_name from pcxn.
  # the top neighbors based on the correlation magnitude |Pathcor|
  #
  # Args
  # path_name: character with pathway name
  # pcxn: a data frame with the pathway correlations, with
  #  the following columns: Pathway.A, Pathway.B, PathCor
  # top_n: number of top neighbors to retrieve
  #
  # Returns
  # a data frame (subset of pcxn) with the correlations 
  #   - between path_name and its top_n neighbors
  #   - the correlations between the top_n neighbors
  # neighbors of path_name
  neigh_ind =  pcxn$Pathway.A %in% path_name | pcxn$Pathway.B %in% path_name
  # return NULL if path_name is not present in pcxn
  if(sum(neigh_ind) == 0){
    warning(paste0(path_name," is not present in data frame provided."))
    return(NULL)
  }
  
  # select neighbors from significant correlation coefficients
  pcxn_neigh = pcxn[neigh_ind,]
  # order by magnitude
  pcxn_neigh = pcxn_neigh[order(abs(pcxn_neigh$PathCor),decreasing=T),]
  
  # switch pathway names: query pathway in Pathway.A, neighbors in Pathway.B
  sw_ind = pcxn_neigh$Pathway.B == path_name 
  pcxn_neigh$Pathway.B[sw_ind] = pcxn_neigh$Pathway.A[sw_ind]
  pcxn_neigh$Pathway.A[sw_ind] = path_name
  
  # get names of all neighbors
  tmp = pcxn_neigh$Pathway.B[!duplicated(pcxn_neigh$Pathway.B)]
  if(top_n > length(tmp)){
    warning(paste0("Only ",length(tmp)," neighbors available"))
    top_n = length(tmp)
  }
  # get top_n neighbors names
  top_neighbors = tmp[1:top_n]
  rm(tmp)
  
  # get correlations between top_n neighbors
  res_ind = pcxn$Pathway.A %in% top_neighbors & pcxn$Pathway.B %in% top_neighbors
  
  # return correlations between query pathway and top_n neighbors
  # AND correlations between top_n neighbors
  return(rbind(
    pcxn_neigh[1:top_n,],
    pcxn[res_ind,]
  ))
}


pcxn2cormat = function(dat){
  # transform the data frame into a correlation matrix
  #
  # Args
  # dat: a data frame with the pathway correlations, with
  #  the following columns: Pathway.A, Pathway.B, PathCor
  #
  # Returns
  # the corresponding correlation matrix
  
  # get names of pathways in network
  pathnames = unique(c(dat$Pathway.A,dat$Pathway.B))
  # matrix to store correlations 
  tmp_mat = matrix(0,ncol=length(pathnames),nrow=length(pathnames))
  dimnames(tmp_mat) = list(pathnames,pathnames)
  # get correlations from data frame (lower triangular entries)
  N=length(pathnames)
  for(j in 1:(N-1)){
    for(i in (j+1):N){
      # cat("(",i,",",j,")\n")
      tmp_ind = dat$Pathway.A %in% pathnames[c(i,j)] & dat$Pathway.B %in% pathnames[c(i,j)]
      if(sum(tmp_ind) == 1){
        tmp_mat[i,j] = dat$PathCor[tmp_ind]
      }
    }
  }
  # fill upper triangular entries
  tmp_mat = tmp_mat + t(tmp_mat)
  # add 1 to diagonal
  diag(tmp_mat) = 1
  return(tmp_mat)
}



pcxn2dist = function(dat){
  # transform the pathway correlation coeffcients into a distance matrix
  # to perform clustering
  #
  # Args
  # dat: a data frame with the pathway correlations, with
  #  the following columns: Pathway.A, Pathway.B, PathCor
  #
  # Returns
  # a dist object with the distance based on the pathway correlations
  # 1-|PathCor|
  
  # get names of pathways in network
  pathnames = unique(c(dat$Pathway.A,dat$Pathway.B))
  # matrix to store correlations 
  tmp_mat = matrix(NA,ncol=length(pathnames),nrow=length(pathnames))
  dimnames(tmp_mat) = list(pathnames,pathnames)
  # get correlations from data frame
  N=length(pathnames)
  for(j in 1:(N-1)){
    for(i in (j+1):N){
      # cat("(",i,",",j,")\n")
      tmp_ind = dat$Pathway.A %in% pathnames[c(i,j)] & dat$Pathway.B %in% pathnames[c(i,j)]
      if(sum(tmp_ind) == 1){
        tmp_mat[i,j] = dat$PathCor[tmp_ind]
      }else{
        tmp_mat[i,j] = 0
      }
    }
  }
  # get distance
  d_mat = as.dist(1-abs(tmp_mat))
  return(d_mat)
}