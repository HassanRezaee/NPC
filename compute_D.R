
compute_D = function(igr,data_ids,knot_ids,degrees, nodes_coords){
  print("Computing path distances started")
  profvis::pause(3)
  n = length(data_ids)
  m = length(knot_ids)
  D_ref = W_ref = matrix(NA,n,m); i=1
  for (i in 1:m) {
    
    print(paste("Computing path distances ", round(100*i/m,1), "% completed...") )
    
    
    paths_all = igraph::get.shortest.paths(igr, from=knot_ids[i], to = data_ids, mode = "all")
    
    neighbs_nb = length(data_ids)
    spds_comp = matrix(0,n,1); j = 5
    for (j in 1:neighbs_nb) {
      path = as.matrix(paths_all[["vpath"]][[j]])

      if (length(path)==1 & knot_ids[i]!=data_ids[j]){
        spds_comp[j,1] = Inf
      }else if (knot_ids[i]==data_ids[j]){
        spds_comp[j,1] = 0
      }else if (length(path)>1 & knot_ids[i]!=data_ids[j]) {
        for (jj in 1:(length(path)-1)){
          d = sqrt((nodes_coords[path[jj+1],1]-nodes_coords[path[jj],1])^2   +
                     (nodes_coords[path[jj+1],2]-nodes_coords[path[jj],2])^2)
          spds_comp[j,1] = spds_comp[j,1]+d
        }
      }

      D_ref[j,i] = spds_comp[j,1]
      path_degrees = degrees[path]; path_degrees[path_degrees==1] = 2
      
      if (length(path)==1 & knot_ids[i]!=data_ids[j]){
        W_ref[j,i] = 0
      }else if (knot_ids[i]==data_ids[j]){
        W_ref[j,i] = 1
      }else if (length(path)>1 & knot_ids[i]!=data_ids[j]) {
        W_ref[j,i] = sqrt(1/prod(path_degrees-1))
      }
      
    }
  }
  
  print("Computing path distances finished")
  
  W_ref[is.infinite(D_ref)] = 0
  D_ref[is.infinite(D_ref)] = Inf
  
  return(list(D_ref,W_ref))
  
}