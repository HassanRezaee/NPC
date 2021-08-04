

# Program to compute the path distances between set of locations on an iGraph graph object

# Inputs:
 # igr: input iGraph graph object
 # data_ids: node IDs corresponding to data locations
 # knot_ids: node IDs corresponding to knot locations
 # nodes_coords: igr nodes coordinates

# Outputs: 
 # W: weight matrix of size nxm
 # D: distance matrix of size nxm



compute_D = function(igr, data_ids, knot_ids, nodes_coords){
  
  print("Computing path distances started")
  
  profvis::pause(1)
  
  # find the intersection degrees
  degrees = as.matrix(igraph::degree(igr, v = igraph::V(igr), mode = "total", loops = F, normalized = FALSE))
  
  
  n = length(data_ids)
  m = length(knot_ids)
  D = W = matrix(NA, n, m)
  
  
  for (i in 1:m) {
    
    print(paste("Computing path distances ", round(100*i/m,1), "% completed...") )
    
    
    paths_all = igraph::get.shortest.paths(igr, from=knot_ids[i], to = data_ids, mode = "all")
    
    spds_comp = matrix(0, n, 1);
    
    for (j in 1:n) {
      
      path = as.matrix(paths_all[["vpath"]][[j]])

      if (length(path)==1 & knot_ids[i]!=data_ids[j]){ # no path exists
        spds_comp[j,1] = Inf
      }else if (knot_ids[i]==data_ids[j]){ # data and knot locations are co-located
        spds_comp[j,1] = 0
      }else if (length(path)>1 & knot_ids[i]!=data_ids[j]) {
        for (jj in 1:(length(path)-1)){
          d = sqrt((nodes_coords[path[jj+1],1]-nodes_coords[path[jj],1])^2   +
                     (nodes_coords[path[jj+1],2]-nodes_coords[path[jj],2])^2)
          spds_comp[j,1] = spds_comp[j,1]+d
        }
      }

      D[j,i] = spds_comp[j,1]
      path_degrees = degrees[path]; path_degrees[path_degrees==1] = 2
      
      if (length(path)==1 & knot_ids[i]!=data_ids[j]){
        W[j,i] = 0
      }else if (knot_ids[i]==data_ids[j]){
        W[j,i] = 1
      }else if (length(path)>1 & knot_ids[i]!=data_ids[j]) {
        W[j,i] = sqrt(1/prod(path_degrees-1))
      }
      
    }
    
    
  }
  
  print("Computing path distances finished")
  
  # disconnected locations are specified with 0 weights
  W[is.infinite(D)] = 0
  D[is.infinite(D)] = Inf
  
  return(list(D,W))
  
}