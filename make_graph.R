

# Program to create an igraph graph object from the discretized road shape files

# Inputs:
 # roads_shp: road network shapefiles
 # delta_l: discretization unit length

# Outputs: 
 # igr: igraph graph object
 # nodes_coords_fine: coordinates of all the nodes in igr



make_graph = function(roads_shp, delta_l){
  
  roads = roads_shp
  roads = SplitLines(roads, split_length = delta_l, return.dataframe = F, plot.results = F)
  Sl_IDs <- sapply(slot(roads, "lines"), function(x) slot(x, "ID"))
  data = data.frame(var=seq(from=1,to=length(roads),by=1))
  rownames(data) = Sl_IDs
  roads <- sp::SpatialLinesDataFrame(sl = roads, data = data, match.ID = T)
  ea.prop = rep(1,length(roads)-1)
  
  network = shp2graph::readshpnw(roads, ELComputed=T, longlat=T, Detailed=F, ea.prop=ea.prop)
  igr <- shp2graph::nel2igraph(network[[2]], network[[3]], weight=network[[4]],Directed = TRUE)
  
  nodelist = network[[2]]; edgelist = network[[3]]
  nodes_coords_fine = shp2graph::Nodes.coordinates(nodelist)

  return(list(igr,nodes_coords_fine))
}