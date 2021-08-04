
# Program to discretize road shape files

# Inputs:
 # spatial_line: road shape files
 # split_length: discretization unit length


# Outputs: 
 # roads: discretized road shape files



SplitLines = function(spatial_line,
                      split_length,
                      return.dataframe = F,
                      plot.results = F) {
  
  require(sp)
  #### Define support functions ####
  # SpatialLines2df extracts start and end point coordinates of each segment of a SpatialLine object
  # spatial_line: an object class SpatialLinesDataFrame of the package sp
  SpatialLines2df = function(spatial_line) {
    df = data.frame(
      id = character(),
      mline_id = character(),
      segment_id = character(),
      fx = numeric(),
      fy = numeric(),
      tx = numeric(),
      ty = numeric(),
      stringsAsFactors = FALSE
    )
    for (i in 1:length(spatial_line)) {
      #print(100*i/length(spatial_line))
      coords = spatial_line@lines[[i]]@Lines[[1]]@coords # For each line takes the coords of the vertex
      row_nums = 1:(nrow(coords) - 1)
      mline_id = formatC(i, width = 9, flag = '0') # Creates id for the line
      segment_id = formatC(row_nums, width = 3, flag = '0') # Creates id for each single segment belonging to the line
      id = paste0(mline_id, '_', segment_id) # Creates a composite id
      for (j in row_nums) {
        # For each segment stores ids and coordinates
        df[nrow(df) + 1, ] = c(id[j],
                               mline_id,
                               segment_id[j],
                               coords[j, 1],
                               coords[j, 2],
                               coords[j + 1, 1],
                               coords[j + 1, 2])
      }
    }
    row.names(df) = NULL
    df$fx = as.numeric(df$fx)
    df$fy = as.numeric(df$fy)
    df$tx = as.numeric(df$tx)
    df$ty = as.numeric(df$ty)
    return(df)
  }
  
  
  require(sp)
  linedf2SpatialLines = function(linedf) {
    sl = list()
    for (i in 1:nrow(linedf)) {
      c1 = cbind(rbind(linedf$fx[i], linedf$tx[i]),
                 rbind(linedf$fy[i], linedf$ty[i]))
      l1 = Line(c1)
      sl[[i]] = Lines(list(l1), ID = linedf$id[i])
    }
    SL = SpatialLines(sl)
    return(SL)
  }
  
  
  #### Split the lines ####
  # Convert the input SpatialLine object into a dataframe and create an empty output dataframe
  linedf = SpatialLines2df(spatial_line)
  df = data.frame(
    id = character(),
    fx = numeric(),
    fy = numeric(),
    tx = numeric(),
    ty = numeric(),
    stringsAsFactors = FALSE
  )
  
  
  profvis::pause(1)
  
  for (i in 1:nrow(linedf)) {
    
    
    print(paste("Discretizing the road network ", round(100*i/nrow(linedf),1), "% completed...") )
    # For each line of the dataframe, corresponding to a single line of the spatial object
    # skips if length is less then split_length
    v_seg = linedf[i, ]
    seg_length = sqrt((v_seg$fx - v_seg$tx) ^ 2 + (v_seg$fy - v_seg$ty) ^
                        2) # segment length
    if (seg_length <= split_length) {
      df[nrow(df) + 1,] = c(paste0(v_seg$id, '_', '0000'),
                            v_seg$fx,
                            v_seg$fy,
                            v_seg$tx,
                            v_seg$ty)
      next()
    }
    # Create a vector of direction as the line and unit length
    # vector v corresponding to the line
    v = c(v_seg$tx  -  v_seg$fx,
          v_seg$ty  -  v_seg$fy)
    # vector of direction v and unit length
    u = c(v[1]  /  sqrt(v[1]  ^  2 + v[2]  ^  2), v[2]  /  sqrt(v[1]  ^  2 + v[2]  ^ 2))
    # Calculates how many segment the line is split into and the leftover
    num_seg = floor(seg_length  /  split_length)
    seg_left = seg_length - (num_seg  *  split_length)
    
    # Add to the output dataframe each segment plus the leftover
    for (i in 0:(num_seg  -  1)) {
      # Add num_seg segments
      df[nrow(df)  +  1,] = c(
        paste0(v_seg$id, '_', formatC(i, width = 4, flag = '0')),
        v_seg$fx + u[1]  *  split_length  *  i,
        v_seg$fy + u[2]  *  split_length  *  i,
        v_seg$fx + u[1]  *  split_length  *  (i  +  1),
        v_seg$fy + u[2]  *  split_length  *  (i  +  1)
      )
    }
    df[nrow(df) + 1,] = c(
      paste0(v_seg$id, '_', formatC(
        num_seg, width = 4, flag = '0'
      )),
      # Add leftover segment
      v_seg$fx + u[1] * split_length * num_seg,
      v_seg$fy + u[2] * split_length * num_seg,
      v_seg$tx,
      v_seg$ty
    )
    
  }
  
  print("Discretizing the road network shape files, Done!")
  
  
  #### Visualise the results to check ####
  if (plot.results) {
    plot(spatial_line)
    coords = cbind(as.numeric(df$fx), as.numeric(df$fy))
    coords = rbind(coords, as.numeric(df$tx[nrow(df)]), as.numeric(df$ty)[nrow(df)])
    sp_points = SpatialPoints(coords)
    plot(sp_points, col = 'red', add = T)
  }
  
  #### Output ####
  df$fx = as.numeric(df$fx)
  df$fy = as.numeric(df$fy)
  df$tx = as.numeric(df$tx)
  df$ty = as.numeric(df$ty)
  if (return.dataframe) {
    return(df)
  } # Return a dataframe
  sl = linedf2SpatialLines(df)
  return(sl) # Return a SpatialLine object
}
