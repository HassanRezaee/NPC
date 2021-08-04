



setwd("C:/Users/Koala/Documents/1_Academic/PostDoc/HEC_McGill/Publication/Methodology_paper/Codes/NPC_code_publication")




remove(list=ls())
while (!is.null(dev.list()))  dev.off()




# load pre-defined functions

{
  source("SplitLines.R")
  source("compute_K.R")
  source("compute_D.R")
  source("make_graph.R")
}





# load data

load("intersection_crash_data_sim.RData")

n = length(y)


# Create graph and compute distances --------------------------------------

{

  # Directize the graph to scatter the knots
  graph_data = make_graph(roads_shp, delta_l = 50); # delta_l is the disretization unit length
  igr = graph_data[[1]]
  nodes_coords = graph_data[[2]]

  degrees = as.matrix(igraph::degree(igr, v = igraph::V(igr), mode = "total", loops = F, normalized = FALSE))
  
  
  # Locate data on the discretized graph
  knns = FNN::get.knnx(query=coords, data=nodes_coords, k=5, algo="kd_tree")
  data_nodes = as.matrix(knns[["nn.index"]][,1])

    
 
  
  # define knot locations
  {
    library(plot3D)
    mx = n*3
    xg = seq(from=min(nodes_coords[,1]), to=max(nodes_coords[,1]), length.out=round(sqrt(mx)))
    yg = seq(from=min(nodes_coords[,2]), to=max(nodes_coords[,2]), length.out=round(sqrt(mx)))
    xy = mesh(xg,yg);xy=cbind(c(xy$x),c(xy$y))
    
    
    m = nrow(xy)
    
    
    knns = FNN::get.knnx(query=xy, data=nodes_coords, k=5, algo="kd_tree")
    knot_ids = as.matrix(knns[["nn.index"]][,1])
  
    
    # knots that are within the convex hull of data locations
    pl = as.matrix(nodes_coords[data_nodes,])
    point_in_ids = sp::point.in.polygon(nodes_coords[knot_ids,1], nodes_coords[knot_ids,2], pl[,1], pl[,2])
    point_in_ids = which(point_in_ids==1);
    knot_ids = knot_ids[point_in_ids]; m=length(knot_ids)
    
    
    
    # display data and knots on the road network
    if (1==1){
      w=3; par(mfrow=c(1,1), mar = c(w, w, w, w));
      plot(roads_shp,axes=T)
      points(coords,col="red",pch=19,cex=0.75); 
      points(nodes_coords[knot_ids,],pch=19,col="blue",cex=1);
      points(nodes_coords[data_nodes,],col="green",cex=1);
    }
    
  }
  
  
  
  
  
}


# Compute path distances --------------------------------------------------

DW = compute_D(igr, data_nodes, knot_ids, nodes_coords); D = DW[[1]]; W = DW[[2]]



# Compute weighted densities ----------------------------------------------

kernel_width = 1500
K = compute_K(3, kernel_width, D, W)



# prepare input data for INLA ---------------------------------------------

inla_data = data.frame("id" = 1:n, "crash" = y, "covar1" = covars[,1], "covar2" = covars[,2], "covar3" = covars[,3])






# Run PR (Poisson regression) ---------------------------------------------

{
  formula_pr <- crash ~ 1 + covar1 + covar2 + covar3
  inla_fit_pr <- INLA::inla(formula_pr, family = "poisson", data = inla_data, control.compute=list(config=T,dic=T,waic=T),
                      control.inla = list(h = 1),
                      verbose=T)
  
  y_hat_pr = inla_fit_pr$summary.fitted.values$mean
  
  y_hat_pr_qts = matrix(NA, n, 2)
  y_hat_pr_qts[,1] = inla_fit_pr[["summary.fitted.values"]][["0.025quant"]]
  y_hat_pr_qts[,2] = inla_fit_pr[["summary.fitted.values"]][["0.975quant"]]
  
  
  beta_hat_pr = inla_fit_pr[["summary.fixed"]][["mean"]]
  beta_hat_pr_qts1 = inla_fit_pr[["summary.fixed"]][["0.025quant"]]
  beta_hat_pr_qts2 = inla_fit_pr[["summary.fixed"]][["0.975quant"]]
}



# Run pCAR ----------------------------------------------------------------
{
  g = INLA::inla.read.graph(adj_mat)
  formula_car <- crash ~ 1 + covar1 + covar2 + covar3 + f(id, model = "besagproper", graph = g, constr= FALSE,
                                                          hyper = list(prec = list(
                                                            prior = "loggamma",
                                                            initial = 0,
                                                            fixed = F,
                                                            param = c(1.316611, 0.3361755))))
  
  inla_fit_car <- INLA::inla(formula_car, family = "poisson", data = inla_data, control.compute=list(config=T,dic=T,waic=T),
                       control.inla = list(h = 1),
                       verbose=T, control.fixed = list(prec.intercept = 0.3, prec=1))
  
  z_hat_car = as.matrix(inla_fit_car[["summary.random"]][["id"]][["mean"]])
  
  
  y_hat_car = as.matrix(inla_fit_car$summary.fitted.values$mean)
  
  beta_hat_car = inla_fit_car[["summary.fixed"]][["mean"]]
  beta_hat_car_qts1 = inla_fit_car[["summary.fixed"]][["0.025quant"]]
  beta_hat_car_qts2 = inla_fit_car[["summary.fixed"]][["0.975quant"]]
  
  z_hat_car_qts = matrix(NA, n, 2)
  z_hat_car_qts[,1] = inla_fit_car[["summary.random"]][["id"]][["0.025quant"]]
  z_hat_car_qts[,2] = inla_fit_car[["summary.random"]][["id"]][["0.975quant"]]
  
  y_hat_car_qts = matrix(NA, n, 2)
  y_hat_car_qts[,1] = inla_fit_car[["summary.fitted.values"]][["0.025quant"]]
  y_hat_car_qts[,2] = inla_fit_car[["summary.fitted.values"]][["0.975quant"]]
  
  print("CAR done")
}




# Run NPC -----------------------------------------------------------------

{
  
  formula_npc <- crash ~ 1 + covar1 + covar2 + covar3 + f(id, model = "z", Z = K,
                                                          hyper = list(prec = list(
                                                            prior = "loggamma",
                                                            initial = 0,
                                                            fixed = F,
                                                            param = c(1.316611, 0.3361755))));
  
  
  inla_fit_npc <- INLA::inla(formula_npc, family = "poisson", data = inla_data, control.compute=list(config=T,dic=T,waic=T),
                       control.inla = list(h=1),
                       verbose=T, control.fixed = list(prec.intercept = 0.3, prec=1))
  
  
  y_hat_npc = inla_fit_npc[["summary.fitted.values"]][["0.5quant"]]
  
  x_hat = inla_fit_npc[["summary.random"]][["id"]][["mean"]][-(1:n)]
  z_hat_npc = K%*%x_hat
  
  beta_hat_npc = inla_fit_npc[["summary.fixed"]][["mean"]]
  beta_hat_npc_qts1 = inla_fit_npc[["summary.fixed"]][["0.025quant"]]
  beta_hat_npc_qts2 = inla_fit_npc[["summary.fixed"]][["0.975quant"]]
  
  x_hat_npc_qts = matrix(NA, m, 2)
  x_hat_npc_qts[,1] = inla_fit_npc[["summary.random"]][["id"]][["0.025quant"]][-(1:n)]#[tr_ids]
  x_hat_npc_qts[,2] = inla_fit_npc[["summary.random"]][["id"]][["0.975quant"]][-(1:n)]#[tr_ids]
  
  
  z_hat_npc_qts = matrix(NA, n, 2)
  z_hat_npc_qts[,1] = inla_fit_npc[["summary.random"]][["id"]][["0.025quant"]][1:n]#[tr_ids]
  z_hat_npc_qts[,2] = inla_fit_npc[["summary.random"]][["id"]][["0.975quant"]][1:n]#[tr_ids]
  
  y_hat_npc_qts = matrix(NA, n, 2)
  y_hat_npc_qts[,1] = inla_fit_npc[["summary.fitted.values"]][["0.025quant"]]
  y_hat_npc_qts[,2] = inla_fit_npc[["summary.fitted.values"]][["0.975quant"]]
}





# display the scatter plot between true crash data and the fitted values

{
  w=4; par(mfrow=c(1,3), mar = c(w, w, w, w));
  m1 = c(0, max(y_hat_pr_qts[,2], y_hat_car_qts[,2], y_hat_npc_qts[,2]))
  
  
  plot(y, y_hat_pr, xlim = m1, ylim=m1, ylab=expression(hat(y)),xlab=expression(y), main="PR");abline(0,1)
  for (i in 1:n){
    segments(x0=y[i], y0=y_hat_pr_qts[i,1], x1 = y[i], y1 = y_hat_pr_qts[i,2], col = "gray", lwd = 1)
  }
  points(y, y_hat_pr, pch=19, cex=1, col="black")
  
  
  plot(y, y_hat_car, xlim = m1, ylim=m1, ylab=expression(hat(y)),xlab=expression(y), main="pCAR");abline(0,1)
  for (i in 1:n){
    segments(x0=y[i], y0=y_hat_car_qts[i,1], x1 = y[i], y1 = y_hat_car_qts[i,2], col = "gray", lwd = 1)
  }
  points(y, y_hat_car, pch=19, cex=1, col="black")
  
  
  
  plot(y, y_hat_npc, xlim = m1, ylim=m1, ylab=expression(hat(y)),xlab=expression(y), main="NPC");abline(0,1)
  for (i in 1:n){
    segments(x0=y[i], y0=y_hat_npc_qts[i,1], x1 = y[i], y1 = y_hat_npc_qts[i,2], col = "gray", lwd = 1)
  }
  points(y, y_hat_npc, pch=19, cex=1, col="black")
  
}





