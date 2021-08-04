compute_K = function(kernel_type, kw, D, W){
  if (kernel_type==1){
    K = W * exp(-0.5 * (D/kw)^2)
  }else if (kernel_type==2){
    K = W * exp(-D/kw)
  }else if (kernel_type==3){
    u = D/kw
    K = matrix(0,dim(D)[1],dim(D)[2])
    K[u<=1] = 0.75*(1-u[u<=1]^2)
    K = W * K
  }
  return(K)
}