apply_weights <- function(w, idx, max = NULL, min = NULL, maxit, accuracy) {
  
  w <- rbind(1,w)
  idx <- rbind(1,idx)
  if(is.null(max)) max <- rep(1, ncol(w))
  if(length(max)==1) max <- rep(max,ncol(w))
  if(is.null(min)) min <- rep(0, ncol(w))
  if(length(min)==1) min <- rep(min, ncol(w))
  
  iter <- 0
  conv <- 1
  while (conv > 0 & iter < maxit) {
    
    cw <- apply(w,2,cumprod)
    const <- rbind(max, min)
    const_nxt <- const
    
    for (i in nrow(w):2) {
      for (j in 1:ncol(w)) {
        if (!is.na(w[i,j])) {
          
          idx0 <- idx[i,j]==idx[i,]
          idx0[is.na(idx0)] <- FALSE
          idx1 <- idx[i-1,j]==idx[i-1,]
          idx1[idx0==1] <- FALSE
          idx1[is.na(idx1)] <- FALSE
          
          min_val <- const[2,j] + ifelse(sum(idx1)>0,max(const[2,idx1]),0)
          max_val <- const[1,j] + ifelse(sum(idx1)>0,min(const[1,idx1]),0)
          
          if (cw[i-1,j] < min_val) {
            const_nxt[c(1:2),idx0|idx1] <- min_val
          }
          if (cw[i-1,j] > max_val) {
            const_nxt[c(1:2),idx0|idx1] <- max_val
          }
          if (cw[i-1,j] <= max_val & cw[i-1,j] >= min_val) {
            const_nxt[c(1:2),idx0|idx1] <- cw[i-1,j]
          }
          if (w[i,j]*const_nxt[2,j] < const[2,j]) {
            w[i,idx0] <- const[2,j] / const_nxt[2,j]
            w[i,idx1] <- 1 - w[i,j]
          }
          if (w[i,j]*const_nxt[1,j] > const[1,j]) {
            w[i,idx0] <- const[1,j] / const_nxt[1,j]
            w[i,idx1] <- 1 - w[i,j]
          }
        }
      }
      const <- const_nxt
    }
    
    w_vec <- apply(w, 2, prod, na.rm = T)
    error <- c(w_vec - max, min - w_vec)
    error <- error[error>0]
    conv <- round(sum(error),accuracy)
    iter <- iter + 1
  }
  return(w_vec)
  
}
