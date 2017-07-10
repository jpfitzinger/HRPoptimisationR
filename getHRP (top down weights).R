      getHRP <- function(cov, corr, max = NULL, min = NULL, return_raw = NULL, robust_cov = F) {
        # Construct a hierarchical portfolio
        if (robust_cov == T) {
          cov <- cov_shrink(return_raw)
          corr <- cov2cor(cov)
        }
        
        # Set the constraint matrix
        if (is.null(max)) max <- rep(1,ncol(cov))
        if (is.null(min)) min <- rep(0,ncol(cov))
        if (length(max)==1) max <- rep(max,ncol(cov))
        if (length(min)==1) min <- rep(min,ncol(cov))
        const <- rbind(max, min)
        
        distmat <- ((1 - corr) / 2)^0.5
        clustOrder <- hclust(dist(distmat), method = 'single')$order
        out <- getRecBipart(cov, clustOrder, const)
        return(out)
      }
      
      getIVP <- function(cov) {
        # get inverse variance portfolio from diagonal of covariance matrix
        invDiag <- 1/diag(as.matrix(cov))
        weights <- invDiag/sum(invDiag)
        return(weights)
      }
      
      getClusterVar <- function(cov, cItems) {
        # compute cluster variance from the inverse variance portfolio above
        covSlice <- cov[cItems, cItems]
        weights <- getIVP(covSlice)
        cVar <- t(weights) %*% as.matrix(covSlice) %*% weights
        return(cVar)
      }
      
      getRecBipart <- function(cov, sortIx, const) {
      
        w <- rep(1, ncol(cov))
        
        # create recursion function within parent function to avoid use of globalenv
        recurFun <- function(cov, sortIx, const) {
          # get first half of sortIx which is a cluster order
          subIdx <- 1:trunc(length(sortIx)/2)
          
          # subdivide ordering into first half and second half
          cItems0 <- sortIx[subIdx]
          cItems1 <- sortIx[-subIdx]
          
          # compute cluster variances of covariance matrices indexed
          # on first half and second half of ordering
          cVar0 <- getClusterVar(cov, cItems0)
          cVar1 <- getClusterVar(cov, cItems1)
          alpha <- 1 - cVar0/(cVar0 + cVar1)
          
          # determining whether weight constraint binds
          w0 <- w[cItems0[1]] * alpha
          w1 <- w[cItems1[1]] * (1-alpha)
          if (w[cItems0[1]] < sum(const[2,c(cItems0,cItems1)])-0.00001 | w[cItems0[1]] > 0.00001+sum(const[1,c(cItems0,cItems1)])) {
            stop("Incompatible weights")
          }
          if (w0 < sum(const[2,cItems0])) {
            alpha0 <- sum(const[2,cItems0]) / w[cItems0[1]]
            alpha <- alpha0
          }
          if (w0 > sum(const[1,cItems0])) {
            alpha0 <- sum(const[1,cItems0]) / w[cItems0[1]]
            alpha <- alpha0
          }
          if (w1 < sum(const[2,cItems1])) {
            alpha1 <- sum(const[2,cItems1]) / w[cItems1[1]]
            alpha <- 1 - alpha1
          }
          if (w1 > sum(const[1,cItems1])) {
            alpha1 <- sum(const[1,cItems1]) / w[cItems1[1]]
            alpha <- 1 - alpha1
          }
          
          w[cItems0] <<- w[cItems0] * rep(alpha, length(cItems0))
          w[cItems1] <<- w[cItems1] * rep((1-alpha), length(cItems1))
          
          # rerun the function on a half if the length of that half is greater than 1
          if(length(cItems0) > 1) {
            recurFun(cov, cItems0, const)
          }
          if(length(cItems1) > 1) {
            recurFun(cov, cItems1, const)
          }
          
        }
        
        # run recursion function
        recurFun(cov, sortIx, const)
        return(w)
      }