      getHRP <- function(cov, corr, return_raw = NULL, robust_cov = F, max = NULL, min = NULL, accuracy = 10, maxit = 500) {
        # Construct a hierarchical portfolio
        if (robust_cov == T) {
          cov <- cov_shrink(return_raw)
          corr <- cov2cor(cov)
        }
        distmat <- ((1 - corr) / 2)^0.5
        clustOrder <- hclust(dist(distmat), method = 'single')$order
        out <- getRecBipart(cov, clustOrder)
        out <- apply_weights(w = out, idx, max, min, maxit, accuracy)
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
      
      getRecBipart <- function(cov, sortIx) {
        # keeping track of weights vector in the global environment
        assign("w", value = matrix(NA, ceiling(log(ncol(cov))/log(2)), ncol(cov)), envir = .GlobalEnv)
        assign("idx", value = matrix(NA, ceiling(log(ncol(cov))/log(2)), ncol(cov)), envir = .GlobalEnv)
        
        # run recursion function
        recurFun(cov, sortIx)
        return(w)
      }
      
      recurFun <- function(cov, sortIx) {
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
        
        # updating weights outside the function using scoping mechanics 
        w[min(which(is.na(w[,cItems0]))),cItems0] <<- alpha
        w[min(which(is.na(w[,cItems1]))),cItems1] <<- (1-alpha)
        idx[min(which(is.na(idx[,cItems0]))),cItems0] <<- cItems0[1]
        idx[min(which(is.na(idx[,cItems1]))),cItems1] <<- cItems1[1]
        
        
        # rerun the function on a half if the length of that half is greater than 1
        if(length(cItems0) > 1) {
          recurFun(cov, cItems0)
        }
        if(length(cItems1) > 1) {
          recurFun(cov, cItems1)
        }
      }
      