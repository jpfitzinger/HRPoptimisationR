#--------------------------------------------------------------------------
# HIERARCHICAL RISK PARITY
# Replication of Thomas Wiecki's Python code
# Author: Johann Pfitzinger

# FUNCTION LIBRARY

#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Hierarchical Risk Parity (by LdP)
#--------------------------------------------------------------------------

      getHRP <- function(cov, corr, max = NULL, min = NULL, return_raw = NULL, robust_cov = F) {
        # Construct a hierarchical portfolio
        if (robust_cov == T) {
          cov <- cov_shrink(return_raw)
          corr <- cov2cor(cov)
        }
        
        # Set the constraint matrix
        if (is.null(max)) max <- rep(1,ncol(cov))
        if (is.null(min)) min <- rep(0,ncol(cov))
        if (length(max)==1) max <- rep(max,ncol(cov)) else if (length(max)<ncol(cov)) stop("Provide correct weights")
        if (length(min)==1) min <- rep(min,ncol(cov)) else if (length(min)<ncol(cov)) stop("Provide correct weights")
        const <- rbind(max, min)
        
        # check constraints
        if (sum(const[1,]) < 1 | sum(const[2,]) > 1) stop("Incompatible weights")
        
        distmat <- ((1 - corr) / 2)^0.5
        clustOrder <- hclust(dist(distmat), method = 'single')$order
        out <- getRecBipart(cov, clustOrder, const)
        return(out)
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
          alpha <- min(sum(const[1,cItems0]) / w[cItems0[1]],
                       max(sum(const[2,cItems0]) / w[cItems0[1]], 
                           alpha))
          alpha <- 1 - min(sum(const[1,cItems1]) / w[cItems1[1]], 
                           max(sum(const[2,cItems1]) / w[cItems1[1]], 
                               1 - alpha))
          
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
      
#--------------------------------------------------------------------------
# Mean-Variance Efficient Portfolio
#--------------------------------------------------------------------------
      
      # Using Quadratic programming, no shorting, no leverage
      # shrink_means = TRUE produces the minimum variance portfolio
      get_mean_variance <- function(cov, return_raw, shrink_means = FALSE, robust_cov = FALSE, max.weight = 1) {
        
        Amat <- cbind(matrix(1,ncol(rets),1), diag(ncol(rets)), -diag(ncol(rets)))
        bvec <- c(1,rep(0,ncol(rets)), rep(-max.weight,ncol(rets)))
        return <- colMeans(return_raw)
        
        if (robust_cov == T) {
          cov <- cov_shrink(return_raw)
        }
        
        if (shrink_means == F) {
      
          N <- 50
          w <- matrix(NA,N,ncol(rets))
          qs = 10^(5 * c(1:N)/N - 1)
          for (i in seq_along(qs)) {
            q <- qs[i]
            w[i,] <- solve.QP(Dmat = q*as.matrix(cov), dvec = return, Amat = Amat, bvec = bvec, meq = 1)$solution
          }
          params <- matrix(NA,N,2)
          params[,1] <- return %*% t(w)
          params[,2] <- sqrt(diag(w %*% cov %*% t(w)))
        
          poly.fit <- lm(params[,2] ~ cbind(params[,1],params[,1]^2))
          poly.coef <- coefficients(poly.fit)
          q.opt <- sqrt(poly.coef[1]/poly.coef[3])
        
          if (!is.nan(q.opt) & !is.na(q.opt) & is.finite(q.opt)) {
            out <- solve.QP(Dmat = q.opt*as.matrix(cov), dvec = return, Amat = Amat, bvec = bvec, meq = 1)$solution
          } else {out <- rep(NA,ncol(rets))}
        } else {
          out <- solve.QP(Dmat = 2*as.matrix(cov), dvec = 0*return, Amat = Amat, bvec = bvec, meq = 1)$solution
        }
                      
        return(out)
      }
      
#--------------------------------------------------------------------------
# Inverse Variance Portfolio
#--------------------------------------------------------------------------
      
      getIVP <- function(cov, return_raw = NULL, robust_cov = F) {
        if (robust_cov == T) {
          cov <- cov_shrink(return_raw)
        }
        ivp <- 1 / diag(as.matrix(cov))
        ivp <- ivp / sum(ivp)
        return(ivp)
      }
      
#--------------------------------------------------------------------------
# Method Library
#--------------------------------------------------------------------------
      
      opt_methods <- list(MVO = c("Mean-Variance weighting","get_mean_variance(cov, return_raw)"),
                          rMVO = c("Robust Mean-Varniance weighting","get_mean_variance(cov, return_raw,robust_cov = T)"),
                          cMVO = c("Constrained Mean-Variance weighting", "get_mean_variance(cov, return_raw, max.weight = max.w)"),
                          minVar = c("Min-Variance weighting","get_mean_variance(cov, return_raw,T)"),
                          rminVar = c("Robust Min-Variance weighting","get_mean_variance(cov, return_raw,T,robust_cov = T)"),
                          invVar = c("Inverse Variance weighting","getIVP(cov)"),
                          rinvVar = c("Robust Inverse Variance weighting","getIVP(cov,return_raw,robust_cov = T)"),
                          Equal = c("Equal weighting","rep(1/ncol(rets),ncol(rets))"),
                          HRP = c("Hierarchical weighting","getHRP(cov, corr)"),
                          cHRP = c("Constrained Hierarchical weighting","getHRP(cov, corr, max = max.w)"),
                          rHRP = c("Robust Hierarchical weighting","getHRP(cov, corr, return_raw, robust_cov = T)")
                         )
      
#--------------------------------------------------------------------------
# Output Function
#--------------------------------------------------------------------------

      optFx <- function(x, calc_int, rebal_per, 
                        methods = c("MVO","rMVO","minVar","rminVar","invVar","rinvVar","Equal","HRP","rHRP")) {
        
        # Set up arguments
        nopt <- length(methods)
        dates <- index(rets)
        eops <- dates[calc_int:length(dates)][endpoints(dates[calc_int:length(dates)], rebal_per)] %>% as.character
        
        # Create empty panels to store rebalancing values
        covs <- array(dim = c(ncol(rets),ncol(rets),length(eops)), 
                      dimnames = list(colnames(rets),colnames(rets),eops))
        corrs <- array(dim = c(ncol(rets),ncol(rets),length(eops)), 
                      dimnames = list(colnames(rets),colnames(rets),eops))
        returns <- array(dim = c(ncol(rets), length(eops)),
                         dimnames = list(colnames(rets),eops))
        w <- array(dim = c(nopt,ncol(rets),length(eops)), 
                   dimnames = list(methods,colnames(rets),eops))
        weights <- array(dim = c(nopt,ncol(rets),length(dates)), 
                         dimnames = list(methods,colnames(rets),as.character(dates)))
        
        # Rebalancing loop
        for (eop in eops) {
          date.series <- seq.Date(max(as.Date(eop) - calc_int, dates[1]), as.Date(eop), "day")
          cov <- rets[date.series,] %>% cov
          corr <- rets[date.series,] %>% cor
          return_raw <- rets[date.series,]
          return <- colMeans(return_raw)
          covs[,,eop] <- cov
          corrs[,,eop] <- corr
          returns[,eop] <- return
          
          # Fetch optimisation from method library
          for (i in seq_along(methods)) {
            fx <- opt_methods[[methods[i]]][2]
            w[i,,eop] <- eval(parse(text = fx))
          }
        }
        
        # Calculate portfolio returns
        weights[,,eops] <- w
        port_returns <- matrix(NA,length(dates),nopt) %>% data.frame(row.names = dates)
        for (i in 1:nopt) {
          weights[i,,] <- t(na.locf(t(weights[i,,]), na.rm = F))
          port_returns[,i] <- colSums(t(as.matrix(rets)) * weights[i,,])
        }
        port_returns[is.na(port_returns)] <- 0
        colnames(port_returns) <- methods
        
        output <- list(port_returns = port_returns, weights = w)
        return(output)
      }
      
     
      
      
      