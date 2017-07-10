#--------------------------------------------------------------------------
# HIERARCHICAL RISK PARITY
# Replication of Thomas Wiecki's Python code
# Author: Johann Pfitzinger

# FUNCTION LIBRARY

#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# Hierarchical Risk Parity (by LdP)
#--------------------------------------------------------------------------

      getClusterVar <- function(cov, cItems) {
        # Compute variance per cluster
        cov <- cov[cItems,cItems]
        if (length(cItems)>1) {
          ivp <- 1 / diag(cov)
          ivp <- ivp / sum(ivp)
          w <- ivp
          cVar <- crossprod(t(t(w)%*%cov),w)
        } else {
          cVar <- cov
        }
        return(cVar)
      }
      getRecBipart <- function(cov, sortIx) {
        w <- matrix(1, length(sortIx),1)
        w <- data.frame(w, row.names = sortIx)
        cItems <- list(sortIx)
        while (length(cItems) > 0) {
          cItemsNew <- vector("list",length(cItems)*2)
          for (i in 1:length(cItems)) {
            if (length(cItems[[i]])>1) {
              items <- cItems[[i]]
              cItems1 <- items[1:floor(length(items)/2)]
              cItems2 <- items[ceiling((length(items)+1)/2):length(items)]
              cVar1 <- getClusterVar(cov, cItems1)
              cVar2 <- getClusterVar(cov, cItems2)
              alpha <- 1 - cVar1 / (cVar1 + cVar2)
              w[cItems1,] <- w[cItems1,] * as.numeric(alpha)
              w[cItems2,] <- w[cItems2,] * as.numeric(1-alpha)
              cItemsNew[[i*2-1]] <- cItems1
              cItemsNew[[i*2]] <- cItems2
            }
          }
          cItemsNew <- cItemsNew[as.logical(1-unlist(lapply(cItemsNew,is.null)))]
          cItems <- cItemsNew
        }
        return(w)
      }
      correlDist <- function(corr) {
        # A distance matrix based on correlation, where 0<=d[i,j]<=1
        # This is a proper distance metric
        dist = ((1 - corr) / 2)^0.5  # distance matrix
        return(dist)
      }
      
      # Generates matrix of weights:
      getHRP <- function(cov, corr, return_raw = NULL, robust_cov = F) {
        # Construct a hierarchical portfolio
        if (robust_cov == T) {
          cov <- cov_shrink(return_raw)
          corr <- cov2cor(cov)
        }
        dist <- correlDist(corr)
        link <- hclust(dist(dist), method = "single")
        sortIx <- names(corr[,1])[link$order]
        hrp <- getRecBipart(cov, sortIx)
        return(hrp[colnames(rets),])
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
        ivp <- 1 / diag(cov)
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
      
     
      
      
      