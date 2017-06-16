#--------------------------------------------------------------------------
# HIERARCHICAL RISK PARITY
# Replication of Thomas Wiecki's Python code
# Author: Johann Pfitzinger

# Notes:
# The below provides an exact replication of HRP, Mean-Variance, Minimum Variance,
# Inverse Variance, and Equal weighting optimisation methods.
# No replication has been provided for CLA minimum variance optimisation.
# Package "tawny" used to obtain shrinkage estimates of the covariance matrix.
# This is distinct from the OAS method used in the source file, of which no 
# direct implementation exists in R. Other options can be explored in future.
#--------------------------------------------------------------------------

      pacman::p_load(dplyr, zoo, reshape2, ggplot2, corrplot, quantmod, tawny)

#--------------------------------------------------------------------------
# Import Data
#--------------------------------------------------------------------------

      # Fetch some JSE shares (customizable). Currently using Google Finance.
      symbols <- c("JSE:AMS","JSE:CFR","JSE:FSR","JSE:KIO","JSE:MTN","JSE:NPN","JSE:SAB","JSE:SBK","JSE:SOL","JSE:VOD")
      getSymbols(symbols, src="google", return.class = "zoo")
      # Merge closing prices
      prices.data <- get(symbols[1])[,4]
      for (i in 2:length(symbols)) {
        prices.data <- merge.zoo(prices.data, get(symbols[i])[,4], all = F)
      }
      colnames(prices.data) <- symbols
      rm(list = symbols)
      # Calculate Return
      rets <- diff(log(prices.data))
      # Export for comparison with Python
      write.csv(as.data.frame(rets), "./Data/rets.csv")
      #python.load(HRP.py)

      # Fetch functions library
      source("functions.R")

#--------------------------------------------------------------------------
# Optimisation Function
      
# This function generates optimal portfolio weights. The function take the following args:
# - x                    = zoo object of returns
# - calc_int             = period length of historicals used to calculate risk parameters
# - rebal_per            = frequency of rebalancing, c("days","months","quarters","years")
# - methods              = optimization method:
#   - MVO / rMVO         = Mean-Variance optimisation (w/ & w/out robust covariance)
#   - invVar / rinvVar   = Inverse Variance
#   - minVar / rminVar   = Minimum Variance
#   - Equal              = Equal weighting
#   - HRP / rHRP         = Hierarchical Risk Parity
#   Robust covariance matrices are calculated using the "tawny" package (not OAS used by Wiecki)
#   Multiple methods can be passed to obtain comparisons.
#--------------------------------------------------------------------------

      p <- optFx(rets, calc_int = 252, rebal_per = "months", methods = c("HRP","MVO"))
  
#--------------------------------------------------------------------------
# Plots to replicate Thomas Wiecki (Python code)
#--------------------------------------------------------------------------

      # Correlation Plot & Clustering
      
      cor(prices.data) %>% corrplot(method = "ellipse")
      cor(prices.data) %>% correlDist %>% dist %>% hclust("single") %>% plot
      
      # Cumulative log returns
      
      l.rets <- p$port_returns %>% log1p %>% cumsum %>% melt
      l.rets$date <- rownames(p$port_returns) %>% as.Date %>% rep(ncol(p$port_returns))
      l.rets$label <- sapply(as.character(l.rets$variable), function(x) opt_methods[[x]][1])
      
      ggplot(data = l.rets, aes(x = date)) + 
        geom_line(aes(y = value, color = label)) +
        scale_color_brewer(palette = "Set1", name = NULL) +
        labs(x = NULL, y = "Cumulative Log Returns") + 
        scale_y_continuous(breaks = seq(-2,2,0.25)) +
        theme_bw()
        
      # Sharpe Ratio & Volatility
      
      sr <- colMeans(p$port_returns)*252 / (apply(p$port_returns,2,sd)*sqrt(252))
      ggplot() + geom_col(aes(sapply(as.character(names(sr)), function(x) opt_methods[[x]][1]), sr, fill = names(sr))) +
        labs(x = NULL, y = "Sharpe Ratio") +
        scale_fill_brewer(palette = "Set1") +
        theme_bw() + theme(legend.position = "none")
      
      vol <- apply(p$port_returns,2,sd)*sqrt(252)
      ggplot() + geom_col(aes(sapply(as.character(names(vol)), function(x) opt_methods[[x]][1]), vol, fill = names(vol))) +
        labs(x = NULL, y = "Portfolio Volatility") +
        scale_fill_brewer(palette = "Set1") +
        theme_bw() + theme(legend.position = "none")
