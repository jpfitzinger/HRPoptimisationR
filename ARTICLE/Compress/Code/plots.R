library(rugarch)
library(dplyr)
library(latex2exp)
library(ggplot2)
library(reshape2)

mat.main <- ggplot(cval2, aes(variable2.f, response)) + theme_classic(base_size = 8, base_family = "Times") +
  geom_tile(aes(fill = value), colour = "gray90", size = 2, stat = "identity") + 
  geom_text(aes(label = text), size = 2) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(name = NULL, low = "darkorange3", mid = "white", midpoint = 0, 
                       high = "dodgerblue", space = "Lab", na.value = "gray90", guide = F) +
  labs(x = NULL, y = NULL, caption = "Colour scale visualizes magnitude of spillover, where insignificant values are greyed out. Columns represent source variables; Rows represent response variables. Significance Level: *** 1%, ** 5%, * 10%.") +
  facet_grid(grouping.y.f~grouping.f, scales = "free", space = "free", switch = "y")


# CB

cb.plot.1 <- ggplot() + theme_classic(base_family = "Times", base_size = 9) +
  geom_line(aes(x = seq.Date(from = as.Date("1990/01/1"), to = as.Date("2017/04/30"), by = "month"), y = log(imat["BOND_FL","FED..BS.",]*100), color = "FED -> BOND FLOW")) +
  geom_line(aes(x = seq.Date(from = as.Date("1990/01/1"), to = as.Date("2017/04/30"), by = "month"), y = log(imat["BOND_FL","ECB..BS.",]*100), color = "ECB -> BOND FLOW")) +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set1", name = NULL) +
  scale_x_date(limits = c(as.Date("2000/01/1"), NA)) +
  #ylim(c(0,0.02)) +
  labs(y = TeX("$log I_t$"),x=NULL)
cb.plot.2 <- ggplot() + theme_classic(base_family = "Times", base_size = 9) +
  geom_line(aes(x = seq.Date(from = as.Date("1990/01/1"), to = as.Date("2017/04/30"), by = "month"), y = log(imat["EQY_FL","FED..BS.",]*100), color = "FED -> EQY FLOW")) +
  geom_line(aes(x = seq.Date(from = as.Date("1990/01/1"), to = as.Date("2017/04/30"), by = "month"), y = log(imat["EQY_FL","ECB..BS.",]*100), color = "ECB -> EQY FLOW")) +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set1", name = NULL) +
  scale_x_date(limits = c(as.Date("2000/01/1"), NA)) +
  #ylim(c(0,0.02)) +
  labs(y = TeX("$log I_t$"),x=NULL)

cb.plot.3 <- ggplot() + theme_classic(base_family = "Times", base_size = 9) +
  geom_line(aes(x = seq.Date(from = as.Date("1990/01/1"), to = as.Date("2017/04/30"), by = "month"), y = log(imat["SA..EQY.","FED..BS.",]*100), color = "FED -> EQY")) +
  geom_line(aes(x = seq.Date(from = as.Date("1990/01/1"), to = as.Date("2017/04/30"), by = "month"), y = log(imat["SA..EQY.","ECB..BS.",]*100), color = "ECB -> EQY")) +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set1", name = NULL) +
  scale_x_date(limits = c(as.Date("2000/01/1"), NA)) +
  #ylim(c(0,0.02)) +
  labs(y = TeX("$log I_t$"),x=NULL)
cb.plot.4 <- ggplot() + theme_classic(base_family = "Times", base_size = 9) +
  geom_line(aes(x = seq.Date(from = as.Date("1990/01/1"), to = as.Date("2017/04/30"), by = "month"), y = log(imat["USD..FX.","FED..BS.",]*100), color = "FED -> USD:ZAR")) +
  geom_line(aes(x = seq.Date(from = as.Date("1990/01/1"), to = as.Date("2017/04/30"), by = "month"), y = log(imat["EUR..FX.","ECB..BS.",]*100), color = "ECB -> EUR:ZAR")) +
  theme(legend.position = "top") +
  scale_color_brewer(palette = "Set1", name = NULL) +
  scale_x_date(limits = c(as.Date("2000/01/1"), NA)) +
  #ylim(c(0,0.02)) +
  labs(y = TeX("$log I_t$"),x=NULL)
