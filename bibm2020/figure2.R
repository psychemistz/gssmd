## Figure 2
rm(list=ls())

## Source
source("C:/Users/sypark/Desktop/Projects/w_Shujaat/0ovltools-GSSMD/2code/src/gssmd_source.R")

## parameter
nol <- c(0, 10, 30, 50, 100, 200, 300)
savename <- paste0("outlier_lnorm", nol/10, "_res.csv")
sigma = 1
xlim = c(-15, 45)
ylim = c(0,0.13) # for sigma=1

for (k in 1:length(nol)){
  # parameter
  m_list <- c(0,1,3,5,10,20,30) 

  # Sample
  set.seed(111)
  target <- rlnorm(1000-nol[k], 0, 1)
  background <- rlnorm(1000, 0, 1)
  outlier <- list(rnorm(nol[k],m_list[1],1), rnorm(nol[k],m_list[2],1), rnorm(nol[k],m_list[3],1), 
                  rnorm(nol[k],m_list[4],1), rnorm(nol[k],m_list[5],1), rnorm(nol[k],m_list[6],1), rnorm(nol[k],m_list[7],1)) 
  
  # generate target matrix
  target_temp <- matrix(rep(0,7000), ncol=7)
  for (i in 1:7){
    if(length(outlier[[i]]) == 0) target_temp[,i] <- target
    else target_temp[,i] <- c(target, outlier[[i]])
  }
  target <- target_temp
  
  # qc measure
  zpfactor <- ssmd <- robust_ssmd <- gcnr <- gssmd <- c(); 
  
  for (i in 1:length(m_list)){
    zpfactor[i] <- qc.measure(target[,i], background, "zpfactor")
    ssmd[i] <- qc.measure(target[,i], background, "ssmd")
    #robust_ssmd[i] <- qc.measure(target[,i], background, "robust_ssmd")
    #gcnr[i] <- qc.measure(target[,i], background, "gcnr")
    gssmd[i] <- qc.measure(target[,i], background, "gssmd")
  }
  
  res <- round(data.frame(zpfactor, ssmd, gssmd),3)
  res <- cbind(outlier_mean=m_list, outlier_percentile=nol[k]/10, res)
  write.csv(res, paste0("C:/Users/sypark/Desktop/Projects/w_Shujaat/0ovltools-GSSMD/3results/", savename[[k]]), row.names = FALSE)
}

# read all results
outlier_res <- list.files("C:/Users/sypark/Desktop/Projects/w_Shujaat/0ovltools-GSSMD/3results/", full.names = TRUE)[grep("outlier_lnorm", list.files("C:/Users/sypark/Desktop/Projects/w_Shujaat/0ovltools-GSSMD/3results/"))]
res <- lapply(outlier_res, function(x) read.csv(x))
res <- do.call(rbind, res)

# print figure.2
library(tidyr)
library(dplyr)
library(ggplot2)
#res_total <- read.csv("C:/Users/sypark/Desktop/2018-Fall-KAIST-Research/3.GCNR/4.Results/GCNR_w_1_outlier/outlier_results.csv")
res_total <- res
res_gather <- res_total %>% gather(key = measure, value = value, zpfactor, ssmd, gssmd)

res_gather %>% group_by(outlier_mean) %>% 
  filter(measure != "zpfactor") %>%
  ggplot(aes(x=outlier_percentile, y=value, col = measure)) + 
  geom_point() + 
  facet_grid(.~outlier_mean) + 
  theme_bw()

# print table.2
library(knitr)
library(magrittr)
library(kableExtra)
res_mean_ordered <- res[order(res$outlier_mean),]

## mean_diff=3
res_3 <- res_mean_ordered[res_mean_ordered$outlier_mean==3,]
res_3 <- res_3[order(res_3$outlier_percentile),]
rownames(res_3) <- c("Diff_mu=0","Diff_mu=1", "Diff_mu=3", "Diff_mu=5", "Diff_mu=10","Diff_mu=20","Diff_mu=30")
kable(res_3[,-(1:2)], format = 'latex', booktabs = TRUE) %>% kable_styling(position = "center")

## mean_diff=5
res_5 <- res_mean_ordered[res_mean_ordered$outlier_mean==5,]
res_5 <- res_5[order(res_5$outlier_percentile),]
rownames(res_5) <- c("Diff_mu=0","Diff_mu=1", "Diff_mu=3", "Diff_mu=5", "Diff_mu=10","Diff_mu=20","Diff_mu=30")
kable(res_5[,-(1:2)], format = 'latex', booktabs = TRUE) %>% kable_styling(position = "center")

## mean_diff=30
res_30 <- res_mean_ordered[res_mean_ordered$outlier_mean==30,]
res_30 <- res_30[order(res_30$outlier_percentile),]
rownames(res_30) <- c("Diff_mu=0","Diff_mu=1", "Diff_mu=3", "Diff_mu=5", "Diff_mu=10","Diff_mu=20","Diff_mu=30")
kable(res_30[,-(1:2)], format = 'latex', booktabs = TRUE) %>% kable_styling(position = "center")

## Distribution Fitting
library(scales)
x <- seq(-40, 40, 0.05)
plot(x, dlnorm(x, 0,sigma), col = "red", lwd = 3, type = 'l', xlim = xlim, ylim = ylim, ylab = "", xlab = "")
lines(x, dlnorm(x, 0,sigma), col = "orange", lwd = 3)
lines(x, dlnorm(x, 1,sigma), col = "cyan", lwd = 3)
lines(x, dlnorm(x, 3,sigma), col = "green", lwd = 3)
lines(x, dlnorm(x, 5,sigma), col = "blue", lwd = 3)
lines(x, dlnorm(x, 10,sigma), col = "yellow", lwd = 3)
lines(x, dlnorm(x, 20,sigma), col = "purple", lwd = 3)
lines(x, dlnorm(x, 30,sigma), col = "grey", lwd = 3)
abline(h=0, lwd=3)
hist(background, pch=20, breaks=25, prob=TRUE, add = TRUE, col = alpha("red", 0.2))
hist(target[,1], pch=20, breaks=25, prob=TRUE, add = TRUE, col = alpha("orange", 0.2))
hist(target[,2], pch=20, breaks=25, prob=TRUE, add = TRUE, col = alpha("cyan", 0.2))
hist(target[,3], pch=20, breaks=25, prob=TRUE, add = TRUE, col = alpha("green", 0.2))
hist(target[,4], pch=20, breaks=25, prob=TRUE, add = TRUE, col = alpha("blue", 0.2))
hist(target[,5], pch=20, breaks=25, prob=TRUE, add = TRUE, col = alpha("yellow", 0.2))
hist(target[,6], pch=20, breaks=25, prob=TRUE, add = TRUE, col = alpha("purple", 0.2))
hist(target[,7], pch=20, breaks=25, prob=TRUE, add = TRUE, col = alpha("grey", 0.2))
box()
# legend("topright", c("Diff_mu=0", "Diff_mu=1", "Diff_mu=3", "Diff_mu=5", "Diff_mu=10","Diff_mu=20","Diff_mu=30"),
#        col = c("red", "orange", "cyan","green","blue","yellow","purple","grey"), lty = c(rep(1,8)), cex=0.8)

