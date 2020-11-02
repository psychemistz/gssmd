## Figure1 (N(mu, sigma))

source("C:/Users/sypark/Desktop/Projects/w_Shujaat/0ovltools-GSSMD/2code/src/gssmd_source.R")

# parameter
sigma = 5
xlim = c(-15, 45)
# xlim = c(-2, 33) for sigma=1
ylim = c(0,0.13) # for sigma=1
savename = "table3.csv"

# Sample
set.seed(123)
m_list <- c(0,1,3,5,10,20,30)
target <- lapply(m_list, FUN = rnorm, n = 1000, sd = sigma)
target <- do.call(cbind, target)
colnames(target) <- c("m0", "m1", "m3", "m5", "m10", "m20", "m30")
background <- rnorm(1000, 0, sigma)

# calculate measures
zpfactor <- ssmd <- gssmd <- c(); 

for (i in 1:length(m_list)){
  zpfactor[i] <- qc.measure(target[,i], background, "zpfactor")
  ssmd[i] <- qc.measure(target[,i], background, "ssmd")
  #robust_ssmd[i] <- qc.measure(target[,i], background, "robust_ssmd")
  #gcnr[i] <- qc.measure(target[,i], background, "gcnr")
  gssmd[i] <- qc.measure(target[,i], background, "gssmd")
}

# print table.1
library(knitr)
library(magrittr)
library(kableExtra)
res <- round(data.frame(zpfactor, ssmd, gssmd),3)
rownames(res) <- c("Diff_mu=0", "Diff_mu=1", "Diff_mu=3", "Diff_mu=5", "Diff_mu=10","Diff_mu=20","Diff_mu=30")
kable(res, format = 'latex', booktabs = TRUE) %>% kable_styling(position = "center")
write.csv(res, paste0("C:/Users/sypark/Desktop/Projects/w_Shujaat/0ovltools-GSSMD/3results/", savename))

# print figure.1
library(scales)
x <- seq(-40, 40, 0.05)

plot(x, dnorm(x, 0,sigma), col = "red", lwd = 3, type = 'l', xlim = xlim, ylim = ylim, ylab = "", xlab = "")
lines(x, dnorm(x, 0,sigma), col = "orange", lwd = 3)
lines(x, dnorm(x, 1,sigma), col = "cyan", lwd = 3)
lines(x, dnorm(x, 3,sigma), col = "green", lwd = 3)
lines(x, dnorm(x, 5,sigma), col = "blue", lwd = 3)
lines(x, dnorm(x, 10,sigma), col = "yellow", lwd = 3)
lines(x, dnorm(x, 20,sigma), col = "purple", lwd = 3)
lines(x, dnorm(x, 30,sigma), col = "grey", lwd = 3)
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

