## Figure 6
library(parallel)

## Estimation of GSSMD from 10^6 number of random sampling with N(0,1). 
gssmd.trial.cl <- function(n, num.cores){
  res <- mclapply(n, function(i){
    target <- rnorm(n, 0, 1)
    background <- rnorm(n, 0, 1)
    nsample <- length(target)
    nbins <- round(1+log2(nsample))
    x <- range(min(target, background), max(target, background))
    x <- seq(x[1], x[2], length.out=nbins)
    mu_t <- mean(target)
    mu_b <- mean(background)
    pdf_t <- hist(target, x, plot=FALSE)
    pdf_b <- hist(background, x, plot=FALSE)
    OVL <- sum(apply((rbind(pdf_t$counts/sum(pdf_t$counts), pdf_b$counts/sum(pdf_b$counts))),2,min))
    GSSMD <- (mu_t - mu_b)/abs(mu_t - mu_b) * (1-OVL)
    return(GSSMD)
  }, mc.cores=num.cores)
  return(res)
}

## define cluster object
cl <- makeCluster(paste0("node", 40:49))
clusterEvalQ(cl, library(parallel))
n <- c(3, 10, 100, 1000, 10000, 100000, 1000000)

## perform experiment
res.ns <- c()
for(i in 1:10000){
  res.ns <- rbind(res.ns, unlist(parLapplyLB(cl, 1000, gssmd.trial.cl, num.cores=12)))
}

close(cl)
closeAllConnections()

## Assigned all results
# res1000000 <- res.ns
# res100000 <- res.ns
# res10000 <- res.ns
# res1000 <- res.ns
# res100 <- res.ns
# res10 <- res.ns
# res3 <- res.ns

res <- cbind(res3, res10, res100, res1000, res10000, res100000, res1000000)
colnames(res) <- n

## draw violin plot
library(ggplot2)
library(tidyr)
gssmd_sample <- tidyr::gather(as.data.frame(res), key="N", value="GSSMD")
gssmd_sample$N <- factor(gssmd_sample$N, levels=unique(gssmd_sample$N))
ggplot(gssmd_sample, aes(x=N, y=GSSMD, fill=N))+
  geom_violin(trim=TRUE) +
  facet_wrap(~N, scales="free") +
  theme_bw()


apply(res, 2, mean)
apply(res, 2, sd)
apply(res, 2, max)
apply(res, 2, min)
