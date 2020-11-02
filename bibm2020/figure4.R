source("C:/Users/sypark/Desktop/Projects/w_Shujaat/0ovltools-GSSMD/2code/src/gssmd_source.R")

# Case study #1
# Application of GSSMD in assay quality control
library(fitdistrplus)
library(cellHTS2)
experimentName <- "KcViab"
dataPath <- system.file(experimentName, package="cellHTS2")
dataPath
rev(dir(dataPath))[1:12]

# platelist.txt: The set of available result fils - filenames, Plate, and Replicate
x <- readPlateList("Platelist.txt", name = experimentName, path = dataPath)
x <- configure(x, 
               descripFile = "Description.txt",
               confFile = "Plateconf.txt",
               logFile = "Screenlog.txt",
               path = dataPath)


xn <- normalizePlates(x,
                      scale = "multiplicative",
                      log = TRUE,
                      method = "median",
                      varianceAdjust = "none")


background1 <- xn@assayData$`Channel 1`[,1][xn@featureData@data$controlStatus == "neg"]
target1 <- xn@assayData$`Channel 1`[,1][xn@featureData@data$controlStatus == "pos"]
sample1 <- xn@assayData$`Channel 1`[,1][xn@featureData@data$controlStatus == "sample"]
background2 <- xn@assayData$`Channel 1`[,2][xn@featureData@data$controlStatus == "neg"]
target2 <- xn@assayData$`Channel 1`[,2][xn@featureData@data$controlStatus == "pos"]
sample2 <- xn@assayData$`Channel 1`[,2][xn@featureData@data$controlStatus == "sample"]

# histogram
fit_b <- fitdist(background1, "norm")
fit_t <- fitdist(target1, "norm")
hist(background1, prob=TRUE, breaks = 8, col = rgb(1,0,0,0.3), ylab = "Probability density", 
     xlab = "Luciferase activity", xlim = c(-3.5,1.5), ylim = c(0,6), main = "")
hist(target1, prob=TRUE, breaks = 8, col = rgb(0,0,1,0.3), add=TRUE)
abline(v=qnorm(0.005,fit_b$estimate[1], fit_b$estimate[2]), lwd=3)
legend("topleft", c(paste("SSMD =", round(qc.measure(target1, background1, "ssmd"),3)), 
                    paste("GSSMD =", round(qc.measure(target1, background1, "gcnr"),3)),
                    paste("simple threshold =", round(qnorm(0.005,fit_b$estimate[1], fit_b$estimate[2]),3)),
                    paste("LR threshold =", round(thres_logit(fit$coefficients[1], fit$coefficients[2], p = 0.5),3))))
box()


# train classifier
# training data
x <- c(background1, target1)
y <- c(rep(0,length(background1)), rep(1,length(target1)))
train_data <- data.frame(x, y)

# test data
x <- c(background2, target2)
y <- c(rep(0,length(background2)), rep(1,length(target2)))
test_data <- data.frame(x, y)

## training and test classical LR classifier
fit <- glm(y ~ x, train_data, family = binomial)

## prediction with classical classifier
pred <- predict(fit, newdata = test_data, type = "response")
label <- pred>0.5
accuracy <- mean(label == test_data$y)
cat("Overall Accuracy: ", accuracy, '\n')

ind0 <- which(test_data$y == 0)
typeI <- mean(label[ind0] != test_data$y[ind0])
cat("Type I error: ", typeI, '\n')

plot(test_data[test_data$y==1,]$x, col = "blue", pch=16, ylim = c(-4, 4), xlab = "sample", ylab = "x1")
points(test_data[test_data$y==0,]$x, col = "red", pch = 16)
legend("topright", legend = c("Class 1", "Class 0"), col = c("blue", "red"), pch = rep(16,2))
legend("topleft", legend = c(paste("Accuracy =", round(accuracy,3)), paste("Type I error =", round(typeI,3))))

thres_logit <- function(a, b, p) {-1/b*(log((1/p)-1)+a)}



## Hit selection based on GSSMD
## Upperbound
library(edfun)
fit_b <- fitdist(background1, "norm")
fit_t <- fitdist(target1, "norm")
plot(ecdf(target1), main = "Target ECDF")
abline(h=qc.measure(background1, target1,"gssmd"), lwd=3, lty = 2, col = rgb(0,1,0,0.5))
#abline(h=0.95, lwd=3, lty = 2, col = rgb(1,0,0,0.5))
#abline(v=edfun::edfun(target1)$qfun(0.95), lwd=3, lty = 2, col = rgb(1,0,0,0.5))
abline(v=edfun::edfun(target1)$qfun(abs(qc.measure(target1, background1,"gssmd"))), lwd=3, lty = 2, col = rgb(0,1,0,0.5))
gssmd_upper = edfun::edfun(target1)$qfun(abs(qc.measure(target1, background1,"gssmd")))
pgssmd_upper = qnorm(qc.measure(background1, target1,"pgssmd"),fit_t$estimate[1], fit_t$estimate[2])

## Lowerbound
plot(ecdf(background1), main = "Background ECDF")
abline(h=1-qc.measure(background1, target1,"gssmd"), lwd=3, lty = 2, col = rgb(0,1,0,0.5))
#abline(h=0.05, lwd=3, lty = 2, col = rgb(1,0,0,0.5))
#abline(v=edfun::edfun(background1)$qfun(0.05), lwd=3, lty = 2, col = rgb(1,0,0,0.5))
abline(v=edfun::edfun(background1)$qfun(1-abs(qc.measure(target1,background1,"gssmd"))), lwd=3, lty = 2, col = rgb(0,1,0,0.5))
gssmd_lower = edfun::edfun(background1)$qfun(1-abs(qc.measure(target1,background1,"gssmd")))
pgssmd_lower = qnorm(1-qc.measure(background1, target1,"pgssmd"),fit_b$estimate[1], fit_b$estimate[2])

## Threshold
thres_gssmd = gssmd_threshold(target1, background1)
thres_pgssmd = (pgssmd_lower+pgssmd_upper)/2
thres_classifier = thres_logit(fit$coefficients[1], fit$coefficients[2], p = 0.5)


## Hit selection based on SSMD
## Upperbound
ssmd_upper = samp_ssmd(sample1, background1,-5, type="threshold")
ssmd_robust_upper = samp_ssmd(sample1, background1,-5, type="threshold_robust")

## lowerbound
ssmd_lower = samp_ssmd(sample1, background1, 0, type="threshold")
ssmd_robust_lower = samp_ssmd(sample1, background1, 0, type="threshold_robust")

## SSMD criteria
ssmd_thres = samp_ssmd(sample1, background1, -3, type="threshold")
ssmd_robust_thres = samp_ssmd(sample1, background1, -3, type="threshold_robust")


## Visualize GSSMD threshold
range <- seq(0,21660, length.out = length(target2))
plot(sample2, col = "grey", ylab = "Luciferase activity", xlab = "Samples", ylim=c(-3.5, 1))
points(range, background2, pch=16, col = "red")
points(range, target2, pch=16, col = "blue")
abline(h=gssmd_upper, lwd=3, lty = 1, col = rgb(0,1,0,0.5))
abline(h=gssmd_lower, lwd=3, lty = 1, col = rgb(0,1,0,0.5))
abline(h=thres_gssmd, lwd=3, lty = 2, col = rgb(0,1,0,0.5))
abline(h=pgssmd_upper, lwd=3, lty = 1, col = rgb(0,0,1,0.5))
abline(h=pgssmd_lower, lwd=3, lty = 1, col = rgb(0,0,1,0.5))
abline(h=thres_pgssmd, lwd=3, lty = 2, col = rgb(0,0,1,0.5))
abline(h=thres_classifier, lwd = 2, lty = 2, col = "black")
legend("bottomright", 
       c(paste("GSSMD_Upper_bound =", round(gssmd_upper,3)),
         paste("GSSMD_Lower_bound =", round(gssmd_lower,3)),
         paste("PGSSMD_Upper_bound =", round(pgssmd_upper,3)),
         paste("PGSSMD_Lower_bound =", round(pgssmd_lower,3)),
         paste("GSSMD_threshold =", round(thres_gssmd,3)),
         paste("PGSSMD_threshold =", round(thres_pgssmd,3)),
         paste("LR threshold =", round(thres_classifier,3))),
       lty = c(rep(c(1,1,2),2), 2), 
       col = c(rep(rgb(0,1,0,0.5),3),rep(rgb(0,0,1,0.5),3), "black"), 
       lwd = rep(3,7))



## Visualize SSMD threshold
fit_b <- fitdist(background1, "norm")
fit_t <- fitdist(target1, "norm")
range <- seq(0,21660, length.out = length(target2))
plot(sample2, col = "grey", ylab = "Luciferase activity", xlab = "Samples", ylim=c(-3.5, 1))
points(range, background2, pch=16, col = "red")
points(range, target2, pch=16, col = "blue")
abline(h=ssmd_upper, lwd=3, lty = 1, col = rgb(0,1,0,0.5))
abline(h=ssmd_lower, lwd=3, lty = 1, col = rgb(0,1,0,0.5))
abline(h=ssmd_thres, lwd=3, lty = 2, col = rgb(0,1,0,0.5))
abline(h=ssmd_robust_upper, lwd=3, lty = 1, col = rgb(0,0,1,0.5))
abline(h=ssmd_robust_lower, lwd=3, lty = 1, col = rgb(0,0,1,0.5))
abline(h=ssmd_robust_thres, lwd=3, lty = 2, col = rgb(0,0,1,0.5))
abline(h=thres_classifier, lwd = 2, lty = 2, col = "black")
legend("bottomright",
       c(paste("ssmd_upper =", round(ssmd_upper,3)),
         paste("ssmd_lower =", round(ssmd_lower,3)),
         paste("ssmd_robust_upper =", round(ssmd_robust_upper,3)),
         paste("ssmd_robust_lower =", round(ssmd_robust_lower,3)),
         paste("ssmd_thres =", round(ssmd_thres,3)),
         paste("ssmd_robust_thres =", round(ssmd_robust_thres,3)),
         paste("LR threshold =", round(thres_classifier, 3))),
       lty = c(rep(c(1,1,2), 2), 2), 
       col = c(rep(rgb(0,1,0,0.5),3),rep(rgb(0,0,1,0.5),3), "black"), 
       lwd = rep(3,7))


