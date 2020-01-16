# Figure 5
source("C:/Users/sypark/Desktop/0PJ2-GSSMD/gssmd_source.R")

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

## sample SSMD
plot(samp_ssmd(sample1, background1, "mean"), col = "grey", ylab = "sample SSMD", xlab = "sample", ylim = c(-150, 50))
points(range, samp_ssmd(background1, background1, "mean"), pch=16, col = "red")
points(range, samp_ssmd(target1, background1, "mean"), pch=16, col = "blue")
abline(h=-5, lwd=3, lty = 2)
legend("topleft", legend = "sample SSMD threshold = 5")

## sample Z-factor
plot(sample1, col = "grey", ylab = "Luciferase activity", xlab = "sample", ylim = c(-5.5, 1.8))
points(range, background1, pch=16, col = "red")
points(range, target1, pch=16, col = "blue")
abline(h=zfactor(sample1, background1), lwd=3, lty = 2)
legend("topleft", legend = paste("z-factor threshold =", round(zfactor(sample1, background1),3)))

# QC measure for replicate 1
fit_b <- fitdist(background1, "norm")
range <- seq(0,21660, length.out = length(target1))
plot(sample1, col = "grey", ylab = "Luciferase activity", xlab = "Samples", ylim=c(-5.5, 1.8))
points(range, background1, pch=16, col = "red")
points(range, target1, pch=16, col = "blue")
abline(h=qnorm(0.005,fit_b$estimate[1], fit_b$estimate[2]), lwd=3, lty = 2, col = rgb(0,0,1,0.5))
abline(h=samp_ssmd(sample1, background1,"threshold"), lwd=3, lty = 2, col = rgb(0,1,0,0.5))
abline(h=zfactor(sample1, background1), lwd=3, lty = 2, col = rgb(1,0,0,0.5))
abline(h=thres_logit(fit$coefficients[1], fit$coefficients[2], p = 0.5), lwd = 3, lty = 2)
legend("topleft", c(paste("SSMD =", round(qc.measure(target1, background1, "ssmd"),3)), 
                    paste("GSSMD =", round(qc.measure(target1, background1, "gcnr"),3))))
legend("topright", c("LR Classifier accuracy = 0.956", "Type I error = 0.018"))
legend("bottomright", 
       c(paste("Z-factor threshold =", round(zfactor(sample1, background1),3)),
         paste("SSMD threshold =", round(samp_ssmd(sample1, background1,"threshold"),3)),
         paste("GSSMD threshold =", round(qnorm(0.005,fit_b$estimate[1], fit_b$estimate[2]),3)),
         paste("LR threshold =", round(thres_logit(fit$coefficients[1], fit$coefficients[2], p = 0.5),3))),
       lty = rep(2,4), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5), rgb(0,0,1,0.5), "black"), lwd = rep(3,4))

# QC measure for replicate 2
fit_b <- fitdist(background2, "norm")
range <- seq(0,21660, length.out = length(target2))
plot(sample2, col = "grey", ylab = "Luciferase activity", xlab = "Samples", ylim=c(-4, 1.8))
points(range, background2, pch=16, col = "red")
points(range, target2, pch=16, col = "blue")
abline(h=qnorm(0.005,fit_b$estimate[1], fit_b$estimate[2]), lwd=3, lty = 2, col = "darkgrey")
abline(h=thres_logit(fit$coefficients[1], fit$coefficients[2], p = 0.5), lwd = 3, lty = 2)
legend("topleft", c(paste("SSMD =", round(qc.measure(target1, background1, "ssmd"),3)), 
                    paste("GSSMD =", round(qc.measure(target1, background1, "gcnr"),3)),
                    paste("simple threshold =", round(qnorm(0.005,fit_b$estimate[1], fit_b$estimate[2]),3)),
                    paste("LR threshold =", round(thres_logit(fit$coefficients[1], fit$coefficients[2], p = 0.5),3))))

