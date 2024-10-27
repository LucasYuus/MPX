rm(list = ls())
setwd("")  # Set working directory
rt = dat2
rt <- read.table("", sep = "\t", check.names = F, stringsAsFactors = F, header = T, row.names = 1)
library(glmnet)
library(survival)

str(rt)
### 3. Model Construction (Prognostic Model)
set.seed(2)  # Set random seed
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$OS.time,rt$OS))
fit=glmnet(x, y, family = "cox",maxit = 10000)

plot(fit, xvar = "lambda", label = TRUE)
cvfit = cv.glmnet(x,y,family="cox",maxit = 10000)

plot(cvfit)
abline(v = log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
coef=coef(fit, s = cvfit$lambda.min)

index = which(coef != 0)
actCoef = coef[index]
lassoGene = row.names(coef)[index]
geneCoef = cbind(Gene = lassoGene, Coef = actCoef)
geneCoef  # View model-related coefficients
write.table(geneCoef, "", row.names = TRUE, col.names = TRUE, sep = ",")

# Modeling (Training Set)
FinalGeneExp = rt[lassoGene]
myFun = function(x) { crossprod(as.numeric(x), actCoef) }
riskScore = apply(FinalGeneExp, 1, myFun)
outCol = c(lassoGene)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
trainRiskOut = cbind(rt[, c("OS.time", "OS", outCol)], riskScore = as.vector(riskScore), risk)
write.table(trainRiskOut, "", row.names = TRUE, col.names = TRUE, sep = ",")

# Validation (Validation Set)
rt1 <- read.table("", sep = "\t", check.names = F, stringsAsFactors = F, header = T, row.names = 1)
FinalGeneExp = rt1[, lassoGene]
myFun = function(x) { crossprod(as.numeric(x), actCoef) }
riskScore = apply(FinalGeneExp, 1, myFun)
outCol = c(lassoGene)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
trainRiskOut = cbind(rt1[, c("", "", outCol)], riskScore = as.vector(riskScore), risk)
write.table(trainRiskOut, "", row.names = TRUE, col.names = TRUE, sep = ",")

# ROC Analysis
# install.packages("pROC")
library(pROC)      
library(ggplot2)  

gfit <- roc(trainRiskOut$OS, trainRiskOut$riskScore)

plot(gfit,
     print.auc = TRUE,
     auc.polygon = TRUE, 
     auc.polygon.col = "white",  # Set ROC curve fill color
     smooth = F,  # Draw unsmoothed curve
     col = "#ff1107",  # Curve color
     legacy.axes = TRUE)  # Set x-axis to range from 0 to 1, representing 1 - specificity

# Time-Dependent ROC Analysis
library(timeROC)
library(survival)
ROC3 <- timeROC(T = trainRiskOut$OS.time,  # Outcome time
                delta = trainRiskOut$OS,   # Outcome indicator
                marker = trainRiskOut$riskScore,  # Predictor variable
                cause = 1,   # Positive outcome indicator value
                weighting = "marginal",   # Calculation method, default is marginal
                times = c(1, 5, 7),  # Time points, 3-year, 5-year, and 7-year survival rates
                iid = TRUE)
plot(ROC3,
     time = 1, col = "blue")   # time specifies the time point, col specifies line color
plot(ROC3,
     time = 5, col = "red", add = TRUE)  # add specifies whether to add to the previous plot
plot(ROC3,
     time = 7, col = "green", add = TRUE)
legend("bottomright",
       c("Year 3", "Year 5", "Year 7"),
       col = c("blue", "red", "green"),
       lwd = 1)