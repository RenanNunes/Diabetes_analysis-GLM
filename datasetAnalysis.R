# read dataset
setwd("/home/renan/Renan/Lisboa - IST/Estatistica Computacional/Projeto") # change for your working directory
data = read.table('Trabalho5_EC.txt', sep = '\t', dec = '.', header = TRUE)

# checking some statistics of the data
summary(data)

# there is a small amount of numeric data missing. its going to inputed the mean for this values
for(i in 1:ncol(data)) {
  if (is.numeric(data[, i])) {
    data[, i][is.na(data[, i])] <- mean(data[, i], na.rm = TRUE)
  }
}
summary(data)

# checking correlation between variables
# first, we want the numerical variables
nums <- unlist(lapply(data, is.numeric))
nums[1] <- FALSE # we don't need to use id for the analysis
cor_data = cor(data[, nums])
cor_matrix = as.matrix(cor_data)
# install.packages("corrplot")
library(corrplot)
corrplot(cor_matrix, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
# the most correlated variables are: SGLU and GHB; H and WHT; H and W (0.83455575 - biggest one); SBP and DSP

# boxplots
boxplot(data[nums], use.cols=TRUE)

# bmi calculation (with lbs and inches)
data$BMI = 703 * data$WHT / (data$HHT*data$HHT)
# adding BMI classification
data$BMI_CAT <- cut(data$BMI, breaks = c(0, 18.5, 25, 30, 100), labels=c("Underweight", "Normal weight", "Overweight", "Obesity"))

# glycosolated hemoglobin (outcome) > 7.0 is usually taken as a positive diagnosis of diabetes
data$DIABETES <- cut(data$GHB, breaks = c(0,7,50), labels = c("Negative", "Positive"))

# ages' bin
data$AGE_CAT <- cut(data$AGE, breaks = c(17,25,35,45,55,65,200), labels=c("18-25","26-35","36-45","46-55","56-65",">=66"))

# summary with the created columns added
summary(data)

# chi square tests for a preliminar analysis
# diabetes x bmi_cat
test_bmi_cat = table(data$DIABETES, data$BMI_CAT)
test_bmi_cat
chisq.test(test_bmi_cat) # may indicate that higher BMI is related to diabetes
# fisher.test(test_bmi_cat)

# diabetes x age
test_age_cat = table(data$DIABETES, data$AGE_CAT)
test_age_cat
chisq.test(test_age_cat)

# diabetes x location
test_location = table(data$DIABETES, data$LOCATION)
test_location
chisq.test(test_location)

# diabetes x gender
test_gender = table(data$DIABETES, data$GENDER)
test_gender
chisq.test(test_gender)

# diabetes x frame
test_frame = table(data$DIABETES, data$FRAME)
test_frame
chisq.test(test_frame)

# some scatter plots with colors to differentiate HDL > 7 from HDL < 7
nums_ext <- unlist(lapply(data, is.numeric))
nums_ext[1] <- FALSE # we don't need to use id for the analysis
diabetes_color <- c("#E7B800", "#FC4E07")  
pairs(data[, nums_ext], col = diabetes_color[data$DIABETES], lower.panel=NULL) # we can see that high values of SGLU may indicate higher HDL. It makes sense, because Glycosylated hemoglobin is the Hemoglobin to which glucose is bound
# zoomed chart of SGLU x GHB
plot(data$SGLU, data$GHB, main="Glycosolated Hemoglobin x  Stabilized Glucose", xlab = "SGLU", ylab = "GHB")
# zoomed chart of CHOL x GHB
plot(data$CHOL, data$GHB, main="Glycosolated Hemoglobin x  Total Cholesterol", xlab = "CHOL", ylab = "GHB")
# zoomed chart of HDL x GHB
plot(data$HDL, data$GHB, main="Glycosolated Hemoglobin x  High Density Lipoprotein", xlab = "HDL", ylab = "GHB")

# histogram
hist(data$GHB, xlab = "GHB", main = "Histogram of GHB", freq = FALSE, breaks = 50)
library(MASS)
fit.params <- fitdistr(data$GHB, "gamma", lower = c(0, 0))
curve(dgamma(x=x, shape=fit.params$estimate['shape'], rate=fit.params$estimate['rate']), col="blue", add=TRUE)
fit.params <- fitdistr(data$GHB, "weibull", lower = c(0, 0))
curve(dweibull(x=x, shape=fit.params$estimate['shape'], scale=fit.params$estimate['scale']), col="green", add=TRUE)

# Generalized linear model (without using prior)
glm_ghb <- glm(GHB ~ CHOL + SGLU + HDL + AGE + SBP + DSP + BMI + W + H, data=data, family=Gamma(link="inverse"))
summary(glm_ghb)
# removing variables with p-value > 0.05 
glm_ghb_2 <- glm(GHB ~ CHOL + SGLU + HDL + AGE, data=data, family=Gamma(link="inverse"))
summary(glm_ghb_2)

confint(glm_ghb_2)

# anova test removing a variable that is important to check the result
glm_ghb_3 <- glm(GHB ~ SGLU + HDL + AGE + SBP + DSP + BMI + W + H, data=data, family=Gamma(link="inverse"))
anova(glm_ghb, glm_ghb_3, test = "Chisq")

# glm with prior
library("rjags")
set.seed(5)
# note: dgamma in JAGS = dgamma(shape, rate)
# model 1
mod1_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dgamma(shape, shape * inv_mu[i])
      inv_mu[i] = (b[1] + b[2]*CHOL[i] + b[3]*SGLU[i] + b[4]*HDL[i] + b[5]*AGE[i] + b[6]*SBP[i] + b[7]*DSP[i] + b[8]*BMI[i])
    }
    b[1] ~ dnorm(0.5, 1.0/1.0e4)
    for (j in 2:8) {
      b[j] ~ dnorm(0.0, 1.0/1.0e4)
    }
    shape ~ dgamma(0.001, 0.001)
} "
params = c("b", "shape")

# normalization
library(dplyr)
columns_models <- select(data, "CHOL", "SGLU", "HDL", "AGE", "SBP", "DSP", "BMI")
scaled_data <- scale(columns_models)

data_jags = list(y=data$GHB, CHOL=scaled_data[,"CHOL"], SGLU=scaled_data[,"SGLU"], HDL=scaled_data[,"HDL"], AGE=scaled_data[,"AGE"], SBP=scaled_data[,"SBP"], DSP=scaled_data[,"DSP"], BMI=scaled_data[,"BMI"])
# data_jags = list(y=data$GHB, CHOL=data[,"CHOL"], SGLU=data[,"SGLU"], HDL=data[,"HDL"], AGE=data[,"AGE"], SBP=data[,"SBP"], DSP=data[,"DSP"], BMI=data[,"BMI"])

mod1 = jags.model(textConnection(mod1_string), data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))

# convergence
plot(mod1_sim)
gelman.diag(mod1_sim)
heidel.diag(mod1_sim)

# residuals
X = cbind(rep(1.0, nrow(data)), scaled_data[,"CHOL"], scaled_data[,"SGLU"], scaled_data[,"HDL"], scaled_data[,"AGE"], scaled_data[,"SBP"], scaled_data[,"DSP"], scaled_data[,"BMI"])
(pm_params1 = colMeans(mod1_csim)) #posterior mean

# vector of predicted values from the model
yhat1 <- (drop(X %*% pm_params1[1:8]))^(-1)
# calculate residuals
resid1 <- data_jags$y - yhat1

plot(resid1)

# summary
summary(mod1_sim)

# DIC
dic.samples(mod1, n.iter = 1e3)

# model 2
mod2_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dgamma(shape, shape * inv_mu[i])
      inv_mu[i] = (b[1] + b[2]*CHOL[i] + b[3]*SGLU[i] + b[4]*HDL[i] + b[5]*AGE[i] + b[6]*DSP[i] + b[7]*BMI[i])
    }
    b[1] ~ dnorm(0.5, 1.0/1.0e4)
    for (j in 2:7) {
      b[j] ~ dnorm(0.0, 1.0/1.0e4)
    }
    shape ~ dgamma(0.001, 0.001)
} "
data_jags_2 = list(y=data$GHB, CHOL=scaled_data[,"CHOL"], SGLU=scaled_data[,"SGLU"], HDL=scaled_data[,"HDL"], AGE=scaled_data[,"AGE"], DSP=scaled_data[,"DSP"], BMI=scaled_data[,"BMI"])

mod2 = jags.model(textConnection(mod2_string), data=data_jags_2, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

# convergence
plot(mod2_sim)
gelman.diag(mod2_sim)
heidel.diag(mod2_sim)

# residuals
X2 = cbind(rep(1.0, nrow(data)), scaled_data[,"CHOL"], scaled_data[,"SGLU"], scaled_data[,"HDL"], scaled_data[,"AGE"], scaled_data[,"DSP"], scaled_data[,"BMI"])
(pm_params2 = colMeans(mod2_csim)) #posterior mean

# vector of predicted values from the model
yhat2 <- (drop(X2 %*% pm_params2[1:7]))^(-1)
# calculate residuals
resid2 <- data_jags_2$y - yhat2

plot(resid2)

# summary
summary(mod2_sim)

# DIC
dic.samples(mod2, n.iter = 1e3)

# model 3
mod3_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dgamma(shape, shape * inv_mu[i])
      inv_mu[i] = (b[1] + b[2]*CHOL[i] + b[3]*SGLU[i] + b[4]*HDL[i] + b[5]*AGE[i] + b[6]*BMI[i])
    }
    b[1] ~ dnorm(0.5, 1.0/1.0e4)
    for (j in 2:6) {
      b[j] ~ dnorm(0.0, 1.0/1.0e4)
    }
    shape ~ dgamma(0.001, 0.001)
} "
data_jags_3 = list(y=data$GHB, CHOL=scaled_data[,"CHOL"], SGLU=scaled_data[,"SGLU"], HDL=scaled_data[,"HDL"], AGE=scaled_data[,"AGE"], BMI=scaled_data[,"BMI"])

mod3 = jags.model(textConnection(mod3_string), data=data_jags_3, n.chains=3)
update(mod3, 1e3)

mod3_sim = coda.samples(model=mod3,
                        variable.names=params,
                        n.iter=5e3)
mod3_csim = as.mcmc(do.call(rbind, mod3_sim))

# convergence
plot(mod3_sim)
gelman.diag(mod3_sim)
heidel.diag(mod3_sim)

# residuals
X3 = cbind(rep(1.0, nrow(data)), scaled_data[,"CHOL"], scaled_data[,"SGLU"], scaled_data[,"HDL"], scaled_data[,"AGE"], scaled_data[,"BMI"])
(pm_params3 = colMeans(mod3_csim)) #posterior mean

# vector of predicted values from the model
yhat3 <- (drop(X3 %*% pm_params3[1:6]))^(-1)
# calculate residuals
resid3 <- data_jags_3$y - yhat3

plot(resid3)

# summary
summary(mod3_sim)

# DIC
dic.samples(mod3, n.iter = 1e3)

# model 4
mod4_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dgamma(shape, shape * inv_mu[i])
      inv_mu[i] = (b[1] + b[2]*CHOL[i] + b[3]*SGLU[i] + b[4]*AGE[i] + b[5]*BMI[i])
    }
    b[1] ~ dnorm(0.5, 1.0/1.0e4)
    for (j in 2:5) {
      b[j] ~ dnorm(0.0, 1.0/1.0e4)
    }
    shape ~ dgamma(0.001, 0.001)
} "
data_jags_4 = list(y=data$GHB, CHOL=scaled_data[,"CHOL"], SGLU=scaled_data[,"SGLU"], AGE=scaled_data[,"AGE"], BMI=scaled_data[,"BMI"])

mod4 = jags.model(textConnection(mod4_string), data=data_jags_4, n.chains=3)
update(mod4, 1e3)

mod4_sim = coda.samples(model=mod4,
                        variable.names=params,
                        n.iter=5e3)
mod4_csim = as.mcmc(do.call(rbind, mod4_sim))

# convergence
plot(mod4_sim)
gelman.diag(mod4_sim)
heidel.diag(mod4_sim)

# residuals
X4 = cbind(rep(1.0, nrow(data)), scaled_data[,"CHOL"], scaled_data[,"SGLU"], scaled_data[,"AGE"], scaled_data[,"BMI"])
(pm_params4 = colMeans(mod4_csim)) #posterior mean

# vector of predicted values from the model
yhat4 <- (drop(X4 %*% pm_params4[1:5]))^(-1)
# calculate residuals
resid4 <- data_jags_4$y - yhat4

plot(resid4)

# summary
summary(mod4_sim)

# DIC
dic.samples(mod4, n.iter = 1e3)

# model 5
mod5_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dgamma(shape, shape * inv_mu[i])
      inv_mu[i] = (b[1] + b[2]*CHOL[i] + b[3]*SGLU[i] + b[4]*HDL[i] + b[5]*AGE[i])
    }
    b[1] ~ dnorm(0.5, 1.0/1.0e4)
    for (j in 2:5) {
      b[j] ~ dnorm(0.0, 1.0/1.0e4)
    }
    shape ~ dgamma(0.001, 0.001)
} "
data_jags_5 = list(y=data$GHB, CHOL=scaled_data[,"CHOL"], SGLU=scaled_data[,"SGLU"], HDL=scaled_data[,"HDL"], AGE=scaled_data[,"AGE"])

mod5 = jags.model(textConnection(mod5_string), data=data_jags_5, n.chains=3)
update(mod5, 1e3)

mod5_sim = coda.samples(model=mod5,
                        variable.names=params,
                        n.iter=5e3)
mod5_csim = as.mcmc(do.call(rbind, mod5_sim))

# convergence
plot(mod5_sim)
gelman.diag(mod5_sim)
heidel.diag(mod5_sim)

# residuals
X5 = cbind(rep(1.0, nrow(data)), scaled_data[,"CHOL"], scaled_data[,"SGLU"], scaled_data[,"HDL"], scaled_data[,"AGE"])
(pm_params5 = colMeans(mod5_csim)) #posterior mean

# vector of predicted values from the model
yhat5 <- (drop(X5 %*% pm_params5[1:5]))^(-1)
# calculate residuals
resid5 <- data_jags_5$y - yhat5

plot(resid5)

# summary
summary(mod5_sim)

# DIC
dic.samples(mod5, n.iter = 1e3)

# model 6
mod6_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dgamma(shape, shape * inv_mu[i])
      inv_mu[i] = (b[1] + b[2]*CHOL[i] + b[3]*SGLU[i] + b[4]*AGE[i])
    }
    b[1] ~ dnorm(0.5, 1.0/1.0e4)
    for (j in 2:4) {
      b[j] ~ dnorm(0.0, 1.0/1.0e4)
    }
    shape ~ dgamma(0.001, 0.001)
} "
data_jags_6 = list(y=data$GHB, CHOL=scaled_data[,"CHOL"], SGLU=scaled_data[,"SGLU"], AGE=scaled_data[,"AGE"])

mod6 = jags.model(textConnection(mod6_string), data=data_jags_6, n.chains=3)
update(mod6, 2e3)

mod6_sim = coda.samples(model=mod6,
                        variable.names=params,
                        n.iter=5e3)
mod6_csim = as.mcmc(do.call(rbind, mod6_sim))

# convergence
plot(mod6_sim)
gelman.diag(mod6_sim)
heidel.diag(mod6_sim)

# residuals
X6 = cbind(rep(1.0, nrow(data)), scaled_data[,"CHOL"], scaled_data[,"SGLU"], scaled_data[,"AGE"])
(pm_params6 = colMeans(mod6_csim)) #posterior mean

# vector of predicted values from the model
yhat6 <- (drop(X6 %*% pm_params6[1:4]))^(-1)
# calculate residuals
resid6 <- data_jags_6$y - yhat6

plot(resid6)

# summary
summary(mod6_sim)

# DIC
dic.samples(mod6, n.iter = 1e3)

# model 7
mod7_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dgamma(shape, shape * inv_mu[i])
      inv_mu[i] = (b[1] + b[2]*SGLU[i] + b[3]*AGE[i])
    }
    b[1] ~ dnorm(0.5, 1.0/1.0e4)
    for (j in 2:3) {
      b[j] ~ dnorm(0.0, 1.0/1.0e4)
    }
    shape ~ dgamma(0.001, 0.001)
} "
data_jags_7 = list(y=data$GHB, SGLU=scaled_data[,"SGLU"], AGE=scaled_data[,"AGE"])

mod7 = jags.model(textConnection(mod7_string), data=data_jags_7, n.chains=3)
update(mod7, 1e3)

mod7_sim = coda.samples(model=mod7,
                        variable.names=params,
                        n.iter=5e3)
mod7_csim = as.mcmc(do.call(rbind, mod7_sim))

# convergence
plot(mod7_sim)
gelman.diag(mod7_sim)
heidel.diag(mod7_sim)

# residuals
X7 = cbind(rep(1.0, nrow(data)), scaled_data[,"SGLU"], scaled_data[,"AGE"])
(pm_params7 = colMeans(mod7_csim)) #posterior mean

# vector of predicted values from the model
yhat7 <- (drop(X7 %*% pm_params7[1:3]))^(-1)
# calculate residuals
resid7 <- data_jags_7$y - yhat7

plot(resid7)

# summary
summary(mod7_sim)

# DIC
dic.samples(mod7, n.iter = 1e3)
