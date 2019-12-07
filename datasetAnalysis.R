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

# chi square tests
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
fit.params
curve(dgamma(x=x, shape=8.8312777, rate=1.5799002), col="blue", add=TRUE)



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
