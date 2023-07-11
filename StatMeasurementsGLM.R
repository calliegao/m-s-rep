library(ggplot2)

Data_CNvspMCI_t0 <- read.csv("E:/new1104/revision/round2/form/Data/Data_CNvspMCI_t0.csv",header=FALSE, sep=",")
Data_CNvspMCI_t1 <- read.csv("E:/new1104/revision/round2/form/Data/Data_CNvspMCI_t1.csv",header=FALSE, sep=",")
Data_CNvspMCI_t2 <- read.csv("E:/new1104/revision/round2/form/Data/Data_CNvspMCI_t2.csv",header=FALSE, sep=",")
Data_CNvssMCI_t0 <- read.csv("E:/new1104/revision/round2/form/Data/Data_CNvssMCI_t0.csv",header=FALSE, sep=",")
Data_CNvssMCI_t1 <- read.csv("E:/new1104/revision/round2/form/Data/Data_CNvssMCI_t1.csv",header=FALSE, sep=",")
Data_CNvssMCI_t2 <- read.csv("E:/new1104/revision/round2/form/Data/Data_CNvssMCI_t2.csv",header=FALSE, sep=",")
Data_sMCIvspMCI_t0 <- read.csv("E:/new1104/revision/round2/form/Data/Data_sMCIvspMCI_t0.csv",header=FALSE, sep=",")
Data_sMCIvspMCI_t1 <- read.csv("E:/new1104/revision/round2/form/Data/Data_sMCIvspMCI_t1.csv",header=FALSE, sep=",")
Data_sMCIvspMCI_t2 <- read.csv("E:/new1104/revision/round2/form/Data/Data_sMCIvspMCI_t2.csv",header=FALSE, sep=",")
Ran_CNvspMCI <- read.csv("E:/new1104/revision/round2/form/Data/LongMeatures_CNvspMCI.csv",header=FALSE, sep=",")
Ran_CNvssMCI <- read.csv("E:/new1104/revision/round2/form/Data/LongMeatures_CNvssMCI.csv",header=FALSE, sep=",")
Ran_sMCIvspMCI <- read.csv("E:/new1104/revision/round2/form/Data/LongMeatures_sMCIvspMCI.csv",header=FALSE, sep=",")

pvalue <- matrix(0, nrow = 29, ncol = 1)
pvalue <- as.array(pvalue)
#fvalue <- matrix(0, nrow = 29, ncol = 1)
#fvalue <- as.array(fvalue)
d <- matrix(0, nrow = 29, ncol = 1)
d <- as.array(d)
means_Group2 <- matrix(0, nrow = 29, ncol = 1)
means_Group2 <- as.array(means_Group2)
means_Group1 <- matrix(0, nrow = 29, ncol = 1)
means_Group1 <- as.array(means_Group1)
sd_Group2 <- matrix(0, nrow = 29, ncol = 1)
sd_Group2 <- as.array(sd_Group2)
sd_Group1 <- matrix(0, nrow = 29, ncol = 1)
sd_Group1 <- as.array(sd_Group1)

for (i in 2:30) {
   Features <- Ran_sMCIvspMCI[,i]
   Group <- Ran_sMCIvspMCI[,1]
   Age <- Ran_sMCIvspMCI[,31] 
   Sex <- Ran_sMCIvspMCI[,32] 
   fit <- glm(Features ~ Group + Age + Sex, data = Ran_sMCIvspMCI)
   #summary(fit)
   mean_2 <- mean(Features[Group == 2])
   mean_1 <- mean(Features[Group == 1])
   sd_all <- sd(Features)
   d[i-1,] <- (mean_1 - mean_2) / sd_all
   pvalue[i-1,] <- summary(fit)$coefficients[2, 4]
   #fvalue[i-1,] <- summary(fit)$fstatistic[1]
   means_Group2[i-1,] <- mean_2
   means_Group1[i-1,] <- mean_1
   sd_Group2[i-1,] <- sd(Features[Group == 2])
   sd_Group1[i-1,] <- sd(Features[Group == 1])
}
adjustedpvalue <- p.adjust(pvalue, method = "fdr")

stat <- data.frame(means_Group1, sd_Group1, means_Group2, sd_Group2, adjustedpvalue, d)
write.csv(stat, file = "E:/new1104/revision/round2/form/Data/stat_Long_sMCIvspMCI.csv", row.names = FALSE)

