#Plot validation loss

hist1 <- read.csv(file ="history1.csv", sep=',', header=TRUE)
hist1 <- hist1[,-1]
hist1$Epoch <- 1:nrow(hist1)
hist1$Model <- "Adam_0.001_256"

hist2 <- read.csv(file ="history2.csv", sep=',', header=TRUE)
hist2 <- hist2[,-1]
hist2$Epoch <- 1:nrow(hist2)
hist2$Model <- "Adam_0.01_256"

hist3 <- read.csv(file ="history3.csv", sep=',', header=TRUE)
hist3 <- hist3[,-1]
hist3$Epoch <- 1:nrow(hist3)
hist3$Model <- "Adam_0.0001_256"

hist4 <- read.csv(file ="history4.csv", sep=',', header=TRUE)
hist4 <- hist4[,-1]
hist4$Epoch <- 1:nrow(hist4)
hist4$Model <- "RMSProp_0.001_256"

hist5 <- read.csv(file ="history5.csv", sep=',', header=TRUE)
hist5 <- hist5[,-1]
hist5$Epoch <- 1:nrow(hist5)
hist5$Model <- "RMSProp_0.0001_256"

hist6 <- read.csv(file ="history6.csv", sep=',', header=TRUE)
hist6 <- hist6[,-1]
hist6$Epoch <- 1:nrow(hist6)
hist6$Model <- "SGD_0.001_256"

hist7 <- read.csv(file ="history7.csv", sep=',', header=TRUE)
hist7 <- hist7[,-1]
hist7$Epoch <- 1:nrow(hist7)
hist7$Model <- "SGD_0.0001_256"

hist8 <- read.csv(file ="history8.csv", sep=',', header=TRUE)
hist8 <- hist8[,-1]
hist8$Epoch <- 1:nrow(hist8)
hist8$Model <- "Adam_0.001_512"

hist9 <- read.csv(file ="history9.csv", sep=',', header=TRUE)
hist9 <- hist9[,-1]
hist9$Epoch <- 1:nrow(hist9)
hist9$Model <- "RMPSProp_0.001_512"

hist10 <- read.csv(file ="history10.csv", sep=',', header=TRUE)
hist10 <- hist10[,-1]
hist10$Epoch <- 1:nrow(hist10)
hist10$Model <- "SGD_0.001_512"


dfAdam <- rbind(hist1, hist2, hist3, hist8)

dfRMS <- rbind(hist4, hist5, hist9)

dfSGD <- rbind(hist6, hist7, hist10)

dfBest <- rbind(hist1, hist9, hist6)

library(ggplot2)
ggplot(dfAdam, aes(x = Epoch, y = val_loss, color = Model, linetype = Model)) +
  geom_line() +
  xlab("epoch") +
  ylab("loss") +
  ggtitle("Model Validation Loss") +
  theme(plot.title = element_text(hjust = 0.5))

library(ggplot2)
ggplot(dfRMS, aes(x = Epoch, y = val_loss, color = Model, linetype = Model)) +
  geom_line() +
  xlab("epoch") +
  ylab("loss") +
  ggtitle("Model Validation Loss") +
  theme(plot.title = element_text(hjust = 0.5))

library(ggplot2)
ggplot(dfSGD, aes(x = Epoch, y = val_loss, color = Model, linetype = Model)) +
  geom_line() +
  xlab("epoch") +
  ylab("loss") +
  ggtitle("Model Validation Loss") +
  theme(plot.title = element_text(hjust = 0.5))

library(ggplot2)
ggplot(dfBest, aes(x = Epoch, y = val_loss, color = Model, linetype = Model)) +
  geom_line() +
  xlab("epoch") +
  ylab("loss") +
  ggtitle("Best Performance by Optimizer") +
  theme(plot.title = element_text(hjust = 0.5))


