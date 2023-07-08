genetics <- both.nights.fert.bit <- read.csv("C:/Users/parsonsee/Desktop/genetic-analysis/data/both.nights.fert.bit.csv")

genetics

fit <- lm(genetics$Fert~genetics$Pairwise.Distance)
fit
summary(fit)
plot(genetics$Fert~genetics$Pairwise.Distance) +

#add linear trend
lines(predict(lm(genetics$Fert~genetics$Pairwise.Distance)),col='green')

lines(genetics$Fert, fitted(fit), col="blue")

night.1 <- subset(genetics, night=="Night 1")
fit1 <- lm(night.1$Fert~night.1$Pairwise.Distance)
summary(fit1)
plot1 <- plot(night.1$Fert~night.1$Pairwise.Distance, main = "Night 1 Crosses: Fertilization results and Pairwise distance", xlab = "Pairwise distance", ylab = "Fertilization"
)


night.2 <- subset(genetics, night=="Night 2")
fit2 <- lm(night.2$Fert~night.2$Pairwise.Distance, )
summary(fit2)
plot(night.2$Fert~night.2$Pairwise.Distance, main = "Night 2 Crosses: Fertilization results and Pairwise distance", xlab = "Pairwise distance", ylab = "Fertilization"
)

