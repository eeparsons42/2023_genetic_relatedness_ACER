
library(ggplot2)

#Reading in dataset of bitwise distance with fertilization data
genetics <- both.nights.fert.bit <- read.csv("C:/Users/parsonsee/Desktop/genetic-analysis/data/both.nights.fert.bit.csv")

genetics


#fitting all data to a linear model 
fit <- lm(genetics$Fert~genetics$Pairwise.Distance)
fit
summary(fit)



#plot with linear regression line and 95% confidence interval 
bothnightsplot <- ggplot(genetics, aes(x=Pairwise.Distance, y=Fert)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line 
#  (by default includes 95% confidence region)
bothnightsplot



#subsetting data to include only night 1
night.1 <- subset(genetics, night=="Night 1")
fit1 <- lm(night.1$Fert~night.1$Pairwise.Distance)
summary(fit1)

#plot with linear regression line and 95% confidence interval 
night1plot <- ggplot(night.1, aes(x=Pairwise.Distance, y=Fert)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line 
#  (by default includes 95% confidence region)
night1plot


#subsetting data to include only night 2 
night.2 <- subset(genetics, night=="Night 2")
fit2 <- lm(night.2$Fert~night.2$Pairwise.Distance)
summary(fit2)

#plot with linear regression line and 95% confidence interval 
night2plot <- ggplot(night.2, aes(x=Pairwise.Distance, y=Fert)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line 
#  (by default includes 95% confidence region)
night2plot

# Set color by cond
ggplot(genetics, aes(x=Pairwise.Distance, y=Fert, color=dir)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line 
#  (by default includes 95% confidence region)




dir1 <- ggplot(night.1, aes(x=Pairwise.Distance, y=Fert, color=dir)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line 
#  (by default includes 95% confidence region)
dir1

dir2 <- ggplot(night.2, aes(x=Pairwise.Distance, y=Fert, color=dir)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line 
#  (by default includes 95% confidence region)
dir2


