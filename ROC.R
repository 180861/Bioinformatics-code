install.packages("pROC")
library(pROC)
roc1 <- roc(Expression ~ Group, smooth = TRUE, percent = TRUE, mydata)
plot(roc1)