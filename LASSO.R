set.seed(2024)
library(glmnet)
train <- read.csv("ASlasso.csv",row.names = 1,
                  as.is = F)
x <- as.matrix(train[,-1])
y <- ifelse(train$Group == "AS", 1,0)
fit = glmnet(x, y, family = "binomial", alpha = 1)
pdf("1.pdf",width = 7,height = 6)
plot(fit,xvar="lambda",label = F)
dev.off()


cvfit = cv.glmnet(x, y,
                  nfold=10, 
                  family = "binomial",type.measure="class",alpha=1)
pdf("2.pdf",width = 8,height = 6)
plot(cvfit)
dev.off()


cvfit$lambda.min

cvfit$lambda.1se

coef_lasso <- coef(fit,s=cvfit$lambda.min)
coef_lasso
rownames(coef_lasso)[which(coef_lasso!=0)]

