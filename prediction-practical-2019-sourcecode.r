## ----globals -------------------------------------------------------------
packages <- c("pander")
lapply(packages, require, character.only=T)

## ----dir-------------------------------------------------------------
## for example:
setwd("O:/teaching/advancedEpi_shortCourse")
#rm(list=ls())

## ----load-------------------------------------------------------------
load("dataset.rda")

## ----qc1 -------------------------------------------------------------
str(samples)
table(samples$smoking)

## ----qc2 -------------------------------------------------------------
pander(table(samples$smoking))

## ----never.smoke -------------------------------------------------------------
samples$never.smoke <- sign(samples$smoking=="never")

## be sure it matches the counts we'd expected
table(samples$never.smoke)

## ----ahrr -------------------------------------------------------------
samples$cg05575921 <- meth["cg05575921", ]

## ----roc1 -------------------------------------------------------------
## load the pROC package
require("pROC")

## use the formula-based syntax of the package
roc.out <- roc(never.smoke ~ cg05575921, data = samples)

## ----roc2 -------------------------------------------------------------
roc.out # man!! 

## ----plot.roc -------------------------------------------------------------
plot.roc(roc.out)

## ----load.joehanes -------------------------------------------------------------
load("joehanes2016_st2_bonf.rda")

## ----joehanes_str-------------------------------------------------------------
str(joehanes)

## ----subse.meth -------------------------------------------------------------
meth <- meth[joehanes$probe.id, ]

## ----make.x -------------------------------------------------------------
X <- t(meth)

## ----make_coefs -------------------------------------------------------------
coefs <- joehanes$effect 
names(coefs) <- joehanes$probe.id

## ----answers1_1 -------------------------------------------------------------
y.hat <-   X %*% coefs 

## ----answers1_2 -------------------------------------------------------------
samples$y.hat <- as.vector(y.hat)

## ----answers1_3 -------------------------------------------------------------
roc.out.again <- roc(never.smoke ~ y.hat, data = samples)
roc.out.again$auc

## ----answers1_4 -------------------------------------------------------------
plot.roc(roc.out)
lines.roc(roc.out.again, col="red")

## ----answers1_5 -------------------------------------------------------------
roc.out$auc
roc.out.again$auc


## ----folds1-------------------------------------------------------------
require(caret)

set.seed(20180314)
Y <- samples$never.smoke
flds <- createFolds(Y, k = 10, list = TRUE, 
			returnTrain = FALSE)

## ----folds2-------------------------------------------------------------
str(flds)

## ----test_train -------------------------------------------------------------
test <- flds$Fold01
train <- (-test)


## ----lasso-------------------------------------------------------------
require(glmnet)
set.seed(420)
fit.lasso <- cv.glmnet(X[train,], Y[train], 
	family='binomial', alpha=1, standardize=TRUE, 
	type.measure='auc', nfolds = 10)

## ----coef.lasso -------------------------------------------------------------
coef.lasso <- as.matrix(coef(fit.lasso, s = "lambda.min"))
coef.lasso <- data.frame(betas=coef.lasso[,1], 
				stringsAsFactors = F)
coef.lasso <- subset(coef.lasso, betas!=0)

str(coef.lasso)

## ----coef-table-students -------------------------------------------------------------
coef.lasso

## ----pred.lasso -------------------------------------------------------------
pred.lasso <- as.vector(predict(fit.lasso, 
	newx = X[test,], type = "response", s = "lambda.min"))

## ----roc.lasso-------------------------------------------------------------
roc.out.l <- roc(Y[test], pred.lasso)
roc.out.l$auc

## ---- -------------------------------------------------------------
#get the coefficients from the lasso
coefs.lasso <- coef(fit.lasso, s=fit.lasso$lambda.min)

## ----answers2_1 -------------------------------------------------------------
X <- as.data.frame(X)
# new folds
test <- flds$Fold02
train <- (-test)

## Model 1: glm using a random subset of cpgs
set.seed(86)
cpgs <- sample(colnames(X), 20) # sample 20 cpgs

fit.glm <- glm(Y[train] ~. , 
	data=X[train,cpgs], family='binomial')

pred.glm <- predict(fit.glm, 
	newx = X[test,], type = "response")

## ----answers2_2 -------------------------------------------------------------

## Model 2: stepwise glm choosing from 20 random cpgs
fit.step <- step(fit.glm, direction = "backward", 
			trace = 1, k = 2)

pred.step <- predict(fit.step, newdata = X[test,], 
				type = "response")

## ----answers2_2_roc -------------------------------------------------------------
roc.step <- pROC::roc(Y[test], pred.step)
roc.step$auc

## ----answers2_3 -------------------------------------------------------------

## Model 3: knn
fit.knn <- class::knn(train = X[train,], test = X[test,], 
	k = 10, cl = Y[train], prob = TRUE)
pred.knn <- (as.numeric(fit.knn) - 1) * attr(fit.knn, "prob") +
    (1 - (as.numeric(fit.knn) - 1)) * (1 - attr(fit.knn,
        "prob"))

roc.knn <- pROC::roc(Y[test], pred.knn)
roc.knn$auc

## ----answers2_4 -------------------------------------------------------------

## Model 3: bayes regression
fit.bayes <- arm::bayesglm(Y[train] ~., data = X[train, cpgs], family = 'binomial')
pred.bayes <- predict(fit.bayes, newdata = X[test,], 
				type = "response")

roc.bayes <- pROC::roc(Y[test], pred.bayes)
roc.bayes$auc

## ----answers2_5-------------------------------------------------------------
plot.roc(roc.step, col="red")
lines.roc(roc.knn, col="blue")
lines.roc(roc.bayes, col="purple")


