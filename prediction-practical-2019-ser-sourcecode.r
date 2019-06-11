## ----globals -------------------------------------------------------------
packages <- c("pander")
lapply(packages, require, character.only=T)

## ----load-------------------------------------------------------------
load("dataset.rda")

## ----ls -------------------------------------------------------------
ls()

## ----data.save -------------------------------------------------------------

## rm(list=c("cc", "sva.age", "sva.smoking"))
## ever.smoke <- sign(!samples$smoking=="never")
## samples <- with(samples, data.frame(gsm = as.character(gsm),
## 				gse, age, sex, smoking, ever.smoke, stringsAsFactors = F))
## save(list = c("samples", "meth"), file = "dataset.rda")

## ----qc1 -------------------------------------------------------------
str(samples)
summary(samples)

table(samples$smoking)
table(samples$ever.smoke)

## ----ahrr -------------------------------------------------------------
samples$ahrr <- meth["cg05575921", ]

## ----roc1 -------------------------------------------------------------
## load the pROC package
library("pROC")

## use the formula-based syntax of the package
roc(ever.smoke ~ ahrr, data = samples)

## ----plot.roc -------------------------------------------------------------
roc.out <- roc(ever.smoke ~ ahrr, data = samples)
plot(roc.out)
plot(roc.out, 
	print.thres="best",
	print.thres.best.method="closest.topleft")

## ----load.joehanes -------------------------------------------------------------
load("joehanes2016_st2_bonf.rda")

## ----joehanes_str-------------------------------------------------------------
str(joehanes)

## ----subset.meth -------------------------------------------------------------
X <- meth[joehanes$probe.id, ]
## transpose to make columns = methylation site variables,
##					rows = subjects/observations
X <- t(X) 

## ----make_coefs -------------------------------------------------------------
coefs <- joehanes$effect 
names(coefs) <- joehanes$probe.id

## ----apply_coefs -------------------------------------------------------------
y.hat <-   X %*% coefs 

summary(y.hat)

## ----add_yhat -------------------------------------------------------------
samples$y.hat <- as.vector(y.hat)

## ----plot.roc.again -------------------------------------------------------------
roc.out.again <- roc(ever.smoke ~ y.hat, data = samples)
roc.out.again

plot(roc.out)
lines.roc(roc.out.again, col="red")

coords(roc.out.again, "best",
					 best.method="closest.topleft", 
					 ret=c("threshold", "accuracy"))


## ----comp_roc -------------------------------------------------------------
roc.out$auc
roc.out.again$auc



## ----data.partition -------------------------------------------------------------
library(caret)

set.seed(138) # makes random processes reproducible
Y <- samples$ever.smoke
in.train <- createDataPartition(
  y = samples$ever.smoke,
  ## the outcome data are needed
  p = .75,
  ## The percentage of data in the
  ## training set
  list = FALSE
)

## ----data.subset -------------------------------------------------------------
training <- samples[ in.train,]
testing  <- samples[-in.train,]

nrow(training)
nrow(testing)

## ----lasso-------------------------------------------------------------
library(glmnet)
set.seed(20)
fit.lasso <- cv.glmnet(y = training$ever.smoke, 
						x = X[in.train,],  
						family='binomial', 
						alpha=1)

## ----pred.lasso -------------------------------------------------------------
pred.lasso <- predict(fit.lasso, newx = X[-in.train,], 
					s = "lambda.min", 
					type = "response")

## ----roc.lasso-------------------------------------------------------------
roc.lasso <- roc(testing$ever.smoke, as.vector(pred.lasso))
roc.lasso
plot(roc.lasso)

## ----confusion.lasso -------------------------------------------------------------
pred.lasso <- predict(fit.lasso, newx = X[-in.train,], 
					s = "lambda.min", 
					type = "class")

caret::confusionMatrix(as.factor(pred.lasso), as.factor(testing$ever.smoke))



## ----rf -------------------------------------------------------------
set.seed(825)
rf.fit <- caret::train(y = as.factor(training$ever.smoke), 
				x = X[in.train,], 
                method = "ranger")

## ----pred.rf -------------------------------------------------------------
pred.rf <- predict(rf.fit, newdata = X[-in.train,])
confusionMatrix(pred.rf, as.factor(testing$ever.smoke))


## ----caret.models -------------------------------------------------------------
## list names of all caret models
names(getModelInfo())


## ----scratch -------------------------------------------------------------


rf.fit <- train(y = as.factor(training$ever.smoke), 
				x = X[in.train,], 
                method = "svmRadial") 
"gbm"


set.seed(825)
svmFit <- train(Class ~ ., data = training, 
                 method = "svmRadial", 
                 trControl = fitControl, 
                 preProc = c("center", "scale"),
                 tuneLength = 8,
                 metric = "ROC")
svmFit   



rf_fit <- train(Species ~ ., 
                data = iris, 
                method = "ranger")


rf_fit <- train(Species ~ ., 
                data = iris, 
                method = "glmnet")



e_fit <- train(y = as.factor(training$ever.smoke), 
				x = X[in.train,], 
                method = "glmnet")


## list names of all caret models
names(getModelInfo())


pred.enet <- predict(e_fit, newdata = X[-in.train,])
confusionMatrix(pred.enet, as.factor(testing$ever.smoke))

roc(testing$ever.smoke, as.numeric(as.character(pred.enet)))


	newx = X[test,], type = "response", s = "lambda.min")



plot(roc.out, print.thres="best", print.thres.best.method="closest.topleft")
result.coords <- coords(roc.out, "best", best.method="closest.topleft", ret=c("threshold", "accuracy"))
print(result.coords)#to get threshold and accuracy

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


## ----lasso-old-------------------------------------------------------------
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

## ----pred.lasso-old -------------------------------------------------------------
pred.lasso <- as.vector(predict(fit.lasso, 
	newx = X[test,], type = "response", s = "lambda.min"))

## ----roc.lasso.old-------------------------------------------------------------
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


