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

## ----kfolds -------------------------------------------------------------
flds <- createFolds(
          y = samples$ever.smoke, 
          k = 10) ## number of folds 

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

## ----confusion.pred.lasso -------------------------------------------------------------
pred.lasso <- predict(fit.lasso, newx = X[-in.train,], 
					s = "lambda.min", 
					type = "class")

## ----confusion.lasso -------------------------------------------------------------
caret::confusionMatrix(as.factor(pred.lasso), as.factor(testing$ever.smoke))


## ----load.fits -------------------------------------------------------------
load("fit.rf.rda")
load("fit.svm.rda")

## ----rf -------------------------------------------------------------
set.seed(825)
fit.rf <- caret::train(y = as.factor(training$ever.smoke), 
				x = X[in.train,], 
                method = "ranger") # 2.2 minutes

## ----pred.rf -------------------------------------------------------------
pred.rf <- predict(fit.rf, newdata = X[-in.train,])
confusionMatrix(pred.rf, as.factor(testing$ever.smoke))


## ----svm -------------------------------------------------------------
set.seed(637)
fit.svm <- caret::train(y = as.factor(training$ever.smoke), 
				x = X[in.train,], 
                method = "svmRadial") ## 1.6 minutes 


## ----pred.svm -------------------------------------------------------------
pred.svm <- predict(fit.svm, newdata = X[-in.train,])
confusionMatrix(pred.svm, as.factor(testing$ever.smoke))


## ----caret.models -------------------------------------------------------------
## list names of all caret models
names(getModelInfo())




## ----scratch -------------------------------------------------------------

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



pred.enet <- predict(e_fit, newdata = X[-in.train,])
confusionMatrix(pred.enet, as.factor(testing$ever.smoke))

roc(testing$ever.smoke, as.numeric(as.character(pred.enet)))

