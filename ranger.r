## R CMD BATCH --vanilla 'ranger.r' &

##  also attempted:
# method = "gbm_h2o"
# method = "avNNet"
# method = "gbm"

## ----load-------------------------------------------------------------
load("dataset.rda")

## ----load.joehanes -------------------------------------------------------------
load("joehanes2016_st2_bonf.rda")

## ----subset.meth -------------------------------------------------------------
X <- meth[joehanes$probe.id, ]
## transpose to make columns = methylation site variables,
##					rows = subjects/observations
X <- t(X) 

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

system.time(rf.fit <- train(y = as.factor(training$ever.smoke), 
				x = X[in.train,], 
                method = "glmnet"))

rf.fit

a<-proc.time(); a;
set.seed(825)
fit.rf <- caret::train(y = as.factor(training$ever.smoke), 
				x = X[in.train,], 
                method = "ranger")
b<-proc.time()-a; b;
(time <- paste0("Computation time = ", round(b[3]/60, 2)," minutes"))

save(list = c("fit.rf"), file = "fit.rf.rda")
