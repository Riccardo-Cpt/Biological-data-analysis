setwd("/home/riccardo/Desktop/progetto lauria/")

library("GEOquery")
library(genefilter)

geo <- getGEO("GSE43458", getGPL = FALSE)

#rna expression in roots treated with 5 different signalling molecules (5! combinations)
gse <- geo[[1]]
ex = exprs(gse)
dim(ex)
colnames(ex)
row.names(ex)

ex <- na.omit(as.matrix(ex)) #clean NA's
boxplot(ex)
hist(ex) #small differences between values-->we want to normalize them
min(ex);max(ex)

#median normalization


########### hist and boxplot graphycally good
hist(ex, xlab = "norm. expr. value", main =  "")
boxplot(as.data.frame(ex),
        xlab="patient",ylab="Log Mean Signal",
        axes=F)
axis(1,labels=1:ncol(ex),at=1:ncol(ex))
axis(2)
box()

#divide columns
smoker.tumor = grep("Smoker", gse$source_name_ch1)
normal = grep("Normal", gse$characteristics_ch1)
never.smoker.tumor = setdiff(grep("Never-smoker", gse$source_name_ch1), normal)

ex2 = ex[, c(smoker.tumor,never.smoker.tumor, normal)]

pca <- prcomp(t(ex2))
summary(pca)
screeplot(pca)

grpcol <- c(rep("blue",40), rep("orange",40), rep("green",30))
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2", type="p", pch=10, col=grpcol)
    text(pca$x[,1], pca$x[,2], rownames(pca$data), cex=0.75)


### k-means
k <- 2
kmeans_result <- kmeans(t(ex2), k)
table(kmeans_result$cluster)
    
library(useful)
plot(kmeans_result, data=t(ex2))


#hierarchical
dist_matrix <- dist(t(ex2))
hc_result <- hclust(dist_matrix, method = "mcq")
k <- 2
groups <- cutree(hc_result, k=k)
table(groups)
plot(hc_result, hang <- -1, labels=groups)
rect.hclust(hc_result, k = 3, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) 


#########################################################################
#random forest application
# to obtain relevant subset of data
# filter genes on raw scale, then return to log scale;
# also, keep only 30 arrays JUST for computational convenience
library("genefilter");
  
group <- c(rep('A',35),rep('B',35), rep("C", 25)) # classification, in order
group2 = c(rep('T',70), rep("C", 25))
library(randomForest)
set.seed(1357)
print(date()) #prints the current time
rf <- randomForest(x = t(ex2[,c(1:35, 41:75, 81:105)]), y = as.factor(group), ntree = 2000)
predict(rf, t(ex2[,c(36:40, 76:80, 106:110)]))

rf2 <- randomForest(x = t(ex2[,c(1:70, 81:105)]), y = as.factor(group2), ntree = 800)
predict(rf2, t(ex2[,c(71:80,106:110)]))

imp <- abs(rf$importance[,])
imp <- abs(rf2$importance[,]) # how many type a gene appears in the tree we built 
# Aka (attributes important for the classification)
t <- order(imp, decreasing = T)
plot(c(1:nrow(ex)),imp[t], log = "x", cex.lab=1.5,
      pch=16,main='ALL subset results')

important.genes = sort(rf$importance, decreasing = T)[1:200]
important.genes = sort(rf2$importance, decreasing = T)[1:200]
plot(important.genes)
probe.names = rownames(rf2$importance)
top200 = probe.names[order(rf2$importance, decreasing = T)[1:200]]

write.table(top200, "probes-top200-1tum,1cont.txt", quote = F,
            row.names = F, col.names = F)



###### comparing all tumors with normal #######
group <- c(rep('T',70),rep('N',25)) # classification, in order

library(randomForest)
set.seed(157)
print(date()) #prints the current time
rf <- randomForest(x = t(ex2[,c(1:70, 81:105)]), y = as.factor(group), ntree = 2000)
predict(rf, t(ex2[,c(71:80, 106:110)]))

imp <- abs(rf$importance[,]) # how many type a gene appears in the tree we built 
# Aka (attributes important for the classification)
t <- order(imp, decreasing = T)
plot(c(1:nrow(ex)),imp[t], log = "x", cex.lab=1.5,
     pch=16,main='ALL subset results')
#############################################################
important.genes = sort(rf$importance, decreasing = T)[1:200]
plot(important.genes)
probe.names = rownames(rf$importance)
top200 = probe.names[order(rf$importance, decreasing = T)[1:200]]

write.table(top200, "probes-top200.txt", quote = F,
            row.names = F, col.names = F)

######################################################################

# LDA

### Differences between tumor in smoker and non smoker

f <- factor(c(rep("Smoker tumor",40), rep("Non-smoker tumor",40)))
small.df = ex2[,1:80]
library("genefilter")
ttest <- rowttests(small.df,f); ttest; summary(ttest) # feature selection performing T-test filtering 
                                                      # on each row

keepers <- which(p.adjust(ttest$p.value)<0.001)
ex3 <-small.df[keepers,]
dim(ex3)
tex3 <- t(ex3)
dat <- cbind(as.data.frame(tex3),f)
colnames(dat)[ncol(dat)] <- "Individual"
smoker <- 40
non.smoker <- 40

dat$Individual #Creates a new column that labels the type of tissue

train <- sample(1:(smoker), (smoker-5)) #takes 35 numbers without replacement from the pool of 40 
test <- setdiff(1:(smoker),train) # takes the left out 5 numbers
test<- c(test, test+40)  
train <- c(train, train+40)
length(train); length(test)

library("MASS")
mod <- lda(Individual ~ ., data=dat, prior = c(0.5,0.5), 
           subset = train)
plot(mod)
View(mod)
# To test how much the model is able to separate variables

mod.values <- predict(mod, dat[test,]) #as row because previously we have tranposed
mod.values$class
plot(mod.values$x[,1], ylab=c("LDA Axis"))
text(mod.values$x[,1],
     col=c(as.numeric(dat[test,"Individual"])+10))


preds<-predict(mod, dat[test,])
preds$class
table(as.numeric(preds$class),
      as.numeric(dat[test, "Individual"]) ) #confusion matrix

library("pROC")
roc_lda <- plot.roc(as.numeric(preds$class),
                    as.numeric(dat[test, "Individual"]) )


#
# Run algorithms using 10-fold cross validation
#

library("caret")
library("e1071")

### Evaluate perfomrances of LDA and RF with caret cross validation
control <- trainControl(method="cv", number=10) # number = number of resampling iterations
metric <- "ROC"
fit.lda <- train(Individual~., data=dat, method="lda", # y = t(w) * x
                 metric=metric, trControl=control)
fit.rf <- train(Individual~., data=dat, method="rf",
                metric=metric, trControl=control)

results <- resamples(list(LDA=fit.lda, RF=fit.rf))
summary(results)
ggplot(results) + labs(y = "Accuracy")

cm = confusionMatrix(fit.lda)

#

###### 10-fold cross validation ########
# --> split the data in 10 blocks, uses 9 for train and 1 for testing and repeat 10 times 

control <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
fit.lda.2 <- train(Individual~., data=dat, method="lda",
                   metric=metric, trControl=control)
fit.rf.2 <- train(Individual~., data=dat, method="rf",
                  metric=metric, trControl=control)

results <- resamples(list(LDA=fit.lda.2, RF=fit.rf.2))
ggplot(results) + labs(y = "Accuracy")


#### Normal and tumor differences

dim(ex2)
f2 <- factor(c(rep("Tumor",80), rep("Normal",30)))
library("genefilter")
ttest <- rowttests(ex2,f2); ttest; summary(ttest) # feature selection performing T-test filtering 
# on each row

keepers <- which(p.adjust(ttest$p.value)<0.0001)
ex4 <-ex2[keepers,]
dim(ex4)
tex4 <- t(ex4)
dat2 <- cbind(as.data.frame(tex4),f2)
colnames(dat2)[ncol(dat2)] <- "Individual"


library("caret")
library("e1071")

metric <- "Accuracy"
control <- trainControl(method="cv", number=10) # number = number of resampling iterations
fit.lda <- train(Individual~., data=dat2, method="lda", # y = t(w) * x
                 metric=metric, trControl=control)
fit.rf <- train(Individual~., data=dat2, method="rf",
                metric=metric, trControl=control)

results <- resamples(list(LDA=fit.lda, RF=fit.rf))
summary(results)
ggplot(results) + labs(y = "Accuracy")

cm = confusionMatrix(fit.lda)

###############################################

#Lasso

small.df = ex2[,1:80]

library("glmnet")
dim(small.df) #only 2 tissues taken into account
dat <- t(small.df)

y <- c(rep(0,40),rep(1,40))
f <- factor(y) #factor created before

fit=glmnet(dat,y,standardize=FALSE,family="binomial") #fits a lasso with various lambda
plot(fit, xvar = "lambda", label=TRUE)

# with cross validation
cfit=cv.glmnet(dat,y,standardize=FALSE,family="binomial", nfolds = 5) #performs cross validation and fits
plot(cfit) #shows the lamda that gives the least error (minimize the function)

cfit$lambda.min # lambda that gives the least error
cfit$lambda.1se

coef(cfit, s=cfit$lambda.min)

# repeat analysis but by using train + test sample subsets
smoker <- 40
non.smoker <- 40

train <- sample(1:(smoker), (non.smoker-5))
test <- setdiff(1:(smoker),train)

test<- c(test, test+40)
train <- c(train, train+40)

fit=glmnet(dat[train,],y[train],standardize=FALSE,family="binomial")
plot(fit)

cfit=cv.glmnet(dat[train,],y[train],standardize=FALSE,family="binomial")
plot(cfit)
predict(fit,dat[test,], type="class", s= cfit$lambda.min)

# plot ROCR curve
library("ROCR")
pred2 <- predict(fit,dat[test,], type="response", s=cfit$lambda.min)
plot(performance(prediction(pred2, y[test]), 'tpr', 'fpr'))
# compute Area Under the Curve (AUC)
auc.tmp <- performance(prediction(pred2, y[test]),"auc")
auc <- as.numeric(auc.tmp@y.values)

#########################################################################à
#Scudo

y <- c(rep(0,40),rep(1,40),rep(2,30))
f <- factor(y, labels = c("Smoker-tumor",
                          "Non-smoker-tumor", "Normal"))
library("caret")
set.seed(123)
inTrain <- createDataPartition(f, list = FALSE)

library("genefilter")
ttest <- rowttests(ex2,f); ttest; summary(ttest)
keepers <- which(p.adjust(ttest$p.value)<0.01)
ex3 <-small.df[keepers,]

trainData <- ex2[, inTrain]
testData <- ex2[, -inTrain]


# analyze training set
library("rScudo")
trainRes <- scudoTrain(trainData, groups = f[inTrain],
                       nTop = 25, nBottom = 25, alpha = 0.5)

trainRes

# inspect signatures
upSignatures(trainRes)[1:25,1:25]
consensusUpSignatures(trainRes)[1:25, ]

# generate and plot map of training samples
trainNet <- scudoNetwork(trainRes, N = 0.20)
scudoPlot(trainNet, vertex.label = NA)


# perform validation using testing samples
testRes <- scudoTest(trainRes, testData, f[-inTrain],
                     nTop = 50, nBottom = 50)
testNet <- scudoNetwork(testRes, N = 0.2)
scudoPlot(testNet, vertex.label = NA)

# identify clusters on map  
library("igraph")
testClust <- igraph::cluster_spinglass(testNet, spins = 2)
plot(testClust, testNet, vertex.label = NA)
# perform classification
classRes <- scudoClassify(trainData, testData, N = 0.25,
                          nTop = 50, nBottom = 50,
                          trainGroups = f[inTrain], alpha = 0.5)
caret::confusionMatrix(classRes$predicted, f[-inTrain])

install.packages("RCys3")
library("RCys3")
scudoCytoscape(testNet, title="Scudo GRaph", collection = "SCUDO")
get.edgelist(testNet)


renamed_200 = read.table("out_david.txt", sep = "\t")
renamed_200 = renamed_200[,2]

write(renamed_200, file = "file_name.txt", ncolumns = 1)
  

