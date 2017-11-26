# This is an example using CCSOGL for multinomial classification. 

# CCSOGL is able to select class-specific feature groups. For example, select pathways consisting of genes in
#> multi-class cancer classification. This type of feature selection is useful, since intuitively different
#> cancer types depend on different gene pathways. In this way, CCSOGL is of good interpretability.

# The data in this example is the classification of breast cancer subclasses. It can be downloaded from:
#>     http://portals.broadinstitute.org/cgi-bin/cancer/datasets.cgi

# The pathway information used is also accessible: 
#>     http://software.broadinstitute.org/gsea/msigdb/index.jsp

# We provide the preprossed data in the .RData file.
# First load the processed data for breast cancer subclass classification stored in the BCS1 C6.Rdata

require(caret)
fold = createFolds(BCS1.y, 4)

train.bcs1 = bcs1.x[-fold$Fold1,] # training data
test.bcs1 = bcs1.x[fold$Fold1,] # testing data
train.y = bcs1.y[-fold$Fold1,] # training data labels using one-hot encoding
dup.train = vector("list", length(BCS1set)) # feature duplication
for(i in 1:length(BCS1set)){
  dup.train[[i]] = train.bcs1[, bcs1.grpMembership[,i] == 1]
}

fit = cppMultiSOGL(train.bcs1, train.y, bcs1.grpMembership, dup.train, 0.001, 0.002)
plot(1:fit$iteration, fit$cost[1:fit$iteration], type = 'l', col = 'red')

pre = cbind(rep(1,length(fold$Fold1)), test.bcs1) %*% rbind(fit$beta0, fit$beta) 
pre = unname(apply(pre, 1, function(v){which(v == max(v))})) # prediction
accu = mean(pre == (BCS1.y[fold$Fold1] )) # prediction accuracy
accu

