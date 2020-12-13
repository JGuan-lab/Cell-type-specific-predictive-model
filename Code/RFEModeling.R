###Use RFE to find important features and model them.
library(pROC)
library(caret)
name <- list('AST-FB','AST-PP','Endothelial','IN-PV','IN-SST','IN-SV2C','IN-VIP','L2_3','L4','L5_6'
             ,'L5_6-CC','Microglia','Neu-mat','Neu-NRGN-I','Neu-NRGN-II','Oligodendrocytes','OPC')
performance<-data.frame()
for (i in 1:length(name))
{
  ###Load data.
  print(name[i])
  url_training <- paste0('/Data/CellTypeExpressionData/',name[i],'/','data','/','training.txt')
  url_testing <- paste0('/Data/CellTypeExpressionData/',name[i],'/','data','/','testing.txt')
  training<-read.table(url_training,sep=',')
  testing<-read.table(url_testing,sep=',')
  
  
  ###Use RFE to find important features
  caretFuncs$summary<-twoClassSummary
  ctrl <- trainControl(
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
  )
  rfectrl <- rfeControl(
    functions=caretFuncs,
    method = "repeatedcv",
    repeats = 10,
    number = 10,
    verbose = TRUE,
    rerank = FALSE
  )
  k<-seq(100,1000,100)
  set.seed(1)
  ldaProfile <- rfe(
    group ~ .,
    data = training,
    method = "pls",
    tuneLength = 15,
    trControl = ctrl,
    rfeControl = rfectrl,
    metric = "ROC",
    sizes=k,
  )
  
  
  ###Use important features to model and calculate model performance.
  if(ldaProfile$bestSubset==12036)
  {
    li<-ldaProfile$results
    li<-li[1:10,]
    opt<-li$Variables[which.max(li$ROC)]
    optvar<-update(ldaProfile,x=training,y=training$group,size = opt)
    optvar<-optvar$bestVar
    optvar<-c(optvar,'group')
    training<-training[,optvar]
    testing<-testing[,optvar]
    train_group<-training$group
    test_group<-testing$group
    ctrl <- trainControl(
      method = "repeatedcv",
      repeats = 10,
      number = 10,
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
    )
    set.seed(1)
    plsFit <- train(
      group ~ .,
      data = training,
      method = "pls",
      tuneLength = 15,
      trControl = ctrl,
      metric = "ROC",
    )
    train_prob <- predict(plsFit, newdata = training,type='prob')
    train_prob <- as.numeric(train_prob$ASD)
    test_prob <- predict(plsFit, newdata = testing,type='prob')
    test_prob <- as.numeric(test_prob$ASD)
    roc_train<- roc(train_group, train_prob, percent = FALSE,levels=c("Control", "ASD"))
    th<-coords(roc_train, "best", ret="threshold", transpose = FALSE)
    prob<-th$threshold
    roc_test <- roc(test_group, test_prob, percent = FALSE,levels=c("Control", "ASD"))
    for (p in 1:length(train_group))
    {
      if (train_prob[p]>prob)
      {train_prob[p]<-'ASD'}
      else
      {train_prob[p]<-'Control'}
    }
    train_matrix<-confusionMatrix(data = as.factor(train_prob), train_group)
    for (p in 1:length(test_group))
    {
      if (test_prob[p]>prob)
      {test_prob[p]<-'ASD'}
      else
      {test_prob[p]<-'Control'}
    }
    test_matrix<-confusionMatrix(data = as.factor(test_prob), test_group)
    row<-data.frame(cellcluster=as.character(name[i])
                    ,Accuracytrain=train_matrix$overall[1],Sensitivitytrain=train_matrix$byClass[1]
                    ,Specifitytrain=train_matrix$byClass[2]
                    ,AUCtrain=roc_train$auc,Cutoff=th
                    ,Accuracytest=test_matrix$overall[1]
                    ,Sensitivitytest=test_matrix$byClass[1],Specifitytest=test_matrix$byClass[2],AUCtest=roc_test$auc,optimal=opt,class='rfe&pls',row.names = NULL)
    performance<-rbind(performance,row)
  }
  else
  {
    train_group<-training$group
    test_group<-testing$group
    train_prob <- predict(ldaProfile, newdata = training,type='prob')
    train_prob <- as.numeric(train_prob$ASD)
    test_prob <- predict(ldaProfile, newdata = testing,type='prob')
    test_prob <- as.numeric(test_prob$ASD)
    roc_train<- roc(train_group, train_prob, percent = FALSE,levels=c("Control", "ASD"))
    th<-coords(roc_train, "best", ret="threshold", transpose = FALSE)
    prob<-th$threshold
    roc_test <- roc(test_group, test_prob, percent = FALSE,levels=c("Control", "ASD"))
    for (p in 1:length(train_group))
    {
      if (train_prob[p]>prob)
      {train_prob[p]<-'ASD'}
      else
      {train_prob[p]<-'Control'}
    }
    train_matrix<-confusionMatrix(data = as.factor(train_prob), train_group)
    for (p in 1:length(test_group))
    {
      if (test_prob[p]>prob)
      {test_prob[p]<-'ASD'}
      else
      {test_prob[p]<-'Control'}
    }
    test_matrix<-confusionMatrix(data = as.factor(test_prob), test_group)
    row<-data.frame(cellcluster=as.character(name[i])
                    ,Accuracytrain=train_matrix$overall[1],Sensitivitytrain=train_matrix$byClass[1]
                    ,Specifitytrain=train_matrix$byClass[2]
                    ,AUCtrain=roc_train$auc,Cutoff=th
                    ,Accuracytest=test_matrix$overall[1]
                    ,Sensitivitytest=test_matrix$byClass[1],Specifitytest=test_matrix$byClass[2],AUCtest=roc_test$auc,optimal=ldaProfile$bestSubset,class='rfe',row.names = NULL)
    performance<-rbind(performance,row)
  }
}