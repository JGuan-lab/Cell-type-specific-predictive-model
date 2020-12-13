###Use gene set modeling.
library(pROC)
library(caret)
load('/Data/OtherData/genesetname.RData')
geneinfo <-read.table('/Data/OtherData/asd_sce_hvg_12036\ geneinfo.txt',header = TRUE)
name <- list('AST-FB','AST-PP','Endothelial','IN-PV','IN-SST','IN-SV2C','IN-VIP','L2_3','L4','L5_6'
             ,'L5_6-CC','Microglia','Neu-mat','Neu-NRGN-I','Neu-NRGN-II','Oligodendrocytes','OPC')
performance<-data.frame()
for (i in 1:length(name))
{
  ###Load data and create path.
  print(name[i])
  url_training <- paste0('/Data/CellTypeExpressionData/',name[i],'/','data','/','training.txt')
  url_testing <- paste0('/Data/CellTypeExpressionData/',name[i],'/','data','/','testing.txt')
  url_dir <- paste0('/SaveModel',name[i],'/','geneset')
  url_dirmodel <- paste0('/SaveModel',name[i],'/','geneset','/','model')
  url_diracc <- paste0('/SaveModel',name[i],'/','geneset','/','acc')
  url_dirtable <- paste0('/SaveModel',name[i],'/','geneset','/','table')
  url_excel<-paste0('/SaveModel',name[i],'/','geneset/acc','/','excel.txt')
  excel<-data.frame()
  dir.create(url_dir)
  dir.create(url_dirmodel)
  dir.create(url_diracc)
  dir.create(url_dirtable)
  training<-read.table(url_training,sep=',')
  testing<-read.table(url_testing,sep=',')
  
  
  ###Training model and calculate model performance.
  train_group<-training$group
  test_group<-testing$group
  training<-training[,1:12036]
  testing<-testing[,1:12036]
  names(training)<-geneinfo$symbol
  names(testing)<-geneinfo$symbol
  training$group<-train_group
  testing$group<-test_group
  info<-colnames(testing)
  for (j in 1:nrow(set))
  {
    m<-set[j,2:1999]
    n<-t(m)
    n<-as.factor(n)
    li<-intersect(info,n)
    jiao<-li
    print(length(li))
    if(length(li)>30)
    {
      print(set[j,1])
      li<-c(li,'group')
      trainingset<-training[,li]
      testingset<-testing[,li]
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
        data = trainingset,
        method = "pls",
        tuneLength = 15,
        trControl = ctrl,
        metric = "ROC",
      )
      train_prob <- predict(plsFit, newdata = trainingset,type='prob')
      train_prob <- as.numeric(train_prob$ASD)
      test_prob <- predict(plsFit, newdata = testingset,type='prob')
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
      if(test_matrix$overall[1]>0.7)
      {
        url_model <- paste0('/SaveModel/',name[i],'/','geneset/','model','/',set[j,1],'.RData')
        save(plsFit, file = url_model)
        imp<-varImp(plsFit)
        imp<-imp$importance
        imp<-data.frame(name=rownames(imp),Overall=imp$Overall)
        imp<-imp[order(imp$Overall,decreasing = TRUE),]
        imp<-imp[1:5,]
        imp<-imp$name
        gene<-set[j,]
        gene[gene==""]<-NA
        gene<-gene[,-1]
        n<-t(gene)
        n<-as.factor(n)
        w<-0
        for( q in 1:length(n))
        {
          if(is.na(n[q])!=TRUE)
          {
            w<-w+1;
          }
        }
        row<-data.frame(cellcluster=as.character(name[i]),geneset=set[j,1],interect=paste0(length(jiao),'/',w),TOP1=imp[1],TOP2=imp[2]
                        ,TOP3=imp[3],TOP4=imp[4],TOP5=imp[5]
                        ,Accuracytrain=train_matrix$overall[1],Sensitivitytrain=train_matrix$byClass[1]
                        ,Specifitytrain=train_matrix$byClass[2]
                        ,AUCtrain=roc_train$auc,Cutoff=th
                        ,Accuracytest=test_matrix$overall[1]
                        ,Sensitivitytest=test_matrix$byClass[1],Specifitytest=test_matrix$byClass[2],AUCtest=roc_test$auc,row.names = NULL)
        performance<-rbind(performance,row)
      }
    }
  }
  write.table(performance,url_excel,sep=',')
}

