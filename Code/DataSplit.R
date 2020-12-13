###Divide the data into training and test sets at the ratio of 7:3.
library(caret)
name <- list('AST-FB','AST-PP','Endothelial','IN-PV','IN-SST','IN-SV2C','IN-VIP','L2_3','L4','L5_6'
             ,'L5_6-CC','Microglia','Neu-mat','Neu-NRGN-I','Neu-NRGN-II','Oligodendrocytes','OPC')
for (i in 1:length(name))
{
  ###Load data and create path.
  print(name[i])
  url_coldata <- paste0('/Data/CellTypeExpressionData/',name[i],'/',name[i],'_cellcluster_supplyinfo_hvg.txt')
  url_logdata<-paste0('/Data/CellTypeExpressionData/',name[i],'/',name[i],'_cellcluster_logc_hvg.txt')
  url_training <- paste0('/Data/CellTypeExpressionData/',name[i],'/','data','/','training.txt')
  url_testing <- paste0('/Data/CellTypeExpressionData/',name[i],'/','data','/','testing.txt')
  logcount <-read.table(file=url_logdata,sep="",header=T,check.names = T,row.names = 1)
  col <- read.table(file=url_coldata,sep="",header=T,check.names = T,row.names = 1)
  
  
  ###Divide the data.
  group<-col[,7]
  logcount<-t(logcount)
  logcount<-data.frame(logcount)
  logcount$group<-group
  set.seed(1)
  inTrain <- createDataPartition(
    y = logcount$group,
    ## the outcome data are needed
    p = .7,
    ## The percentage of data in the
    ## training set
    list = FALSE
  )
  training <- logcount[ inTrain,]
  testing  <- logcount[-inTrain,]
  write.table(training,file = url_training,sep=',')
  write.table(testing,file = url_testing,sep=',')
}