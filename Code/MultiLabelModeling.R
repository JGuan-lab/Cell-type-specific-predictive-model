###Multilabel modeling
library(mlr)
library(parallelMap)


###Load data
traindata<-read.table('/Data/MultiLabelData/traindata.txt',sep=',')
trainlabel<-read.table('/Data/MultiLabelData/trainlabel.txt',sep=',')
testdata<-read.table('/Data/MultiLabelData/testdata.txt',sep=',')
testlabel<-read.table('/Data/MultiLabelData/testlabel.txt',sep=',')


###Training model.
training<-cbind(trainlabel,traindata)
testing<-cbind(testlabel,testdata)
labels = colnames(training)[1:18]
data.task = makeMultilabelTask(id = "multi", data = training, target = labels)
lrn.br = makeLearner("classif.plsdaCaret", predict.type = "prob")
lrn.br = makeMultilabelBinaryRelevanceWrapper(lrn.br)
parallelStartSocket(14)
parallelStop()
discrete_ps = makeParamSet(makeDiscreteParam("ncomp", values = seq(1:15)))
ctrl = makeTuneControlGrid()
rdesc = makeResampleDesc(method = "RepCV" ,reps=5,folds=5)
res = tuneParams(lrn.br, task = data.task, resampling = rdesc,
                 par.set = discrete_ps, control = ctrl)
lrn = setHyperPars(makeLearner("classif.plsdaCaret",predict.type = "prob"), ncomp = res$x$ncomp)
lrn = makeMultilabelBinaryRelevanceWrapper(lrn)
mod=mlr::train(lrn,data.task)


###Calculate model performance.
pred_train = predict(mod, task=data.task)
trainperformance<-performance(pred_train, measures = list(multilabel.subset01, multilabel.hamloss,multilabel.acc,multilabel.f1, timepredict))
trainget<-getMultilabelBinaryPerformances(pred_train, measures = list(acc, mmce, auc,tnr,tpr))
pred_test = predict(mod, newdata =testing)
test_performance<-performance(pred_test, measures = list(multilabel.subset01, multilabel.hamloss,multilabel.acc,multilabel.f1, timepredict))
test_get<-getMultilabelBinaryPerformances(pred_test, measures = list(acc, mmce, auc,tnr,tpr))







