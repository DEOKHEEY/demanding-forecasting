#################################################################################
library("depmixS4")
library("zoo")
library("TTR")
library("xts")
library("quantmod")
library("devEMF")
library("RHmm")
library("parallel")
library("caret")
library("mixtools")

Start.time<-Sys.time()

TSData_function<-function(dataset){
  dataset<-as.data.frame(dataset)
  Date<-as.character(dataset[,1])
  Data.date<-as.Date.character(Date,format = "%Y-%m-%d")
  TSData<-data.frame(dataset[,2],row.names = Data.date)
  colnames(TSData)<-c("SB50")
  TSData<-as.xts(TSData)
  return(TSData)
}
state_function<-function(dataset,dist){
  
  op.state.AIC<-0;op.state.BIC<-0;op.state.LLH<-0
  compare.AIC<-10000;compare.BIC<-10000;compare.Loglik<-10000
  
  for(i in 2:6){
    hmm_temp <- HMMFit(obs = dataset, nStates = i, dis = dist)
    param <- c(hmm_temp$AIC,hmm_temp$BIC,hmm_temp$LLH)
    print(param)
    print(i)
    if(is.nan(sum(param))){
      cat("state",i,"is not fit for this model","\n")
    }else{
      #AIC
      if(compare.AIC>hmm_temp$AIC){
        compare.AIC<-hmm_temp$AIC
        op.state.AIC<-i
      }
      #BIC
      if(compare.BIC>hmm_temp$BIC){
        compare.BIC<-hmm_temp$BIC
        op.state.BIC<-i
      }
      hmm_temp$LLH
      #Log-likelihood
      if(compare.Loglik>hmm_temp$LLH){
        compare.Loglik<-hmm_temp$LLH
        op.state.LLH<-i
      }
      
    }
    state.compare<-c(op.state.AIC,op.state.BIC,op.state.LLH)
    param.compare<-c(compare.AIC,compare.BIC,compare.Loglik)
    
    print(state.compare)
    print(param.compare)
    cat("\n")
  }
  
  #optimal state
  if(op.state.AIC==op.state.BIC && op.state.BIC==op.state.LLH){
    op.state<-state.compare[1]
  }else{
    cat("Not Number of Optimal States \n")
    #  stop()
    op.state<-2
  }
  cat("optimal state : ", op.state,"\n")
  op.state<-2
  return(op.state)
  
}

#Data split - training set & test set
#Data Handling
setwd("~/Desktop/R/Working Directory")
Dataset<-read.csv("2011-2015.csv",header = F,fileEncoding = "euc-kr")
options(warn=1)

Type<-subset(Dataset, select = c(V1,V20))
Type<-Type[-1,]
Name<-c("Date","SB50")
colnames(Type)<-Name
Type$SB50<-as.numeric(as.character(Type$SB50))
Type<-as.data.frame(Type)

#Time serise
Date<-as.character(Type[,1])
Data.date<-as.Date.character(Date,format = "%Y.%m.%d")
TSData.train<-data.frame(Type[,2],row.names = Data.date)
colnames(TSData.train)<-c("SB50")
TSData.train<-as.xts(TSData.train)
train<-TSData.train

################################## Over-sampling #####################################
train.SB50<-train$SB50

kmeans<-kmeans(train.SB50, 2)
kmeans.cluster<-kmeans$cluster
clustering.data<-data.frame(x = train.SB50, cluster = kmeans.cluster)
clust_1<-subset(clustering.data,clustering.data$cluster==1)
clust_1<-clust_1[-2]
clust_2<-subset(clustering.data,clustering.data$cluster==2)
clust_2<-clust_2[-2]
param_1<-c(kmeans$size[1],mean(clust_1$SB50),sd(clust_1$SB50))
param_2<-c(kmeans$size[2],mean(clust_2$SB50),sd(clust_2$SB50))

#oversampling size
oversamp.size<-10

samp.1<-data.frame(rnorm(param_1[1]*oversamp.size,param_1[2],param_1[3]))
samp.2<-data.frame(rnorm(param_1[1]*oversamp.size,param_1[2],param_1[3]))
sampled.train<-rbind(samp.1,samp.2)
sampled.train<-sampled.train[which(sampled.train>=0),]

################################################################################# 
########################### Classified Sales data ###############################
#################################################################################
train.sales<-data.frame(matrix(nrow=length(sampled.train),ncol=2))
for(i in 1:length(sampled.train)){
  if(sampled.train[i]>=0.5){
    train.sales[i,2]<-1
  }else
    train.sales[i,2]<-0
}
train.sales<-train.sales$X2

#Testset
Test.Dataset.Sales<-read.csv("2016sales.csv",header = F,fileEncoding = "euc-kr")

##Test set
Sales.data.2016<-subset(Test.Dataset.Sales, select = c(V1,V20))
Sales.data.2016<-Sales.data.2016[-1,]
Name<-c("Date","SB50")
colnames(Sales.data.2016)<-Name
Sales.data.2016$SB50<-as.numeric(as.character(Sales.data.2016$SB50))
Sales.data.2016<-as.data.frame(Sales.data.2016)

#Time serise
test.sales<-as.character(TSData_function(Sales.data.2016))

#Define states
op.state.sales<-state_function(train.sales,"DISCRETE")

df.sales<-data.frame(matrix(nrow=15,ncol=2))
for(i in 2:15){
  hmm_temp<-HMMFit(obs = train.sales, nStates = i, dis = "DISCRETE")
  df.sales[i,1]<-c(i)
  df.sales[i,2]<-c(hmm_temp$LLH)
}
df.sales
op.state.sales

#find the model for the given oberv
hm_model.sales <- HMMFit(obs = train.sales, nStates = op.state.sales, dis = "DISCRETE",
                         Levels = c("0","1"))
hm_model.sales

#Viterbi
VitPath.Sales<-viterbi(hm_model.sales, train.sales)
VitPath.Sales

# add a new colum "Pred"
testset<-cbind(test.sales, Pred = 0)

# number of rows of test set data
rows = nrow(testset)
train.sales
classified<-as.data.frame.numeric(test.sales)

################################## Classified #####################################
sales.dist_set<-distributionSet(dis="DISCRETE",hm_model.sales$HMM$distribution$proba,
                                labels=NULL)

hm_sales_set<-HMMSet(hm_model.sales$HMM$initProb,
                     hm_model.sales$HMM$transMat,
                     sales.dist_set)
simul.sales<-HMMSim(length(test.sales),hm_sales_set)
simul.sales
simul.sales <- c(simul.sales, HMMSim(length(test.sales), hm_sales_set, 
                                     simul.sales$states[length(test.sales)]))
simul.sales$obs

classified.samp<-data.frame()
for(i in 1:length(test.sales)){
  if(simul.sales$obs[i]=="p 1"){
    classified.samp[i,1]<-"0"
  }else
    classified.samp[i,1]<-"1"
}
test.actual<-as.data.frame.numeric(test.sales)

colnames(classified.samp)<-"Sampling"
colnames(test.actual)<-"Actual"
classified<-data.frame(test.actual,classified.samp)

#Misclassification rate
count<-0
for(k in 1:length(test.sales)){
  if(classified[k,1]==classified[k,2])
    count<-count+1
}

classification.rate<-count/length(test.sales)
cat("classification rate =",classification.rate,"\n")
write.csv(classification.rate,file="/Users/Deokhee/Desktop/misclassification.rate_SB50.csv")  

as.matrix(classified)
xtab <- table(classified$Sampling, classified$Actual)
write.csv(xtab,file="/Users/Deokhee/Desktop/misclassification.table_SB50.csv") 


################################################################################# 
################################## Prediction ###################################
#################################################################################
#Number of Mixture gaussian distribution
library(mclust)
xy<-Mclust(sampled.train)
xy
op.nc<-xy$G
op.nc<-2

#Number of Hidden markov states
df.pred<-data.frame(matrix(nrow=15,ncol=2))
for(i in 2:10){
  hmm_temp<-HMMFit(obs = sampled.train, nStates = i, nMixt = op.nc, dis = "MIXTURE")
  df.pred[i,1]<-c(i)
  df.pred[i,2]<-c(hmm_temp$LLH)
}
df.pred
op.state

#find the model for the given oberv
train<-sampled.train
hm_model <- HMMFit(obs = train, nStates = op.state, nMixt = op.nc, dis = "MIXTURE",
                   control=list(tol=1e-10),
                   asymptCov=TRUE)
hm_model

#Viterbi
VitPath.RHmm<-viterbi(hm_model, train)
VitPath.RHmm

#Testset
Test.dataset<-read.csv("2016.csv",header = F,fileEncoding = "euc-kr")
test<-subset(Test.dataset, select = c(V1,V20))
test<-test[-1,]
colnames(test)<-Name
test$SB50<-as.numeric(as.character(test$SB50))
test<-as.data.frame(test)
test

#Timeserise
Date.test<-as.character(test[,1])
Data.date.test<-as.Date.character(Date.test,format = "%Y.%m.%d")
TSData.test<-data.frame(test[,2],row.names = Data.date.test)
colnames(TSData.test)<-c("SB50")
TSData.test<-as.xts(TSData.test)
test<-TSData.test

# add a new colum "Pred"
testset<-cbind(test, Pred = 0)

# number of rows of test set data
rows = nrow(testset)

#predict obs (Criteria Single set)
for (i in 1: rows) {
  
  if(i != 0) {
    testrow <- testset[i, ]
  }
  
  
  if(classified[i,2]==0){
    change.RHmm<-0
  }else{
    # predict the closing value of today
    change.RHmm <- sum(hm_model$HMM$transMat[last(VitPath.RHmm$states),] *
                         .colSums((matrix(unlist(hm_model$HMM$distribution$mean), nrow=op.nc,ncol=op.state)) *
                                    (matrix(unlist(hm_model$HMM$distribution$proportion), nrow=op.nc,ncol=op.state)), 
                                  m=op.nc,n=op.state))
  }
  
  print(hm_model$HMM$distribution)
  
  pred.RHmm<-testrow$Pred + change.RHmm
  
  
  # update predicted value
  testset[i, ]$Pred <- pred.RHmm
  print(testset[i,])
  temp<-as.numeric(testrow$SB50)
  temp<-as.data.frame(temp)
  train<-as.data.frame(train)
  colnames(temp)<-"train"
  temp
  train
  train<-rbind(train, temp)
  
  if(classified[i,1]==1){
    # update HMM with the new data
    # Baum-Welch Algorithm to find the model for the given observations
    hm_model <- HMMFit(obs = train, nStates = op.state, nMixt = op.nc, dis = "MIXTURE")
    
    # Viterbi Algorithm to find the most probable state sequence
    VitPath.RHmm <- viterbi (hm_model, train)
    
  }
}

write.csv(testset,file="/Users/Deokhee/Desktop/capstone.csv")  
End.Time<-Sys.time();print(Start.time);print(End.Time)

