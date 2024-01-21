library(janitor)
library(rpart)
library(randomForest)
library(foreign)
library(ggplot2)
library(MASS)
library(Hmisc)
library(reshape2)
library(dplyr)
source("sample-equal.R")

args <- commandArgs(trailingOnly = TRUE)
arg <- args[1]

columns<-c(args)

formule = reformulate(termlabels = columns, response='Annotation')

#On crée le dataset 'équilibré'
df = create_sample_equal_dataset(read.csv(file = './data/data.csv'), c(SS = 0.16, S = 0.16, VS = 0.16, D = 0.52))

df$Annotation = df$Annotation!='D'
df = subset(df, select = -c(ProcessA, ProcessB))
test = sample(1:length(df$Annotation),as.integer(length(df$Annotation)*0.3))
train = -test
train = df[train, ]
test = df[test,]

get_sensitivity <- function(input_tab, sensivity_of) {
    df_sensitivity  = as.data.frame.matrix(input_tab)
    diff = grep(sensivity_of, colnames(df_sensitivity))[1]
    if(sensivity_of=="D"){
        tp = df_sensitivity$D[diff]
        fn = df_sensitivity$S[diff]+df_sensitivity$C[diff]+df_sensitivity$VS[diff]
    }
    if(sensivity_of=="SS"){
        tp = df_sensitivity$C[diff]
        fn = df_sensitivity$S[diff]+df_sensitivity$D[diff]+df_sensitivity$VS[diff]
    }
    if(sensivity_of=="S"){
        if(diff==grep("VS", colnames(df_sensitivity))[1]){
            diff=grep("S", colnames(df_sensitivity))[2]
        }
        tp = df_sensitivity$S[diff]
        fn = df_sensitivity$C[diff]+df_sensitivity$D[diff]+df_sensitivity$VS[diff]
    }
    if(sensivity_of=="VS"){
        tp = df_sensitivity$VS[diff]
        fn = df_sensitivity$S[diff]+df_sensitivity$C[diff]+df_sensitivity$D[diff]
    }
    sensivity = tp/(tp+fn)
    sensivity
}


#REG
alpha = 0.5
model.train=glm(formula=formule, data=train, family=binomial)

pred.glm=rep(FALSE,length(test$Annotation))
options(warn=-1)      #turn off warnings
probs.glm = predict(model.train, test, type="response")
options(warn=1)      #turn warnings back on
pred.glm[probs.glm>alpha]= TRUE # classifieur de Bayes

mat = as.data.frame.matrix(table(test$Annotation,pred.glm))  
tp = mat$"TRUE"[2]
tn = mat$"FALSE"[1]
fn = mat$"FALSE"[2]
fp = mat$"TRUE"[1]

sensivity = tp/(tp+fn)

cat("reg ER ", (tp+fp)/(length(pred.glm)), "\n")
cat("reg sensitivity_binary ", sensivity, "\n")


#==========================
#CART
#==========================
#On crée le dataset 'équilibré'
df = create_sample_equal_dataset(read.csv(file = './data/data.csv'), c(SS = 0.25, S = 0.25, VS = 0.25, D = 0.25))


df$Annotation = factor(df$Annotation)
df = subset(df, select = -c(ProcessA, ProcessB))


test = sample(1:length(df$Annotation),as.integer(length(df$Annotation)*0.3))
train = -test
train = df[train, ]
test = df[test,]


tmax=rpart(formula=formule, data=train, control=rpart.control(cp=0, minsplit=1))
cart = rpart(formula=formule, data=train, control=rpart.control(cp=tmax$cptable[which.min(tmax$cptable[,4]), 1]))
pred.cart = predict(cart, newdata=test, type="class")


#mat = as.data.frame.matrix(table(test$Annotation,pred.cart)) 
 
sensitivity_D = get_sensitivity(table(test$Annotation,pred.cart), "D")
sensitivity_SS = get_sensitivity(table(test$Annotation,pred.cart), "SS")
sensitivity_S = get_sensitivity(table(test$Annotation,pred.cart), "S")
sensitivity_VS = get_sensitivity(table(test$Annotation,pred.cart), "VS")

cat("cart ER ", 1-(sum(diag(table(test$Annotation,pred.cart)))/length(pred.cart)), "\n")
cat("cart sensitivity_D ", sensitivity_D, "\n")
cat("cart sensitivity_C ", sensitivity_SS, "\n")
cat("cart sensitivity_S ", sensitivity_S, "\n")
cat("cart sensitivity_VS ", sensitivity_VS, "\n")


#==========================
#Random Forest
#==========================
fit.rf=randomForest(formule, data=train, ntree=100)

pred.rf = predict(fit.rf, newdata=test, type="class")
tab.rf = table(test$Annotation,pred.rf)

sensitivity_D = get_sensitivity( tab.rf, "D")
sensitivity_SS = get_sensitivity( tab.rf, "SS")
sensitivity_S = get_sensitivity( tab.rf, "S")
sensitivity_VS = get_sensitivity(tab.rf, "VS")

cat("random ER ", 1-(sum(diag(table(test$Annotation,pred.rf)))/length(pred.rf)), "\n")
cat("random sensitivity_D ", sensitivity_D, "\n")
cat("random sensitivity_C ", sensitivity_SS, "\n")
cat("random sensitivity_S ", sensitivity_S, "\n")
cat("random sensitivity_VS ", sensitivity_VS, "\n")







