library(rpart)
library(dplyr)
source("sample-equal.R")

#On crée le dataset 'équilibré'
df = create_sample_equal_dataset(read.csv(file = './data/data.csv'), c(SS = 0.25, S = 0.25, VS = 0.25, D = 0.25))

df$Annotation = factor(df$Annotation)
df = subset(df, select = -c(ProcessA, ProcessB))

#On crée l'abre
tmax=rpart(Annotation~., data=df, control=rpart.control(cp=0, minsplit=1))
cart = rpart(Annotation~., data=df, control=rpart.control(cp=tmax$cptable[which.min(tmax$cptable[,4]), 1]))

#On affiche l'arbre
print(cart)
print(0)





