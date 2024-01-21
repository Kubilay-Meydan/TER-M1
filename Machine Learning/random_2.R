#Classification par forêts aléatoires
library(randomForest)
library(dplyr)
source("sample-equal.R")

#On crée le dataset 'équilibré'
df = create_sample_equal_dataset(read.csv(file = './data/data.csv'), c(SS = 0.16, S = 0.16, VS = 0.16, D = 0.52))


df$Annotation = factor(df$Annotation)
df = subset(df, select = -c(ProcessA, ProcessB))
#On entraine le modèles
fit.rf=randomForest(Annotation~., data=df, ntree=100)
feature_importance <- importance(fit.rf)
#On affiche les attributs avec le plus de poids
print(feature_importance[order(-feature_importance[, "MeanDecreaseGini"]),])
print(0)





