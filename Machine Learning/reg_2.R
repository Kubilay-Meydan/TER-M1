library(dplyr)
source("sample-equal.R")

#On crée le dataset 'équilibré'
df = create_sample_equal_dataset(read.csv(file = './data/data.csv'), c(SS = 0.16, S = 0.16, VS = 0.16, D = 0.52))

#Simplifying the annotations (to only have binary values similar/not similar)
#Pour la classification linéraire (logistique) -> on peut seulement 2 valeurs binaires
df$Annotation = df$Annotation!='D'

df = subset(df, select = -c(ProcessA, ProcessB))
alpha = 0.5
#On entraine le modèle
model.train=glm(Annotation~., data=df, family=binomial)
#On affiche les résultats avec le plus de poids
print(coef(summary(model.train))[,4])
cat(0)

