library(corrplot)
library(Hmisc)
df <- read.csv("D:/StatisticalAnalysis.csv", header = T, row.names = 1)
df <- df[-1]
cor <- cor(df)
cor_test_mat <- rcorr(as.matrix(df))$P
corrplot(cor, p.mat = cor_test_mat, insig = "p-value")
corrplot(cor, p.mat = cor_test_mat, sig.level = 0.01, insig = "blank")
corrplot(cor, p.mat = cor_test_mat, insig = "label_sig", sig.level = c(.001, .01, .05),
         pch.cex=1.4, tl.cex = .9, number.font=2, number.cex=0.8)

citation('corrplot')
citation()
