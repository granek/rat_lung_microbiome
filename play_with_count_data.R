# https://www.r-bloggers.com/do-not-log-transform-count-data-bitches/

# QUASI-POISSON VS. NEGATIVE BINOMIAL REGRESSION: HOW SHOULD WE MODEL OVERDISPERSED COUNT DATA? http://dx.doi.org/10.1890/07-0043.1
# http://stackoverflow.com/questions/21016819/anova-with-count-data-glm
# Generalized linear mixed models: a practical guide for ecology and evolution http://dx.doi.org/10.1016/j.tree.2008.10.008

library(multcomp)
### set up all pair-wise comparisons for count data
data(Titanic)
mod <- glm(Survived ~ Class, data = as.data.frame(Titanic), weights = Freq, family = binomial)
x = as.data.frame(Titanic)
### specify all pair-wise comparisons among levels of variable "Class"
### Note, Tukey means the type of contrast matrix.  See ?contrMat
glht.mod <- glht(mod, mcp(Class = "Tukey"))

###summaryize information
###applying the false discovery rate adjustment
###you know, if that's how you roll
summary(glht.mod, test=adjusted("fdr"))

