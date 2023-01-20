# title: "Script Longueur"
# author: "Juju"
# date: "03/08/2021"

rm(list= ls())
library(ggplot2) 
library(lme4)
library(lattice)
library(reshape2)
library(gridExtra)
library(dfoptim)
library(optimx)
library(nloptr)
library(knitr)
library(MASS)
library(car)

## Stratégie analytique : 
#
#  Y =        Longueur
#  Espèces =  Homard et Tourteau
#  Mois =     Juin et Septembre

## IV.      Y = Longueur
# Chargement du jeu de donnees :
## Definition de la racine : 
#racine   <- "C:/Users/jnormand/Documents/general/boulot/homards/Etude_complementaire"
setwd("N:/02_OBSERVATION_SURVEILLANCE/05_IGA/20_CNPE_FLAMANVILLE/56=ETUDES COMPLEMENTAIRES/2019_CANTOFLAM")

## je charge le jeu de donn?es complet des CPUE : toutes les esp?ces, tous les mois
TAB <- read.table("data/data_analyses_explo/long_et_masse/allspecies_allmonths_long_et_masse.csv", 
                  header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
TAB                                        <- droplevels(TAB)
tab2  <- read.table("data/data_analyses_explo/dist.csv", header=TRUE, sep=";", na.strings="NA", dec=".", strip.white=TRUE)

tab2$POINT <- as.factor(tab2$POINT)

TAB$ZONE     <- as.factor(TAB$DR)
TAB$PERIODE  <- as.factor(TAB$APR)
TAB$CAMP     <- as.factor(TAB$AN)
TAB$POINT    <- as.factor(TAB$POINT)
TAB <- merge(x= TAB, y= tab2, by= "POINT")

#TAB[which(TAB$ZONE == "DEDANS"), "dist"] <- 0
moy  <- mean(TAB$dist, na.rm= T)
et   <- sd(TAB$dist, na.rm= T)
TAB$scale_dist  <- (TAB$dist- moy)/et
TAB$scale_annee <- (TAB$AN - mean(TAB$AN)) / sd(TAB$AN)

summary(TAB)
TAB <- TAB[, c("PERIODE", "CAMP", "ZONE", "POINT", "MOIS", "ESPECE", "TAILLE", "POIDS", "scale_dist")]
rm(tab2, et, moy)

colnames(TAB)[5:8] <- c("mois", "espece", "long", "masse")
TAB$espece <- as.factor(TAB$espece)
TAB$mois <- as.factor(TAB$mois)
summary(TAB)
######################################################################
#####////////////////       HOMARD     \\\\\\\\\\\\\\\\\\\\\\\\\\#####
######################################################################
i=3
tab <- droplevels(TAB[which(TAB$espece == levels(TAB$espece)[i]),])
print(levels(TAB$espece)[i])

######################## Mois= Juin ##################################
j=1
tab1 <- droplevels(tab[which(tab$mois == levels(tab$mois)[j]),])
print(levels(tab$mois)[j])

# Formula is: Longueur  ~ ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP)
lmm1 <- lmer(long ~ ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP), data = tab1)
#Estimation de la variance des intercepts ZONE:PERIODE par CAMP et par POINT, GLMM_1
print(summary(lmm1)) ## Ici, on r?cup?re la variance et la Std Dev. pour les effets random ET les Estimate et Std Error associ?s aux effets fixes
# test des effets random
lmm1_2 <- update(lmm1, . ~ . - (1 | POINT))
lmm1_3 <- update(lmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(lmm1, lmm1_2)
test_CAMP  <- anova(lmm1, lmm1_3)
print(test_POINT) ## Ici, test de l'effet Random POINT, on r?cup?re Df (=6/7), Chisq et Pr(>Chisq)
print(test_CAMP)
# test des effets fixes
print(Anova(lmm1))  ## Ici, on r?cup?re Chisq, Df et Pr(>Chsiq)
rm(lmm1_2, lmm1_3)


#Plot estimates (ranef) and standard-error
newdat                    <- expand.grid(ZONE = unique(tab1$ZONE), PERIODE = unique(tab1$PERIODE))
newdat$long               <- predict(lmm1, newdat, re.form= NA, type = "response")
mm                        <- model.matrix(terms(lmm1), newdat)
pvar1                     <- diag(mm %*% tcrossprod(vcov(lmm1),mm))                   
tvar1                     <- pvar1 + as.vector(VarCorr(lmm1)$CAMP) + as.vector(VarCorr(lmm1)$POINT)
cmult                     <- 1.96

newdat  <- data.frame(
  newdat
  , plo = newdat$long-cmult*sqrt(pvar1)
  , phi = newdat$long+cmult*sqrt(pvar1)
  , tlo = newdat$long-cmult*sqrt(tvar1)
  , thi = newdat$long+cmult*sqrt(tvar1)
  )
newdat$PERIODE<-factor(newdat$PERIODE,c("AVANT","APRES"))
newdat$sp <- "lobster"
newdat$month <- "june"
stock <- newdat
######################## Mois= Sept ##################################
j=2
tab1 <- droplevels(tab[which(tab$mois == levels(tab$mois)[j]),])
print(levels(tab$mois)[j])

#Aa1. Modèle LMM_1

# Formula is: Longueur  ~ ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP)
lmm1 <- lmer(long ~ ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP), data = tab1)
#Estimation de la variance des intercepts ZONE:PERIODE par CAMP et par POINT, GLMM_1
print(summary(lmm1)) ## Ici, on r?cup?re la variance et la Std Dev. pour les effets random ET les Estimate et Std Error associ?s aux effets fixes
# test des effets random
lmm1_2 <- update(lmm1, . ~ . - (1 | POINT))
lmm1_3 <- update(lmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(lmm1, lmm1_2)
test_CAMP  <- anova(lmm1, lmm1_3)
print(test_POINT) ## Ici, test de l'effet Random POINT, on r?cup?re Df (=6/7), Chisq et Pr(>Chisq)
print(test_CAMP)
# test des effets fixes
print(Anova(lmm1))  ## Ici, on r?cup?re Chisq, Df et Pr(>Chsiq)


#Plot estimates (ranef) and standard-error
newdat                    <- expand.grid(ZONE = unique(tab1$ZONE), PERIODE = unique(tab1$PERIODE))
newdat$long               <- predict(lmm1, newdat, re.form= NA)
mm                        <- model.matrix(terms(lmm1), newdat)
pvar1                     <- diag(mm %*% tcrossprod(vcov(lmm1),mm))                   
tvar1                     <- pvar1 + as.vector(VarCorr(lmm1)$CAMP) + as.vector(VarCorr(lmm1)$POINT)
cmult                     <- 1.96

newdat  <- data.frame(
  newdat
  , plo = newdat$long-cmult*sqrt(pvar1)
  , phi = newdat$long+cmult*sqrt(pvar1)
  , tlo = newdat$long-cmult*sqrt(tvar1)
  , thi = newdat$long+cmult*sqrt(tvar1)
  )
newdat$PERIODE<-factor(newdat$PERIODE,c("AVANT","APRES"))
newdat$sp <- "lobster"
newdat$month <- "sept."
stock <- rbind(newdat, stock)
######################################################################
#####////////////////     TOURTEAU     \\\\\\\\\\\\\\\\\\\\\\\\\\#####
######################################################################
i=4
tab <- droplevels(TAB[which(TAB$espece == levels(TAB$espece)[i]),])
print(levels(TAB$espece)[i])

######################## Mois= Juin ##################################
j=1
tab1 <- droplevels(tab[which(tab$mois == levels(tab$mois)[j]),])
print(levels(tab$mois)[j])

#Aa1. Modèle LMM_1

# Formula is: Longueur  ~ ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP)
lmm1 <- lmer(long ~ ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP), data = tab1)

# Formula is: Longueur  ~ ZONE + (1 | POINT) + (1 | CAMP)
lmm1 <- lmer(long ~ ZONE + (1 | POINT) + (1 | CAMP), data = tab1)


#Estimation de la variance des intercepts ZONE:PERIODE par CAMP et par POINT, GLMM_1
print(summary(lmm1)) ## Ici, on r?cup?re la variance et la Std Dev. pour les effets random ET les Estimate et Std Error associ?s aux effets fixes
# test des effets random
lmm1_2 <- update(lmm1, . ~ . - (1 | POINT))
lmm1_3 <- update(lmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(lmm1, lmm1_2)
test_CAMP  <- anova(lmm1, lmm1_3)
print(test_POINT) ## Ici, test de l'effet Random POINT, on r?cup?re Df (=6/7), Chisq et Pr(>Chisq)
print(test_CAMP)
# test des effets fixes
print(Anova(lmm1))  ## Ici, on r?cup?re Chisq, Df et Pr(>Chsiq)
rm(lmm1_2, lmm1_3)

#Plot estimates (ranef) and standard-error
newdat                    <- expand.grid(ZONE = unique(tab1$ZONE))
newdat$long               <- predict(lmm1, newdat, re.form= NA, type = "response")
mm                        <- model.matrix(terms(lmm1), newdat)
pvar1                     <- diag(mm %*% tcrossprod(vcov(lmm1),mm))                   
tvar1                     <- pvar1 + as.vector(VarCorr(lmm1)$CAMP) + as.vector(VarCorr(lmm1)$POINT)
cmult                     <- 1.96

newdat  <- data.frame(
  newdat
  , plo = newdat$long-cmult*sqrt(pvar1)
  , phi = newdat$long+cmult*sqrt(pvar1)
  , tlo = newdat$long-cmult*sqrt(tvar1)
  , thi = newdat$long+cmult*sqrt(tvar1)
  )

newdat$PERIODE<-factor(c(NA,NA),c("AVANT","APRES"))
newdat$sp <- "edible crab"
newdat$month <- "june"
newdat <- newdat[, c("ZONE", "PERIODE", "long", "plo", "phi", "tlo", "thi", "sp", "month")]
stock <- rbind(stock, newdat)

######################## Mois= Sept ##################################
j=2
tab1 <- droplevels(tab[which(tab$mois == levels(tab$mois)[j]),])
print(levels(tab$mois)[j])

#Aa1. Modèle LMM_1

# Formula is: Longueur  ~ ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP)
lmm1 <- lmer(long ~ ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP), data = tab1)

#Estimation de la variance des intercepts ZONE:PERIODE par CAMP et par POINT, GLMM_1
print(summary(lmm1)) ## Ici, on r?cup?re la variance et la Std Dev. pour les effets random ET les Estimate et Std Error associ?s aux effets fixes
# test des effets random
lmm1_2 <- update(lmm1, . ~ . - (1 | POINT))
lmm1_3 <- update(lmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(lmm1, lmm1_2)
test_CAMP  <- anova(lmm1, lmm1_3)
print(test_POINT) ## Ici, test de l'effet Random POINT, on r?cup?re Df (=6/7), Chisq et Pr(>Chisq)
print(test_CAMP)
# test des effets fixes
print(Anova(lmm1))  ## Ici, on r?cup?re Chisq, Df et Pr(>Chsiq)
rm(lmm1_2, lmm1_3)



#Plot estimates (ranef) and standard-error
newdat                    <- expand.grid(ZONE = unique(tab1$ZONE), PERIODE = unique(tab1$PERIODE))
newdat$long               <- predict(lmm1, newdat, re.form= NA)
mm                        <- model.matrix(terms(lmm1), newdat)
pvar1                     <- diag(mm %*% tcrossprod(vcov(lmm1),mm))                   
tvar1                     <- pvar1 + as.vector(VarCorr(lmm1)$CAMP) + as.vector(VarCorr(lmm1)$POINT)
cmult                     <- 1.96

newdat  <- data.frame(
  newdat
  , plo = newdat$long-cmult*sqrt(pvar1)
  , phi = newdat$long+cmult*sqrt(pvar1)
  , tlo = newdat$long-cmult*sqrt(tvar1)
  , thi = newdat$long+cmult*sqrt(tvar1)
  )
newdat$PERIODE<-factor(newdat$PERIODE,c("AVANT","APRES"))
newdat$sp <- "edible crab"
newdat$month <- "sept."
stock <- rbind(stock, newdat)

######################################################################
#####////////////////    ARAIGNEE      \\\\\\\\\\\\\\\\\\\\\\\\\\#####
######################################################################
i=1
tab <- droplevels(TAB[which(TAB$espece == levels(TAB$espece)[i]),])
print(levels(TAB$espece)[i])

######################## Mois= Juin ##################################
j=1
tab1 <- droplevels(tab[which(tab$mois == levels(tab$mois)[j]),])
print(levels(tab$mois)[j])

#Aa1. Modèle LMM_1

# Formula is: Longueur  ~ ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP)
lmm1 <- lmer(long ~ ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP), data = tab1)

lmm1 <- lmer(long ~ 1 + (1 | POINT) + (1 | CAMP), data = tab1)


#Estimation de la variance des intercepts ZONE:PERIODE par CAMP et par POINT, GLMM_1
print(summary(lmm1)) ## Ici, on récupère la variance et la Std Dev. pour les effets random ET les Estimate et Std Error associés aux effets fixes
# test des effets random
lmm1_2 <- update(lmm1, . ~ . - (1 | POINT))
lmm1_3 <- update(lmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(lmm1, lmm1_2)
test_CAMP  <- anova(lmm1, lmm1_3)
print(test_POINT) ## Ici, test de l'effet Random POINT, on récupère Df (num/denom), Chisq et Pr(>Chisq)
print(test_CAMP)

######################## Mois= Sept ##################################
j=2
tab1 <- droplevels(tab[which(tab$mois == levels(tab$mois)[j]),])
print(levels(tab$mois)[j])

#Aa1. Modèle LMM_1

# Formula is: Longueur  ~ ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP)
lmm1 <- lmer(long ~ ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP), data = tab1)
#Estimation de la variance des intercepts ZONE:PERIODE par CAMP et par POINT, GLMM_1
print(summary(lmm1)) ## Ici, on récupère la variance et la Std Dev. pour les effets random ET les Estimate et Std Error associés aux effets fixes
# test des effets random
lmm1_2 <- update(lmm1, . ~ . - (1 | POINT))
lmm1_3 <- update(lmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(lmm1, lmm1_2)
test_CAMP  <- anova(lmm1, lmm1_3)
print(test_POINT) ## Ici, test de l'effet Random POINT, on récupère Df (num/denom), Chisq et Pr(>Chisq)
print(test_CAMP)
# test des effets fixes
print(Anova(lmm1))  ## Ici, on récupère Chisq, Df et Pr(>Chsiq)
rm(lmm1_2, lmm1_3)

#Plot estimates (ranef) and standard-error
newdat                    <- expand.grid(ZONE = unique(tab1$ZONE), PERIODE = unique(tab1$PERIODE))
newdat$long               <- predict(lmm1, newdat, re.form= NA)
mm                        <- model.matrix(terms(lmm1), newdat)
pvar1                     <- diag(mm %*% tcrossprod(vcov(lmm1),mm))                   
tvar1                     <- pvar1 + as.vector(VarCorr(lmm1)$CAMP) + as.vector(VarCorr(lmm1)$POINT)
cmult                     <- 1.96

newdat  <- data.frame(
  newdat
  , plo = newdat$long-cmult*sqrt(pvar1)
  , phi = newdat$long+cmult*sqrt(pvar1)
  , tlo = newdat$long-cmult*sqrt(tvar1)
  , thi = newdat$long+cmult*sqrt(tvar1)
)
newdat$PERIODE<-factor(newdat$PERIODE,c("AVANT","APRES"))
newdat$sp <- "spider crab"
newdat$month <- "sept."
stock <- rbind(stock, newdat)

write.table(stock, "./article/tabs/estimates.csv", row.names = F)
write.table(stock, "C:/Users/jnormand/Documents/boulot/56=ETUDES COMPLEMENTAIRES/article/tabs/estimates.csv", row.names = F)


# graphique 2 pour l'article
library(tidyr)
library(dplyr)
head(stock)
stock <- stock %>% mutate(ZONE2 = ifelse(ZONE == "DEHORS", "outside", "inside"), 
                          PERIODE2 = ifelse(PERIODE == "AVANT", "1before", "2after"))
stock$sp <- as.factor(stock$sp)
stock$month <- as.factor(stock$month)
levels(stock$sp) <- c("Edible crab", "Lobster", "Spider crab")
stock$sp <- factor(stock$sp, levels= c("Lobster", "Edible crab", "Spider crab"))
levels(stock$month) <- c("June", "September")
plot1<- ggplot(stock, aes(x= ZONE2, y= long, colour= PERIODE2))+
  geom_pointrange(aes(ymin= tlo, ymax= thi), position = position_dodge(0.1))+
  theme(panel.background = element_rect(fill= "white", colour = "black", size=1), 
        panel.grid.minor = element_blank(), 
        axis.title= element_text(size=14, colour= "black"), axis.text=element_text(size=12, colour= "black"),
        legend.title = element_text(size=14), legend.text = element_text(size=12),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 14))+
  scale_colour_discrete(name = "Period", labels = c("Before", "After", "No effect for period"))+
  labs(x= "Protected area", y= "Length (mm)")+facet_grid(sp~month, scales= "free")

tiff("//portenbessin/serha/02_OBSERVATION_SURVEILLANCE/05_IGA/20_CNPE_FLAMANVILLE/56=ETUDES COMPLEMENTAIRES/2019_CANTOFLAM/article/figs/Long_mm_221018.tif",
     units="in", width=8, height=8, res=300)
plot1
dev.off()
