# author: "Julien Normand"
# date: "3 aout 2021"

# packages 
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


## Loading dataset : all species and months
TAB <- read.table("~/git/CANTOFLAM/DATA/allspecies_allmonths_CPUE.csv", 
                  header=TRUE, sep=";", na.strings="NA", dec=".", strip.white=TRUE)

TAB[which(TAB$mois== "septembre"), "mois"] <- "sept"
TAB                                        <- droplevels(TAB)

tab2  <- read.table("~/git/CANTOFLAM/DATA/dist.csv",  header=TRUE, sep=";", na.strings="NA", dec=".", strip.white=TRUE)

# Formatting
TAB$ZONE     <- as.factor(substr(TAB$CANTO, 7,12))
TAB$PERIODE  <- as.factor(substr(TAB$CANTO, 1,5))
TAB$CAMP     <- as.factor(TAB$AN)
TAB$POINT    <- as.factor(TAB$POINT)
TAB          <- TAB[which(TAB$NB_CAS!= 0),]
TAB <- merge(x= TAB, y= tab2, by= "POINT")
TAB$DEHORS <- 0
TAB[which(TAB$ZONE == "DEHORS"), "DEHORS"] <- 1
TAB$DEDANS <- 0
TAB[which(TAB$ZONE == "DEDANS"), "DEDANS"] <- 1
TAB$DEHORS <- as.factor(TAB$DEHORS)
TAB$DEDANS <- as.factor(TAB$DEDANS)
#TAB[which(TAB$ZONE == "DEDANS"), "dist"] <- 0
moy  <- mean(TAB$dist, na.rm= T)
et   <- sd(TAB$dist, na.rm= T)
TAB$scale_dist  <- (TAB$dist- moy)/et
TAB$scale_annee <- (TAB$AN - mean(TAB$AN)) / sd(TAB$AN)
TAB$CPUE <- TAB$NB_CAPT / TAB$NB_CAS
TAB[which(TAB$CPUE >= 10000), "CPUE"] <- 0
rm(tab2, et, moy)
TAB$espece <- as.factor(TAB$espece)
TAB$mois <- as.factor(TAB$mois)
summary(TAB)

#### Model fitting by species ###
############################################################################################################
#sp = lobster 
i=3
tab <- droplevels(TAB[which(TAB$espece == levels(TAB$espece)[i]),])
print(levels(TAB$espece)[i])
#month= june
j=1
tab1 <- droplevels(tab[which(tab$mois == levels(tab$mois)[j]),])
print(levels(tab$mois)[j])
glmm1 <- glmer.nb(NB_CAPT ~ offset(log(NB_CAS)) + ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP), data = tab1)
#Interceps variance estimation  ZONE:PERIODE by CAMP(survey) and POINT (sampling location)
print(summary(glmm1))
# test des effets random
glmm1_2 <- update(glmm1, . ~ . - (1 | POINT))
glmm1_3 <- update(glmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(glmm1, glmm1_2)
test_CAMP  <- anova(glmm1, glmm1_3)
print(test_POINT)
print(test_CAMP)
# Fixed effects test
print(Anova(glmm1))  ## Extract Chisq, Df et Pr(>Chsiq)
rm(glmm1_2, glmm1_3)
#Plot estimates (ranef) and standard-error
newdat          <- expand.grid(ZONE = unique(tab1$ZONE), PERIODE = unique(tab1$PERIODE), NB_CAPT = 1, NB_CAS = 1)
newdat$NB_CAPT  <- predict(glmm1, newdat, type = "response", re.form= NA)
mm              <- model.matrix(terms(glmm1), newdat)
pvar1           <- diag(mm %*% tcrossprod(vcov(glmm1),mm))                   
tvar1           <- pvar1 + as.vector(VarCorr(glmm1)$CAMP) + as.vector(VarCorr(glmm1)$POINT)
cmult   <- 1.96
newdat  <- data.frame(
  newdat
  , plo = newdat$NB_CAPT-cmult*sqrt(pvar1)
  , phi = newdat$NB_CAPT+cmult*sqrt(pvar1)
  , tlo = newdat$NB_CAPT-cmult*sqrt(tvar1)
  , thi = newdat$NB_CAPT+cmult*sqrt(tvar1)
  )
newdat$PERIODE<-factor(newdat$PERIODE,c("AVANT","APRES"))
newdat$sp <- "lobster"
newdat$month <- "june"
stock <- newdat
g0 <- ggplot(newdat, aes(x= ZONE, y= NB_CAPT, colour= PERIODE)) 
g0 <- g0 + geom_pointrange(aes(ymin= tlo, ymax= thi), position = position_dodge(0.1))
g0 <- g0 + theme(panel.background = element_rect(fill= "white", colour = "black", linewidth =1), 
                 panel.grid.minor = element_blank(), 
                 axis.title= element_text(size=14, colour= "black"), axis.text=element_text(size=12, colour= "black"),
                 legend.title = element_text(size=14), legend.text = element_text(size=12)) 
g0 <- g0 + labs(x= "ZONE", y= "CPUE")
# png("rapport_acti_oct2019/figs/fig2_I1A.png", width= 300, height = 300)
# g0
# dev.off()
######################################################################################################
#Month= september
j=2
tab1 <- droplevels(tab[which(tab$mois == levels(tab$mois)[j]),])
print(levels(tab$mois)[j])
glmm1 <- glmer.nb(NB_CAPT ~ offset(log(NB_CAS)) + ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP), data = tab1)
#Interceps variance estimation  ZONE:PERIODE by CAMP(survey) and POINT (sampling location)
print(summary(glmm1))
# Random effect test
glmm1_2 <- update(glmm1, . ~ . - (1 | POINT))
glmm1_3 <- update(glmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(glmm1, glmm1_2)
test_CAMP  <- anova(glmm1, glmm1_3)
print(test_POINT) ## Ici, test de l'effet Random POINT, on r?cup?re Df (=6/7), Chisq et Pr(>Chisq)
print(test_CAMP)
#  Fixed effects test
print(Anova(glmm1))  ## Ici, on r?cup?re Chisq, Df et Pr(>Chsiq)
rm(glmm1_2, glmm1_3)
#Plot estimates (ranef) and standard-error
newdat          <- expand.grid(ZONE = unique(tab1$ZONE), PERIODE = unique(tab1$PERIODE), NB_CAPT = 1, NB_CAS = 1)
newdat$NB_CAPT  <- predict(glmm1, newdat, re.form= NA, type = "response")
mm              <- model.matrix(terms(glmm1), newdat)
pvar1           <- diag(mm %*% tcrossprod(vcov(glmm1),mm))                   
tvar1           <- pvar1 + as.vector(VarCorr(glmm1)$CAMP) + as.vector(VarCorr(glmm1)$POINT)
cmult   <- 1.96
newdat  <- data.frame(
  newdat
  , plo = newdat$NB_CAPT-cmult*sqrt(pvar1)
  , phi = newdat$NB_CAPT+cmult*sqrt(pvar1)
  , tlo = newdat$NB_CAPT-cmult*sqrt(tvar1)
  , thi = newdat$NB_CAPT+cmult*sqrt(tvar1)
  )
newdat$PERIODE<-factor(newdat$PERIODE,c("AVANT","APRES"))
newdat$sp <- "lobster"
newdat$month <- "sept."
stock <- rbind(stock, newdat)
g1 <- ggplot(newdat, aes(x= ZONE, y= NB_CAPT, colour= PERIODE)) 
g1 <- g1 + geom_pointrange(aes(ymin= tlo, ymax= thi), position = position_dodge(0.1))
g1 <- g1 + theme(panel.background = element_rect(fill= "white", colour = "black", size=1), 
                 panel.grid.minor = element_blank(), 
                 axis.title= element_text(size=14, colour= "black"), axis.text=element_text(size=12, colour= "black"),
                 legend.title = element_text(size=14), legend.text = element_text(size=12)) 
g1 <- g1 + labs(x= "ZONE", y= "CPUE")
# png("rapport_acti_oct2019/figs/fig2_I1B.png", width= 300, height = 300)
# g1
# dev.off()
################################################################################################
# sp = edible crab 
i=4
tab <- droplevels(TAB[which(TAB$espece == levels(TAB$espece)[i]),])
print(levels(TAB$espece)[i])
#Month = June
j=1
tab1 <- droplevels(tab[which(tab$mois == levels(tab$mois)[j]),])
print(levels(tab$mois)[j])
glmm1 <- glmer.nb(NB_CAPT ~ offset(log(NB_CAS)) + ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP), data = tab1)
#Interceps variance estimation  ZONE:PERIODE by CAMP(survey) and POINT (sampling location)
print(summary(glmm1))
# test random effect
glmm1_2 <- update(glmm1, . ~ . - (1 | POINT))
glmm1_3 <- update(glmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(glmm1, glmm1_2)
test_CAMP  <- anova(glmm1, glmm1_3)
print(test_POINT) ##  Random POINT effect, extract Df (=6/7), Chisq et Pr(>Chisq)
print(test_CAMP)
# Fixed effect test 
print(Anova(glmm1))  ## Extract Chisq, Df et Pr(>Chsiq)
rm(glmm1_2, glmm1_3)

#Plot estimates (ranef) and standard-error
newdat          <- expand.grid(ZONE = unique(tab1$ZONE), PERIODE = unique(tab1$PERIODE), NB_CAPT = 1, NB_CAS = 1)
newdat$NB_CAPT  <- predict(glmm1, newdat, re.form= NA, type = "response")
mm              <- model.matrix(terms(glmm1), newdat)
pvar1           <- diag(mm %*% tcrossprod(vcov(glmm1),mm))                   
tvar1           <- pvar1 + as.vector(VarCorr(glmm1)$CAMP) + as.vector(VarCorr(glmm1)$POINT)
cmult   <- 1.96
newdat  <- data.frame(
  newdat
  , plo = newdat$NB_CAPT-cmult*sqrt(pvar1)
  , phi = newdat$NB_CAPT+cmult*sqrt(pvar1)
  , tlo = newdat$NB_CAPT-cmult*sqrt(tvar1)
  , thi = newdat$NB_CAPT+cmult*sqrt(tvar1)
  )
newdat$PERIODE<-factor(newdat$PERIODE,c("AVANT","APRES"))
newdat$sp <- "edible crab"
newdat$month <- "june"
stock <- rbind(stock, newdat)
g2 <- ggplot(newdat, aes(x= ZONE, y= NB_CAPT, colour= PERIODE)) 
g2 <- g2 + geom_pointrange(aes(ymin= tlo, ymax= thi), position = position_dodge(0.1))
g2 <- g2 + theme(panel.background = element_rect(fill= "white", colour = "black", size=1), 
                 panel.grid.minor = element_blank(), 
                 axis.title= element_text(size=14, colour= "black"), axis.text=element_text(size=12, colour= "black"),
                 legend.title = element_text(size=14), legend.text = element_text(size=12)) 
g2 <- g2 + labs(x= "ZONE", y= "CPUE")
# png("rapport_acti_oct2019/figs/fig2_I2A.png", width= 300, height = 300)
# g2
# dev.off()
######################################################################
#Month= september
j=2
tab1 <- droplevels(tab[which(tab$mois == levels(tab$mois)[j]),])
print(levels(tab$mois)[j])
glmm1 <- glmer.nb(NB_CAPT ~ offset(log(NB_CAS)) + ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP), data = tab1)
#Interceps variance estimation  ZONE:PERIODE by CAMP(survey) and POINT (sampling location)
print(summary(glmm1)) 
# test des effets random
glmm1_2 <- update(glmm1, . ~ . - (1 | POINT))
glmm1_3 <- update(glmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(glmm1, glmm1_2)
test_CAMP  <- anova(glmm1, glmm1_3)
print(test_POINT)
print(test_CAMP)
# test fixed effect
print(Anova(glmm1))
rm(glmm1_2, glmm1_3)

newdat          <- expand.grid(ZONE = unique(tab1$ZONE), PERIODE = unique(tab1$PERIODE), NB_CAPT = 1, NB_CAS = 1)
newdat$NB_CAPT  <- predict(glmm1, newdat, re.form= NA, type = "response")
mm              <- model.matrix(terms(glmm1), newdat)
pvar1           <- diag(mm %*% tcrossprod(vcov(glmm1),mm))                   
tvar1           <- pvar1 + as.vector(VarCorr(glmm1)$CAMP) + as.vector(VarCorr(glmm1)$POINT)
cmult   <- 1.96
newdat  <- data.frame(
  newdat
  , plo = newdat$NB_CAPT-cmult*sqrt(pvar1)
  , phi = newdat$NB_CAPT+cmult*sqrt(pvar1)
  , tlo = newdat$NB_CAPT-cmult*sqrt(tvar1)
  , thi = newdat$NB_CAPT+cmult*sqrt(tvar1)
  )
newdat$PERIODE<-factor(newdat$PERIODE,c("AVANT","APRES"))
newdat$sp <- "edible crab"
newdat$month <- "sept."
stock <- rbind(stock, newdat)

g3 <- ggplot(newdat, aes(x= ZONE, y= NB_CAPT, colour= PERIODE)) 
g3 <- g3 + geom_pointrange(aes(ymin= tlo, ymax= thi), position = position_dodge(0.1))
g3 <- g3 + theme(panel.background = element_rect(fill= "white", colour = "black", size=1), 
                 panel.grid.minor = element_blank(), 
                 axis.title= element_text(size=14, colour= "black"), axis.text=element_text(size=12, colour= "black"),
                 legend.title = element_text(size=14), legend.text = element_text(size=12)) 
g3 <- g3 + labs(x= "ZONE", y= "CPUE")

# png("rapport_acti_oct2019/figs/fig2_I2B.png", width= 300, height = 300)
# g3
# dev.off()
################################################################################################
#    sp = spider crab 
i=1
tab <- droplevels(TAB[which(TAB$espece == levels(TAB$espece)[i]),])
print(levels(TAB$espece)[i])
#Month= June
j=1
tab1 <- droplevels(tab[which(tab$mois == levels(tab$mois)[j]),])
print(levels(tab$mois)[j])
glmm1 <- glmer.nb(NB_CAPT ~ offset(log(NB_CAS)) + ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP), data = tab1)
print(summary(glmm1)) 
# test random effects
glmm1_2 <- update(glmm1, . ~ . - (1 | POINT))
glmm1_3 <- update(glmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(glmm1, glmm1_2)
test_CAMP  <- anova(glmm1, glmm1_3)
print(test_POINT) ## 
print(test_CAMP)
# test fixed effects
print(Anova(glmm1)) 
rm(glmm1_2, glmm1_3)
## Interaction ZONE:PERIODE NS and so Model reduction
glmm1 <- glmer.nb(NB_CAPT ~ offset(log(NB_CAS)) + ZONE + PERIODE + (1 | POINT) + (1 | CAMP), data = tab1)

#Interceps variance estimation  ZONE:PERIODE by CAMP(survey) and POINT (sampling location)
print(summary(glmm1))
# rest random effects
glmm1_2 <- update(glmm1, . ~ . - (1 | POINT))
glmm1_3 <- update(glmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(glmm1, glmm1_2)
test_CAMP  <- anova(glmm1, glmm1_3)
print(test_POINT) 
print(test_CAMP)
# test fixed effects
print(Anova(glmm1))  
rm(glmm1_2, glmm1_3)

#Plot estimates (ranef) and standard-error
newdat          <- expand.grid(ZONE = unique(tab1$ZONE), PERIODE = unique(tab1$PERIODE), NB_CAPT = 1, NB_CAS = 1)
newdat$NB_CAPT  <- predict(glmm1, newdat, re.form= NA, type = "response")
mm              <- model.matrix(terms(glmm1), newdat)
pvar1           <- diag(mm %*% tcrossprod(vcov(glmm1),mm))                   
tvar1           <- pvar1 + as.vector(VarCorr(glmm1)$CAMP) + as.vector(VarCorr(glmm1)$POINT)
cmult   <- 1.96
newdat  <- data.frame(
  newdat
  , plo = newdat$NB_CAPT-cmult*sqrt(pvar1)
  , phi = newdat$NB_CAPT+cmult*sqrt(pvar1)
  , tlo = newdat$NB_CAPT-cmult*sqrt(tvar1)
  , thi = newdat$NB_CAPT+cmult*sqrt(tvar1)
)
newdat$PERIODE<-factor(newdat$PERIODE,c("AVANT","APRES"))
newdat$sp <- "spider crab"
newdat$month <- "june"
stock <- rbind(stock, newdat)

g4 <- ggplot(newdat, aes(x= ZONE, y= NB_CAPT, colour= PERIODE)) 
g4 <- g4 + geom_pointrange(aes(ymin= tlo, ymax= thi), position = position_dodge(0.1))
g4 <- g4 + theme(panel.background = element_rect(fill= "white", colour = "black", size=1), 
                 panel.grid.minor = element_blank(), 
                 axis.title= element_text(size=14, colour= "black"), axis.text=element_text(size=12, colour= "black"),
                 legend.title = element_text(size=14), legend.text = element_text(size=12)) 
g4 <- g4 + labs(x= "ZONE", y= "CPUE")
 # png("rapport_acti_oct2019/figs/fig2_I3A.png", width= 300, height = 300)
# g4
# dev.off()

######################################################################
#Month= september
j=2
tab1 <- droplevels(tab[which(tab$mois == levels(tab$mois)[j]),])
print(levels(tab$mois)[j])

glmm1 <- glmer.nb(NB_CAPT ~ offset(log(NB_CAS)) + ZONE + PERIODE + ZONE:PERIODE + (1 | POINT) + (1 | CAMP), data = tab1)

#Interceps variance estimation  ZONE:PERIODE by CAMP(survey) and POINT (sampling location)
print(summary(glmm1)) 
# test des effets random
glmm1_2 <- update(glmm1, . ~ . - (1 | POINT))
glmm1_3 <- update(glmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(glmm1, glmm1_2)
test_CAMP  <- anova(glmm1, glmm1_3)
print(test_POINT) 
print(test_CAMP)
# test des effets fixes
print(Anova(glmm1))  ## Ici, on r?cup?re Chisq, Df et Pr(>Chsiq)
rm(glmm1_2, glmm1_3)

## PERIODE & ZONE:PERIODE = NS, so model reduction
glmm1 <- glmer.nb(NB_CAPT ~ offset(log(NB_CAS)) + ZONE + (1 | POINT) + (1 | CAMP), data = tab1)

#Estimation de la variance des intercepts ZONE:PERIODE par CAMP et par POINT, GLMM_1
print(summary(glmm1)) 
# test des effets random
glmm1_2 <- update(glmm1, . ~ . - (1 | POINT))
glmm1_3 <- update(glmm1, . ~ . - (1 | CAMP))
test_POINT <- anova(glmm1, glmm1_2)
test_CAMP  <- anova(glmm1, glmm1_3)
print(test_POINT) 
print(test_CAMP)
# test des effets fixes
print(Anova(glmm1))  
rm(glmm1_2, glmm1_3)

#Plot estimates (ranef) and standard-error
newdat          <- expand.grid(ZONE = unique(tab1$ZONE), NB_CAPT = 1, NB_CAS = 1)
newdat$NB_CAPT  <- predict(glmm1, newdat, re.form= NA, type = "response")
mm              <- model.matrix(terms(glmm1), newdat)
pvar1           <- diag(mm %*% tcrossprod(vcov(glmm1),mm))                   
tvar1           <- pvar1 + as.vector(VarCorr(glmm1)$CAMP) + as.vector(VarCorr(glmm1)$POINT)
cmult   <- 1.96
newdat  <- data.frame(
  newdat
  , plo = newdat$NB_CAPT-cmult*sqrt(pvar1)
  , phi = newdat$NB_CAPT+cmult*sqrt(pvar1)
  , tlo = newdat$NB_CAPT-cmult*sqrt(tvar1)
  , thi = newdat$NB_CAPT+cmult*sqrt(tvar1)
)
newdat$PERIODE<-factor(c(NA,NA),c("AVANT","APRES"))
newdat$sp <- "spider crab"
newdat$month <- "sept."
newdat <- newdat[, c("ZONE", "PERIODE", "NB_CAPT", "NB_CAS", "plo", "phi", "tlo", "thi", "sp", "month")]
stock <- rbind(stock, newdat)
#write.table(stock, "./article/tabs/estimates.csv", row.names = F)
#write.table(stock, "C:/Users/jnormand/Documents/boulot/56=ETUDES COMPLEMENTAIRES/article/tabs/estimates.csv", row.names = F)

#stock %>% filter(sp== "edible crab", PERIODE== "APRES") %>% View()
g5 <- ggplot(newdat, aes(x= ZONE, y= NB_CAPT, colour= "")) 
g5 <- g5 + geom_pointrange(aes(ymin= tlo, ymax= thi), position = position_dodge(0.1) )
g5 <- g5 + theme(panel.background = element_rect(fill= "white", colour = "black", size=1), 
                 panel.grid.minor = element_blank(), 
                 axis.title= element_text(size=14, colour= "black"), axis.text=element_text(size=12, colour= "black"),
                 legend.position = "none") 
g5 <- g5 + labs(x= "ZONE", y= "CPUE")
# png("rapport_acti_oct2019/figs/fig2_I3B.png", width= 300, height = 300)
# g5
# dev.off()

#########################################
####         PAPER                    ###
#########################################

# PAPER : Figure 1
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

plot1 <- ggplot(stock, aes(x= ZONE2, y= NB_CAPT, colour= PERIODE2))+
  geom_pointrange(aes(ymin= tlo, ymax= thi), position = position_dodge(0.1))+
  theme(panel.background = element_rect(fill= "white", colour = "black", size=1), 
        panel.grid.minor = element_blank(), 
        axis.title= element_text(size=14, colour= "black"), axis.text=element_text(size=12, colour= "black"),
        legend.title = element_text(size=14), legend.text = element_text(size=12),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 14))+
  #geom_text(aes(y = pos_letter, x= 2.3, label = letter), size= 12, colour= "black")+
  scale_colour_discrete(name = "Period", labels = c("Before", "After", "No effect for period"))+
  labs(x= "Protected area", y= "CPUE")+facet_grid(sp~month, scales= "free")


plot1
dev.off()

#######################################
# Computed values present in the paper#
#######################################
#U: used in the paper
#NU: not used in the paper

# lobster

#NU
### PAPER : computation of the lobser CPUE between 2000 and 2005, then between 2000 and 2018 inside the MPA
pipo <- TAB %>% filter(espece== "homard" & AN %in% c(2000,2005,2017) & ZONE== "DEDANS") %>%
  group_by(mois, CAMP) %>% summarise(mean_CPUE= mean(CPUE, na.rm= T)) %>%
  pivot_wider(names_from= "CAMP", values_from = "mean_CPUE")
colnames(pipo)[2:4] <- c("creation", "ans5", "ans17")  
pipo %>%  mutate(effet_5ans = (ans5-creation)*100/creation, effet_18ans= (ans17-creation)*100/creation)

TAB%>% group_by(POINT) %>% summarise(mean_dist= mean(dist)) %>% arrange(mean_dist)

#U
## Paper, Change in lobster CPUE comparison between 2000, 2005 and 2017 inside the MPA 
pipo <- TAB %>% filter(espece== "homard" & AN %in% c(2000,2005, 2017) & ZONE== "DEDANS") %>%
  group_by(mois, CAMP) %>% summarise(mean_CPUE= mean(CPUE, na.rm= T)) %>%
  tidyr::pivot_wider(names_from= c( "CAMP"), values_from = "mean_CPUE")
colnames(pipo)[2:4] <- c("creation", "interm", "final")  

pipo %>%  mutate(effet_interm = (interm-creation)*100/creation, 
                 effet_final= (final-creation)*100/creation)



### Paper: mean lobster CPUE increase between 2000 and 2018 for a pont close to the reserve (9 = 613m)
# and the furthest point (6= 6335m)
TAB %>% filter(ZONE== "DEHORS") %>% group_by(POINT) %>% summarise(mean_dist= mean(dist)) %>% arrange(mean_dist)

pipo <- TAB %>% filter(espece== "homard" & AN %in% c(2000,2017) & POINT %in% c(6,9)) %>%
  group_by(mois, CAMP, POINT) %>% summarise(mean_CPUE= mean(CPUE, na.rm= T)) %>%
  pivot_wider(names_from= c("POINT", "CAMP"), values_from = "mean_CPUE")
colnames(pipo)[2:5] <- c("creation_6", "creation_9", "ans17_6", "ans17_9")  

pipo %>%  mutate(effet_site6 = (ans17_6-creation_6)*100/creation_6, effet_site9= (ans17_9-creation_9)*100/creation_9)


#edible crab 
###PAPER : computation of the e. crab CPUE between 2000 and 2018, for the point closest 
#to the MPA center (8 = 143m) and the point nearest from the MPA boundaries (10 = 1533)
pipo <- TAB %>% filter(espece== "tourteau" & POINT %in% c(8,10)) %>% 
  group_by(mois, POINT) %>% summarise(mean_CPUE= mean(CPUE, na.rm= T)) %>%
  tidyr::pivot_wider(names_from= c("POINT"), values_from = "mean_CPUE")
colnames(pipo)[2:3] <- c("point_8", "point_10")  
pipo %>%  mutate(effet_site6 = (point_10-point_8)*100/point_8)

###PAPER : computation of the e. crab CPUE between 2000 and 2018,
# U
#for the point nearest from the MPA boundaries (10 = 1533) and the point furthest from the MPA boundaries  (6 = 6335)
pipo <- TAB %>% filter(espece== "tourteau" & POINT %in% c(10,6)) %>% 
  group_by(mois, POINT) %>% summarise(mean_CPUE= mean(CPUE, na.rm= T)) %>%
  tidyr::pivot_wider(names_from= c("POINT"), values_from = "mean_CPUE")
colnames(pipo)[2:3] <- c("point_6", "point_10")  
pipo %>%  mutate(effet_site6 = (point_6-point_10)*100/point_10)

#NU
### PAPER: mean edible crab CPUE variation over 2000-2018, period between the MPA border 
# sample point (10 = 1533) and the furthest point (6 = 6335)
pipo <- TAB %>% filter(espece== "tourteau" & AN %in% c(2000,2017)) %>% 
  group_by(mois, AN) %>% summarise(mean_CPUE= mean(CPUE, na.rm= T)) %>%
  tidyr::pivot_wider(names_from= c("AN"), values_from = "mean_CPUE")
colnames(pipo)[2:3] <- c("point_6", "point_10")  
pipo %>%  mutate(effet_site6 = (point_6-point_10)*100/point_10)

#U
#to the MPA center (8 = 143m) and the point nearest from the MPA boundaries (10 = 1533)
pipo <- TAB %>% filter(espece== "tourteau" & POINT %in% c(8,10)) %>% 
  group_by(mois, POINT) %>% summarise(mean_CPUE= mean(CPUE, na.rm= T)) %>%
  tidyr::pivot_wider(names_from= c("POINT"), values_from = "mean_CPUE")
colnames(pipo)[2:3] <- c("point_8", "point_10")  
pipo %>%  mutate(effet_site6 = (point_10-point_8)*100/point_8)

#spider crab
#U
### computation of the s. crab CPUE between 2000 and 2018, for the point closest 
#to the MPA center (8 = 143m) and the point nearest from the MPA boundaries (10 = 1533)
pipo <- TAB %>% filter(espece== "araignee" & mois== "juin" & POINT %in% c(8,10)) %>% 
  group_by(POINT) %>% summarise(mean_CPUE= mean(CPUE, na.rm= T)) %>%
  tidyr::pivot_wider(names_from= c("POINT"), values_from = "mean_CPUE")
colnames(pipo) <- c("point_8", "point_10")  
pipo %>%  mutate(effet_site6 = (point_10-point_8)*100/point_8)

#NU
### computation of the s. crab CPUE between 2000 and 2018, for the point closest 
#to the MPA center (8 = 143m) and and the furthest point (10 = 1533)
pipo <- TAB %>% filter(espece== "araignee" & mois== "sept" & POINT %in% c(8,10)) %>% 
  group_by(POINT) %>% summarise(mean_CPUE= mean(CPUE, na.rm= T)) %>%
  tidyr::pivot_wider(names_from= c("POINT"), values_from = "mean_CPUE")
colnames(pipo) <- c("point_8", "point_10")  
pipo %>%  mutate(effet_site6 = (point_10-point_8)*100/point_10)

#U
### mean spider crab  CPUE variation over 2000-2017 period in june 
pipo <- TAB %>% filter(espece== "araignee" & AN %in% c(2000,2017) & mois == "juin") %>% 
  group_by(mois, AN) %>% summarise(mean_CPUE= mean(CPUE, na.rm= T)) %>%
  tidyr::pivot_wider(names_from= c("AN"), values_from = "mean_CPUE")
colnames(pipo)[2:3] <- c("creation", "final")  
pipo %>%  mutate(effet_reserve = (final-creation)*100/creation)

#U
### mean spider crab  CPUE variation over 2000-2005 period and then 2005-2017 in september
pipo <- TAB %>% filter(espece== "araignee" & AN %in% c(2000,2005,2017) & mois == "sept") %>% 
  group_by(mois, AN) %>% summarise(mean_CPUE= mean(CPUE, na.rm= T)) %>%
  tidyr::pivot_wider(names_from= c("AN"), values_from = "mean_CPUE")
colnames(pipo)[2:4] <- c("creation", "intermediaire", "final")  
pipo %>%  mutate(effet_reserve_int = (intermediaire-creation)*100/creation)
pipo %>%  mutate(effet_reserve_fin = (final-intermediaire)*100/intermediaire)

