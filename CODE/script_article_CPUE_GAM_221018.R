# author: "Juju"
# date: "19 oct. 2022"
# Attention : fais appel à la fonction plotGAM2
rm(list= ls())
library(ggplot2) 
library(lme4)
library(lattice)
library(reshape2)
library(gridExtra)
# library(dfoptim)
# library(optimx)
# library(nloptr)
# library(knitr)
library(MASS)
library(car)
library(tidyr)
library(dplyr)
library(mgcv)
# library(vroom)
# library(forcats)
# library(visreg)
# library(tidymv)
library(ggeffects)
library(voxel)
library(gridExtra)
library(ISLR)
# library(purrr)
# Nombre de captures
setwd("C:/Users/jnormand/Documents/boulot/56=ETUDES COMPLEMENTAIRES")
#setwd("//portenbessin/serha/02_OBSERVATION_SURVEILLANCE/05_IGA/20_CNPE_FLAMANVILLE/56=ETUDES COMPLEMENTAIRES/2019_CANTOFLAM")
## je charge le jeu de donn?es complet des CPUE : toutes les esp?ces, tous les mois
TAB <- read.table("data/data_analyses_explo/CPUE/allspecies_allmonths_CPUE.csv", 
                  header=TRUE, sep=";", na.strings="NA", dec=".", strip.white=TRUE)
TAB[which(TAB$mois== "septembre"), "mois"] <- "sept"
TAB                                        <- droplevels(TAB)
tab2  <- read.table("data/data_analyses_explo/dist.csv",  header=TRUE, sep=";", na.strings="NA", dec=".", strip.white=TRUE)
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
rm(et, moy)
###################################################################################
tab <- TAB %>% filter(espece %in% c("araignee","homard","tourteau"), AN >= 2000) %>%
  dplyr::select(NB_CAPT, NB_CAS, CPUE, espece, POINT, dist, AN, mois)
stock <- TAB %>% select(POINT, dist) %>% group_by(POINT) %>% summarise(dist_2= mean(dist)) %>% ungroup()
#####################////////////////////////////////////Homard - adjust GAMM et estimate Marginal effects pour la table 4
espece         <- "homard"
tab1           <- droplevels(tab[which(tab$espece == espece),])
tab1$mois <- as.factor(tab1$mois)
homard_gam         <- gam(CPUE ~ mois + s(dist, by= mois) + s(AN, by= mois) + s(POINT, bs = 're', by= mois), 
                      data = tab1,
                      method = 'REML',family=nb(link="log"))
pipo <- as.data.frame(ggpredict(homard_gam, terms = c("dist[142.6506, 6334.8674]", "AN[2000, 2018]", "mois")))
pipo$dist_2 <- pipo$x
pipo$espece <- espece
pipo$point <- 8
pipo[which(pipo$dist_2 != 142.6506), "point"] <- 6
stock <- pipo
head(stock)
stock <- stock %>% dplyr::select(predicted, group, facet, espece, point) %>% 
  rename(pred_CPUE = "predicted",
         an = "group",
         saison = "facet")
# evolution des CPUE entre 2000 et 2018, pour un point "dedans" (le point 8) et un point "dehors" (le point 6)
pipo.an <- stock %>% pivot_wider(names_from = c("an", "point"), values_from = c("pred_CPUE")) %>%
  rename(avant.dedans = "2000_8", apres.dedans = "2018_8", avant.dehors = "2018_6", apres.dehors = "2018_6") %>%
  mutate(effet.an.dedans = (apres.dedans - avant.dedans)*100/avant.dedans,
         effet.an.dehors = (apres.dehors - avant.dehors)*100/avant.dehors) %>%
  dplyr::select(saison, espece, effet.an.dedans, effet.an.dehors)
# evolution des CPUE entre le point le plus proche du centroide (le point 8) et le point qui en est le plus 
# eloigne(le point 6) juste après la mise en place du cantonnement (en 2000) et à la toute fin de la série (en 2018)
pipo.dist <- stock %>% pivot_wider(names_from = c("an", "point"), values_from = c("pred_CPUE")) %>%
  rename(avant.dedans = "2000_8", apres.dedans = "2018_8", avant.dehors = "2000_6", apres.dehors = "2018_6") %>%
  mutate(effet.dist.avant = (avant.dehors - avant.dedans)*100/avant.dedans,
         effet.dist.apres = (apres.dehors - apres.dedans)*100/apres.dedans) %>%
  dplyr::select(saison, espece, effet.dist.avant, effet.dist.apres)
stock2 <- full_join(pipo.an, pipo.dist)
#####################////////////////////////////////////Tourteau - adjust GAMM et estimate Marginal effects pour la table 4
espece         <- "tourteau"
tab1           <- droplevels(tab[which(tab$espece == espece),])
tab1$mois <- as.factor(tab1$mois)
tourteau_gam         <- gam(CPUE ~ mois + s(dist, by= mois) + s(AN, by= mois) + s(POINT, bs = 're', by= mois), 
                          data = tab1,
                          method = 'REML',family=nb(link="log"))
pipo <- as.data.frame(ggpredict(tourteau_gam, terms = c("dist[142.6506, 6334.8674]", "AN[2000, 2018]", "mois")))
pipo$dist_2 <- pipo$x
pipo$espece <- espece
pipo$point <- 8
pipo[which(pipo$dist_2 != 142.6506), "point"] <- 6
stock <- pipo
stock <- stock %>% dplyr::select(predicted, group, facet, espece, point) %>% 
  rename(pred_CPUE = "predicted",
         an = "group",
         saison = "facet")
# evolution des CPUE entre 2000 et 2018, pour un point "dedans" (le point 8) et un point "dehors" (le point 6)
pipo.an <- stock %>% pivot_wider(names_from = c("an", "point"), values_from = c("pred_CPUE")) %>%
  rename(avant.dedans = "2000_8", apres.dedans = "2018_8", avant.dehors = "2000_6", apres.dehors = "2018_6") %>%
  mutate(effet.an.dedans = (apres.dedans - avant.dedans)*100/avant.dedans,
         effet.an.dehors = (apres.dehors - avant.dehors)*100/avant.dehors) %>%
  dplyr::select(saison, espece, effet.an.dedans, effet.an.dehors)
# evolution des CPUE entre le point le plus proche du centroide (le point 8) et le point qui en est le plus 
# eloigne(le point 6) juste après la mise en place du cantonnement (en 2000) et à la toute fin de la série (en 2018)
pipo.dist <- stock %>% pivot_wider(names_from = c("an", "point"), values_from = c("pred_CPUE")) %>%
  rename(avant.dedans = "2000_8", apres.dedans = "2018_8", avant.dehors = "2000_6", apres.dehors = "2018_6") %>%
  mutate(effet.dist.avant = (avant.dehors - avant.dedans)*100/avant.dedans,
         effet.dist.apres = (apres.dehors - apres.dedans)*100/apres.dedans) %>%
  dplyr::select(saison, espece, effet.dist.avant, effet.dist.apres)
stock3 <- full_join(pipo.an, pipo.dist)
stock2 <- rbind(stock2, stock3)
#####################////////////////////////////////////Araignée - adjust GAMM et estimate Marginal effects pour la table 4
espece         <- "araignee"
tab1           <- droplevels(tab[which(tab$espece == espece),])
tab1$mois <- as.factor(tab1$mois)
araignee_gam         <- gam(CPUE ~ mois + s(dist, by= mois) + s(AN, by= mois) + s(POINT, bs = 're', by= mois), 
                            data = tab1,
                            method = 'REML',family=nb(link="log"))
pipo <- as.data.frame(ggpredict(araignee_gam, terms = c("dist[142.6506, 6334.8674]", "AN[2000, 2018]", "mois")))
pipo$dist_2 <- pipo$x
pipo$espece <- espece
pipo$point <- 8
pipo[which(pipo$dist_2 != 142.6506), "point"] <- 6
stock <- pipo
stock <- stock %>% dplyr::select(predicted, group, facet, espece, point) %>% 
  rename(pred_CPUE = "predicted",
         an = "group",
         saison = "facet")
# evolution des CPUE entre 2000 et 2018, pour un point "dedans" (le point 8) et un point "dehors" (le point 6)
pipo.an <- stock %>% pivot_wider(names_from = c("an", "point"), values_from = c("pred_CPUE")) %>%
  rename(avant.dedans = "2000_8", apres.dedans = "2018_8", avant.dehors = "2000_6", apres.dehors = "2018_6") %>%
  mutate(effet.an.dedans = (apres.dedans - avant.dedans)*100/avant.dedans,
         effet.an.dehors = (apres.dehors - avant.dehors)*100/avant.dehors) %>%
  dplyr::select(saison, espece, effet.an.dedans, effet.an.dehors)
# evolution des CPUE entre le point le plus proche du centroide (le point 8) et le point qui en est le plus 
# eloigne(le point 6) juste après la mise en place du cantonnement (en 2000) et à la toute fin de la série (en 2018)
pipo.dist <- stock %>% pivot_wider(names_from = c("an", "point"), values_from = c("pred_CPUE")) %>%
  rename(avant.dedans = "2000_8", apres.dedans = "2018_8", avant.dehors = "2000_6", apres.dehors = "2018_6") %>%
  mutate(effet.dist.avant = (avant.dehors - avant.dedans)*100/avant.dedans,
         effet.dist.apres = (apres.dehors - apres.dedans)*100/apres.dedans) %>%
  dplyr::select(saison, espece, effet.dist.avant, effet.dist.apres)
stock3 <- full_join(pipo.an, pipo.dist)
stock2 <- rbind(stock2, stock3)
########################A bib-big plot
# homard
grob_sp <- grobTree(textGrob("Lobster", x=unit(1, "lines"),  y=unit(1, "lines"), hjust=.05, vjust= 1,
                             gp=gpar(col="grey25", fontsize=22, fontface="italic")))
grob_ltrA <- grobTree(textGrob("A", x=unit(1, "lines"),  y=unit(1, "lines"), hjust=-16, vjust = 1.2,
                               gp=gpar(col="grey25", fontsize=25)))
grob_ltrB <- grobTree(textGrob("B", x=unit(1, "lines"),  y=unit(1, "lines"), hjust=-16.2, vjust = 1.2,
                               gp=gpar(col="grey25", fontsize=25)))
gam_A <- plotGAM2(homard_gam, smooth.cov ="dist", groupCovs= "mois") + 
  theme(legend.position = "none") + geom_vline(xintercept = 500, linetype="dashed", size=1) +
  labs(x= "", y="", title= "") + annotation_custom(grob_sp) + annotation_custom(grob_ltrA)
gam_B <- plotGAM2(homard_gam, smooth.cov ="AN", groupCovs= "mois") + 
theme(legend.position = "none")+ annotation_custom(grob_ltrB) + 
labs(x= "", y="", title= "") 
# tourteau
grob_sp <- grobTree(textGrob("Edible crab", x=unit(1, "lines"),  y=unit(1, "lines"), hjust=.05, vjust= 1,
                             gp=gpar(col="black", fontsize=22, fontface="italic")))
grob_ltrC <- grobTree(textGrob("C", x=unit(1, "lines"),  y=unit(1, "lines"), hjust=-14.3, vjust = 1.2,
                               gp=gpar(col="grey25", fontsize=25)))
grob_ltrD <- grobTree(textGrob("D", x=unit(1, "lines"),  y=unit(1, "lines"), hjust=-14.5, vjust = 1.2,
                               gp=gpar(col="grey25", fontsize=25)))
gam_C <- plotGAM2(tourteau_gam, smooth.cov ="dist", groupCovs= "mois") + 
  theme(legend.position = "none") + geom_vline(xintercept = 500, linetype="dashed", size=1) + 
  labs(x= "", title= "") + annotation_custom(grob_sp) + annotation_custom(grob_ltrC) 
gam_D <- plotGAM2(tourteau_gam, smooth.cov ="AN", groupCovs= "mois") + 
  theme(legend.position = "none") + 
  labs(x= "", y="", title= "") + annotation_custom(grob_ltrD) 
# araignee
grob_sp <- grobTree(textGrob("Spider crab", x=unit(1, "lines"),  y=unit(1, "lines"), hjust=.05, vjust= 1,
                             gp=gpar(col="black", fontsize=22, fontface="italic")))
grob_ltrE <- grobTree(textGrob("E", x=unit(1, "lines"),  y=unit(1, "lines"), hjust=-16, vjust = 1.2,
                               gp=gpar(col="grey25", fontsize=25)))
grob_ltrF <- grobTree(textGrob("F", x=unit(1, "lines"),  y=unit(1, "lines"), hjust=-17.2, vjust = 1.2,
                               gp=gpar(col="grey25", fontsize=25)))
gam_E <- plotGAM2(araignee_gam, smooth.cov ="dist", groupCovs= "mois") + 
  theme(legend.position = "none") + geom_vline(xintercept = 500, linetype="dashed", size=1) + 
  labs(y="", x= "Distance", title= "") + annotation_custom(grob_sp) + annotation_custom(grob_ltrE) 
gam_F <- plotGAM2(araignee_gam, smooth.cov ="AN", groupCovs= "mois") + 
  theme(legend.position = "none") + 
  labs(y="", x= "Years", title= "") + annotation_custom(grob_ltrF) 
bigplot<- grid.arrange(gam_A, gam_B, gam_C, gam_D, gam_E, gam_F, ncol= 2)
# graphics.off()
# X11()
# bigplot
tiff("//portenbessin/serha/02_OBSERVATION_SURVEILLANCE/05_IGA/20_CNPE_FLAMANVILLE/56=ETUDES COMPLEMENTAIRES/2019_CANTOFLAM/article/figs/CPUE_gamm_221018.tif",
     units="in", width=10, height=10, res=300)
grid.arrange(gam_A, gam_B, gam_C, gam_D, gam_E, gam_F, ncol= 2)
dev.off()
