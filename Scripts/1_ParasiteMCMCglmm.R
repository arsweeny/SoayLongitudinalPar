#### 1 - MCMCglmm 

rm(list=ls())

library(MCMCglmm)
library(ggregplot)
library(RColorBrewer)
library(tidyverse)
library(grid)
library(INLA)
library(wesanderson)
library(nationalparkcolors)
library(janitor)
library(lubridate)

source("Functions/Model Effects Plot Tidy.R")


# set-up  -----------------------------------------------------------------

SeasonStrongPal <- brewer.pal(7, "Oranges")[3:7] 
SeasonCoccPal <- brewer.pal(7, "Purples")[3:7] 
SeasonCapPal <- brewer.pal(8, "Blues")[3:7]
SeasonNemPal <- brewer.pal(8, "Reds")[3:7]
SeasonStrongyloidesPal <- brewer.pal(8, "Greens")[3:7]
SeasonMoneziaPal <- brewer.pal(8, "PuRd")[3:7]


ParaCol1 <- wes_palette("Darjeeling1", n=5)
ParaCol2 <- park_palette("Badlands", n=5)
ParaCol3 <- park_palette("Arches", n=5)
ParaCol4 <- c(SeasonStrongPal[3], SeasonCoccPal[[3]], SeasonCapPal[3], SeasonNemPal[3])

SeasonPar <- list(SeasonStrongPal, SeasonCoccPal, SeasonCapPal, SeasonNemPal, 
                  SeasonStrongyloidesPal, SeasonMoneziaPal)

SeasonParAdult <- list(SeasonStrongPal, SeasonCoccPal, SeasonCapPal, #SeasonNemPal, 
                  SeasonStrongyloidesPal, SeasonMoneziaPal)

SeasonParLamb <- list(SeasonStrongPal, SeasonCoccPal, SeasonNemPal, # SeasonCapPal, 
                      SeasonStrongyloidesPal, SeasonMoneziaPal)

ParaPal <- SeasonPar %>% map(3) %>%  flatten() %>% unlist()


# Mods w/o Interaction  ----------------------------------------------------------

ParaModDF <- read.csv("Data/ParaModDF.csv")

ParaBin <- ParaFullYear %>% select(ends_with("Bin")) %>% colnames()

resps <- c("Strongyles", "Coccidia", "Capillaria", "Nematodirus")
respsLamb<- c("Strongyles", "Coccidia", "Nematodirus")
respsF <-paste(resps, "~")
respsFLamb <- respsF[c(1,2,4)]

respsBin <- paste(ParaBin)
respsBinF <- paste(ParaBin, "~")[6:7]

covar <- c("Season", "SexRepro", "AgeClass", "SampleTimeHour")


all <- ParaModDF %>% 
  mutate(Season=fct_relevel(Season, "Winter19", "Spring19", "Summer19","Autumn19", "Winter20"))
adults <- ParaModDF %>% filter(AgeClass=="Adult") %>% 
  mutate(Season=fct_relevel(Season, "Winter19", "Spring19", "Summer19","Autumn19", "Winter20")) 
adultsCap <- ParaModDF %>% filter(AgeClass=="Adult", Season %in% c("Winter19", "Spring19")) %>% 
  mutate(Season=fct_relevel(Season, "Winter19", "Spring19"))
adultsMon <- ParaModDF %>% filter(AgeClass=="Adult", Season %in% c("Spring19", "Summer19", "Autumn19")) %>% 
  mutate(Season=fct_relevel(Season, "Spring19", "Summer19","Autumn19", "Winter20"))
lambs <- ParaModDF %>% filter(AgeClass=="Lamb") %>% 
  mutate(Season=fct_relevel(Season, "Summer19","Autumn19", "Winter20"))
lambsStrong <- lambs %>% filter(Season!="Winter20") %>% 
  mutate(Season=fct_relevel(Season, "Summer19","Autumn19"))
  

ParaDFsAdult <- list(adults, adults, adultsCap)
ParaDFsAdultBin <- list(adults, adultsMon)
ParaDFsLambBin <- list(lambsStrong, lambs)

covar_all <- c(covar[c(1,3)], "Sex")  
covar_adult <- c(covar[1],"Sex")
covar_lamb <- c(covar[1], "Sex") 

intcovar_all <- c(covar_all, "SexRepro:Season")
intcovar_adult <- c(covar_adult[1], "SexRepro", "SexRepro:Season", "AgeScale")
intcovar_lamb <- c(covar_lamb, "Sex:Season")



# formulae  ---------------------------------------------------------------
FEC_F_Adult <- list()
for(x in 1:3) {
  FEC_F_Adult[[x]] <- as.formula(paste(respsF[x], paste(c(covar_adult), collapse = " + "))) 
} 

FEC_F_Lamb <- list()
for(x in 1:length(respsFLamb)) {
  FEC_F_Lamb[[x]] <- as.formula(paste(respsFLamb[x], paste(c(covar_lamb), collapse = " + "))) 
} 

BIN_F_Adult <- list()
for(x in 1:length(respsBinF)) {
  BIN_F_Adult[[x]] <- as.formula(paste(respsBinF[x], paste(c(covar_adult), collapse = " + "))) 
} 

BIN_F_Lamb<- list()
for(x in 1:length(respsBinF)) {
  BIN_F_Lamb[[x]] <- as.formula(paste(respsBinF[x], paste(c(covar_lamb), collapse = " + "))) 
}


FEC_F_Adult_Int <- list()
for(x in 1:2) {
  FEC_F_Adult_Int[[x]] <- as.formula(paste(respsF[x], paste(c(intcovar_adult), collapse = " + "))) 
} 
  
FEC_F_Lamb_Int <- list()
for(x in 1:2) {
  FEC_F_Lamb_Int[[x]] <- as.formula(paste(respsF[x], paste(c(intcovar_lamb), collapse = " + "))) 
}   

# Run Models --------------------------------------------------------------

Prior1 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100)))

Prior1b <- list(R = list(V = 10, fix = 1),
                G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1)*1000)))

Prior0b <- list(R = list(V = 10, fix = 1)) 


PriorMult <- list(R = list(V = diag(6), nu = 0, fix = 4), 
               G = list(G1 = list(V = diag(6), nu = 6)))

PriorMultB<-list(R=list(V=diag(6), nu=6.002), 
                    G=list(G1=list(V=diag(6), nu=6,alpha.mu=rep(0,6),alpha.V=diag(6)*100)))


IntensityModsAdult <- list() 
IntensityModsLamb <- list()

PrevModsAdult<-list()
PrevModsLamb<-list()

IntensityModsAdultInt <- list()


# Intensity Mods 
mf<-20

for(x in 1:3) {
  IntensityModsAdult[[x]] <- MCMCglmm(fixed = FEC_F_Adult[[x]] ,
                                 random =~ Id, 
                                 prior = Prior1,
                                 data = ParaDFsAdult[[x]],
                                 family = "poisson",
                                 pr=TRUE, 
                                 nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
                                 thin = 10*mf,burnin=3000*mf)
}

for(x in 1:3){
  saveRDS(IntensityModsAdult[[x]], paste0("Output/IntensityModSplitAdult", resps[x], ".rds"))
}

for(x in 1:3) {
  IntensityModsLamb[[x]] <- MCMCglmm(fixed = FEC_F_Lamb[[x]] ,
                                      random =~ Id, 
                                      prior = Prior1,
                                      data = lambs,
                                      family = "poisson",
                                      pr=TRUE, 
                                      nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
                                      thin = 10*mf,burnin=3000*mf)
}

for(x in 1:3){
  saveRDS(IntensityModsLamb[[x]], paste0("Output/IntensityModSplitLamb", respsLamb[x], ".rds"))
}


## quick look effects 
Efxplot(IntensityModsAdult[1:3], ModelNames = resps[1:3]) +
  theme_bw(base_size = 14) + 
  scale_colour_manual(values=ParaCol2) 

Efxplot(IntensityModsLamb[1:3], ModelNames = respsLamb[1:3]) +
  theme_bw(base_size = 14) + 
  scale_colour_manual(values=ParaCol2) 


# Prev Mods 

mf <-20 

for(x in 1:length(respsBinF)) {
  PrevModsAdult[[x]] <- MCMCglmm(fixed = BIN_F_Adult[[x]] ,
                            random =~ Id, 
                            prior = Prior1b,
                            data = ParaDFsAdultBin[[x]],
                            trunc=TRUE,
                            family = "categorical",
                            nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
                            thin = 10*mf,burnin=3000*mf)
}

for(x in 1:2){
  saveRDS(PrevModsAdult[[x]], paste0("Output/PrevModSplitAdult", respsBin[x+5], ".rds"))
}


for(x in 1:length(respsBinF)) {
  PrevModsLamb[[x]] <- MCMCglmm(fixed = BIN_F_Lamb[[x]] ,
                                 #random =~ Id, 
                                 prior = Prior0b,
                                 data = ParaDFsLambBin[[x]],
                                 trunc=TRUE,
                                 family = "categorical",
                                 nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
                                 thin = 10*mf,burnin=3000*mf)
}


for(x in 1:2){
  saveRDS(PrevModsLamb[[x]], paste0("Output/PrevModSplitLamb", respsBin[x+5], ".rds"))
}



##  quick look 
Efxplot(PrevModsAdult[1:2], ModelNames = ParaBin[6:7]) +
  theme_bw(base_size = 14) + 
  scale_colour_manual(values = ParaCol3) 
#monezis is p dodgy 

Efxplot(PrevModsLamb[1:2], ModelNames = ParaBin[6:7]) +
  theme_bw(base_size = 14) + 
  scale_colour_manual(values = ParaCol3) 


# adult models with interaction 

adults %<>% mutate(SexRepro=fct_relevel(SexRepro, "FemaleNoLamb", "FemaleLamb", "Male"))

for(x in 1:2) {
  IntensityModsAdultInt[[x]] <- MCMCglmm(fixed = FEC_F_Adult_Int[[x]] ,
                                      random =~ Id, 
                                      prior = Prior1,
                                      data = adults, 
                                      family = "poisson",
                                      pr=TRUE, 
                                      nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
                                      thin = 10*mf,burnin=3000*mf)
}

for(x in 1:2){
  saveRDS(IntensityModsAdultInt[[x]], paste0("Output/IntensityModAdultIntNoLambBase", resps[x], ".rds"))
}

Efxplot(IntensityModsAdultInt[1:2], ModelNames = resps[1:2]) +
  theme_bw(base_size = 14) + 
  scale_colour_manual(values=ParaCol4) 



MultiPriorPar<-list(R=list(V=diag(2), nu=2.002), 
                    G=list(G1=list(V=diag(2), nu=2,alpha.mu=rep(0,2),alpha.V=diag(2)*100)))



IntensityModsAdultInt[[3]] <-MCMCglmm(data = adults,
                              cbind(Strongyles, Coccidia) ~ trait -1 + 
                                trait:(Season + SexRepro + SexRepro:Season + AgeScale),
                              prior = MultiPriorPar,
                              random =~ us(trait):Id,
                              rcov =~ us(trait):units,
                              nitt = 13000*mf,
                              thin = 10*mf, burnin = 3000*mf,
                              family = c("poisson", "poisson")
)


saveRDS(IntensityModsAdultInt[[3]], "Output/IntensityModAdultIntBivarNoLambBase.rds")


summary.MCMCglmm(IntensityModsAdultInt[[3]])
plot(IntensityModsAdultInt[[3]])
ggregplot::PhenCorr(IntensityModsAdultInt[[3]], 2)
Efxplot(IntensityModsAdultInt[[3]])


