## 2 Model Outputs 

source("Functions/MCMC Effect Size Comparison.R")
source("Functions/MCMC Posterior Comparisons.R")
source("Functions/CleanMCMC.R")

ParaFullYear <- read.csv("ParasitesAnalysisDF_ARS.csv")

ParaBin <- ParaFullYear %>% select(ends_with("Bin")) %>% colnames()

resps <- c("Strongyles", "Coccidia", "Capillaria")
respsLamb<- c("Strongyles", "Coccidia", "Nematodirus")
respsBin <- paste(ParaBin)[6:7]

IntensitySplitAdult <- list.files("Output/", pattern= "IntensityModSplitAdult")[c(3,2,1)]
IntensitySplitLamb <- list.files("Output/", pattern= "IntensityModSplitLamb")[c(3,1,2)]

PrevSplitAdult <- list.files("Output/", pattern= "PrevModSplitAdult")[c(2,1)]
PrevSplitLamb <- list.files("Output/", pattern= "PrevModSplitLamb")[c(2,1)]

SplitAdult <- c(IntensitySplitAdult, PrevSplitAdult)
SplitLamb <- c(IntensitySplitAdult, PrevSplitLamb)

IntensityInt <- list.files("Output/", pattern = "IntensityModAdultInt")[c(6,3,1)]
IntensityInt2 <- list.files("Output/", pattern = "IntensityModAdultInt")[c(5,4,2)]

IntensityModsAdult <- list() 
IntensityModsLamb <- list()

PrevModsAdult<-list()
PrevModsLamb<-list()

IntensityModsAdultInt <- list()
IntensityModsAdultInt2 <- list()

# read in models for working with output  ---------------------------------

for(x in 1:3){
  IntensityModsAdult[[x]] <- readRDS(paste0("Output/IntensityModSplitAdult", resps[x], ".rds"))
}

for(x in 1:3){
  IntensityModsLamb[[x]] <- readRDS(paste0("Output/IntensityModSplitLamb", respsLamb[x], ".rds"))
}

for(x in 1:2){
  PrevModsAdult[[x]] <- readRDS(paste0("Output/PrevModSplitAdult", respsBin[x], ".rds"))
}

for(x in 1:2){
  PrevModsLamb[[x]] <- readRDS(paste0("Output/IntensityModSplitLamb", respsLamb[x], ".rds"))
}

for(x in 1:3){
  IntensityModsAdultInt[[x]] <- readRDS(paste0("Output/", IntensityInt[x]))
}

for(x in 1:3){
  IntensityModsAdultInt2[[x]] <- readRDS(paste0("Output/", IntensityInt2[x]))
}

RespAdult <- c(resps, respsBin)
RespLamb <- c(respsLamb, respsBin)
RespsInt <- c(resps[1:2], "Multi")

ModsAdult <- c(IntensityModsAdult, PrevModsAdult)
ModsLamb <-  c(IntensityModsLamb, PrevModsLamb)
ModsInt <- IntensityModsAdultInt


# set up aesthetics  ------------------------------------------------------
library(wesanderson); library(nationalparkcolors); library(RColorBrewer); library(ggregplot)


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
ParaCol5 <- c("#ffcce7", "#daf2dc", "81b7d2") # pink green blue 
ParaCol6 <- c("#33539e", "#7facd6", "#A195C9", "#c0b9db", "#e9b7d4", "#a5678e") # blue purple pinks 


SeasonPar <- list(SeasonStrongPal, SeasonCoccPal, SeasonCapPal, SeasonNemPal, 
                  SeasonStrongyloidesPal, SeasonMoneziaPal)

SeasonParAdult <- list(SeasonStrongPal, SeasonCoccPal, SeasonCapPal, #SeasonNemPal, 
                       SeasonStrongyloidesPal, SeasonMoneziaPal)

SeasonParLamb <- list(SeasonStrongPal, SeasonCoccPal, SeasonNemPal, # SeasonCapPal, 
                      SeasonStrongyloidesPal, SeasonMoneziaPal)


ParaPal <- SeasonPar %>% map(3) %>%  flatten() %>% unlist()

# tidy output -------------------------------------------------------------

ModelsOutAdult <- lapply(ModsAdult, clean.MCMC.GLMM) 
ModelsOutLamb <- lapply(ModsLamb, clean.MCMC.GLMM)
ModelsOutInt <- lapply(ModsInt, clean.MCMC.GLMM)

for(x in 1:length(ModelsOutAdult)){ 
  ModelsOutAdult[[x]] %>% 
    mutate(Response=RespAdult[x]) -> ModelsOutAdult[[x]]
}
for(x in 1:length(ModelsOutLamb)){ 
  ModelsOutLamb[[x]] %>% 
    mutate(Response=RespLamb[x]) -> ModelsOutLamb[[x]]
}
for(x in 1:length(ModelsOutInt)){ 
  ModelsOutInt[[x]] %>% 
    mutate(Response=RespsInt[x]) -> ModelsOutInt[[x]]
}

ModelsOutInt[[3]] <- ModelsOutInt[[3]] %>% as.data.frame() %>% 
  mutate(trait=str_split_fixed(as.character(variable), "[[:]]", 2)[,1], variable, 
         trait=str_remove(trait, "trait"), 
         term=str_split_fixed(as.character(variable), "[[:]]", 2)[,2], variable) %>% 
  select(-variable) %>% rename(variable=term)

outRound <- function(x){round(x, digits=2)}

ModelsTidyAdult <- bind_rows(ModelsOutAdult) %>% 
  rename(Estimate=post.mean, Lower=`l.95..CI`, Upper=`u.95..CI`, Term=variable) %>% 
  mutate_at(vars(Estimate, Lower, Upper), outRound) %>%
  mutate(EstimateReport= paste0(Estimate, " (", Lower, "-", Upper, ")")) %>% 
  select(Term, EstimateReport, eff.samp, pMCMC, Response)

ModelTidyLamb <- bind_rows(ModelsOutLamb) %>% 
  rename(Estimate=post.mean, Lower=`l.95..CI`, Upper=`u.95..CI`, Term=variable) %>% 
  mutate_at(vars(Estimate, Lower, Upper), outRound) %>%
  mutate(EstimateReport= paste0(Estimate, " (", Lower, "-", Upper, ")")) %>% 
  select(Term, EstimateReport, eff.samp, pMCMC, Response)

ModelTidyInt <- bind_rows(ModelsOutInt) %>% 
  rename(Estimate=post.mean, Lower=`l.95..CI`, Upper=`u.95..CI`, Term=variable) %>% 
  mutate_at(vars(Estimate, Lower, Upper), outRound) %>%
  mutate(EstimateReport= paste0(Estimate, " (", Lower, " to ", Upper, ")")) %>% 
  select(Term, EstimateReport, eff.samp, pMCMC, Response, trait)

write.csv(ModelsTidyAdult, "Output/AdultModsSplit.csv", row.names = F)
write.csv(ModelTidyLamb, "Output/LambModsSplit.csv", row.names = F)
write.csv(ModelTidyInt, "Output/AdultModsInt.csv", row.names = F)


# EffectPlots -------------------------------------------------------------

#Long Model Format Plot 

ModelOutLong ## lol what was this - bind all the outputa together i guess 

Starlocs <- ModelOutLong %>% group_by(Effect) %>% 
      summarise(max=max(Upper), 
                min=min(Lower), 
                Starloc = max+ min/ -10) %>% 
      select(Effect, Starloc)

ModelOutLong %>% 
  mutate(Sig = as.numeric(Lower*Upper>0)) %>% 
  mutate(Stars = case_when(pMCMC<0.001 ~ "***", pMCMC<0.01 ~ "**", pMCMC < 0.05 ~ "*", TRUE ~ "")) -> EffectsLong

EffectsLong %>% left_join(., Starlocs, by="Effect") %>% as.data.frame() %>% 
  mutate(Starloc=ifelse(Sig==0, NA, Starloc)) -> EffectsLong 

EffectsLong %>% 
  mutate(Effect=fct_relevel(Effect, "Intercept"), 
         Parasite=fct_relevel(Parasite, "Strongyles", "Coccidea", "Capillaria", 
                                          "Nematodirus", "Strongyloides", "Monezia")) %>% 
  filter(Parasite!="Trichuris") %>% 
  ggplot(aes(EffectLevel, Estimate, colour = Parasite)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.6) +
  geom_point(position = position_dodge(w = 0.6)) + 
  #coord_flip() +
  geom_errorbar(
    width = 0.3,
    position = position_dodge(w = 0.6),
    aes(ymin = Lower, ymax = Upper)) +
  geom_text(aes(label = sig, y = Starloc),
            position = position_dodge(w = 0.6), angle=45, vjust=1, show.legend=FALSE) +
  facet_wrap(~Effect, scales="free", ncol=4) +
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        strip.background = element_rect(fill="white"), 
        legend.position = "top") + 
  scale_colour_manual(limits=c("Strongyles", "Coccidea", "Capillaria", 
                               "Nematodirus", "Strongyloides", "Monezia"), 
                      values=ParaPal) -> LongModPlot 

LongModPlot + ggsave("Figures/AllParasiteModEfx.png", units="mm", height=180, width=350, dpi=400)


# Base Models  ------------------------------------------------------------

ParaDFsAdult <- list(adults, adults, adultsCap)
ParaDFsLamb <- list(lambs, lambs, lambs) 
ParaDFsAdultBin <- list(adults, adultsMon)
ParaDFsLambBin <- list(lambsStrong, lambs)

respsAdult <- resps[c(1:3)]
respsLamb <- resps[c(1,2,4)]

SummaryFECAdult<-
  lapply(1:length(respsAdult), function(x){
  require(plotrix)
  df <- ParaDFsAdult[[x]] %>% as.data.frame()
  df$Response <- df[, respsAdult[[x]]] 
  Summary <-list()
  df2 <- df %>% 
    mutate(logY=log(Response+1)) %>% 
    group_by(Season) %>% 
    summarise(mean=mean(logY),
              sd=sd(logY), 
              se=std.error(logY)) %>% as.data.frame() %>% 
    mutate(Response=resps[x])
  return(df2)
}) #%>% bind_rows() 

SummaryFECLamb<-
  lapply(1:length(respsLamb), function(x){
    require(plotrix)
    df <- ParaDFsLamb[[x]] %>% as.data.frame()
    df$Response <- df[, respsLamb[[x]]] 
    Summary <-list()
    df2 <- df %>% 
      mutate(logY=log(Response+1)) %>% 
      group_by(Season) %>% 
      summarise(mean=mean(logY),
                sd=sd(logY), 
                se=std.error(logY)) %>% as.data.frame() %>% 
      mutate(Response=resps[x])
    return(df2)
  }) #%>% bind_rows() 



modsBin <- respsBin[6:7]

SummariesBINAdult<-
  lapply(1:length(modsBin), function(x){
    require(plotrix)
    df <- ParaDFsAdultBin[[x]] %>% as.data.frame()
    df$Response <- df[, modsBin[[x]]] 
    Summary <-list()
    df2 <- df %>% 
      group_by(Season) %>% 
      summarise(n= length(Response), 
               mean=length(Response[Response==1])/ 
                  length(Response),
                se= sqrt(mean*(1-mean)/n)) %>% as.data.frame() %>% 
      mutate(Response=modsBin[x])
    return(df2)
  }) #%>% bind_rows() 

SummariesBINLamb<-
  lapply(1:length(modsBin), function(x){
    require(plotrix)
    df <- ParaDFsLambBin[[x]] %>% as.data.frame()
    df$Response <- df[, modsBin[[x]]] 
    Summary <-list()
    df2 <- df %>% 
      group_by(Season) %>% 
      summarise(n= length(Response), 
                mean=length(Response[Response==1])/ 
                  length(Response),
                se= sqrt(mean*(1-mean)/n)) %>% as.data.frame() %>% 
      mutate(Response=modsBin[x])
    return(df2)
  }) #%>% bind_rows() 

AdultSummaries <- do.call(c, list(SummaryFECAdult, SummariesBINAdult))
LambSummaries <- do.call(c, list(SummaryFECLamb, SummariesBINLamb))

RespsAdults <- c(paste0("log", resps[c(1:3)]), modsBin) 
RespsLambs <- c(paste0("log", respsLamb), modsBin) 

labsAdult <- c("log strongyle count", 
          "log coccidia count", 
          expression("log"~italic("Capillaria")~"count"), 
          #expression("log"~italic("Nematodirus")~"count"), 
          expression(~italic("Strongyloides")~"prevalence"), 
          expression(~italic("Moniezia")~"prevalence"))

labsLamb <- c("log strongyle count", 
              "log coccidia count", 
              #expression("log"~italic("Capillaria")~"count"), 
              expression("log"~italic("Nematodirus")~"count"), 
              expression(~italic("Strongyloides")~"prevalence"), 
              expression(~italic("Moniezia")~"prevalence"))

SeasonPlotsAdult <- list()
SeasonPlotsLamb <- list()

ParaDFsAdultAll <- do.call(c, list(ParaDFsAdult, ParaDFsAdultBin))
ParaDFsLambAll <- do.call(c, list(ParaDFsLamb, ParaDFsLambBin))


lapply(1:length(RespsAdults), function(x){
  require(ggforce)
  df<-ParaDFsAdultAll[[x]] %>% as.data.frame() 
  SummaryDF<- AdultSummaries[[x]]
  df$Response <- df[, RespsAdults[x]]
  ggplot() + 
    geom_sina(data=df, aes(y=Response, x=Season, colour=Season), alpha=0.24, maxwidth=0.5) + 
    geom_point(data=SummaryDF, aes(y=mean, x=Season), size=0.8) + 
    geom_errorbar(data=SummaryDF, aes(ymax=mean+se,ymin=mean-se, x=Season), 
                  width=0.3, size=0.6) + 
    theme_bw(base_size = 16) + 
    scale_colour_manual(values=SeasonParAdult[[x]])+ 
    labs(y=labsAdult[x]) + 
    theme(legend.position="none", 
          axis.text.x = (element_text(angle=45, hjust=1, colour="black")), 
          axis.title.x = element_blank())
}) -> SeasonPlotsAdult

lapply(1:length(RespsLambs), function(x){
  require(ggforce)
  df<-ParaDFsLambAll[[x]] %>% as.data.frame() 
  SummaryDF<- LambSummaries[[x]]
  df$Response <- df[, RespsLambs[x]]
  ggplot() + 
    geom_sina(data=df, aes(y=Response, x=Season, colour=Season), alpha=0.24, maxwidth=0.5) + 
    geom_point(data=SummaryDF, aes(y=mean, x=Season), size=0.8) + 
    geom_errorbar(data=SummaryDF, aes(ymax=mean+se,ymin=mean-se, x=Season), 
                  width=0.3, size=0.6) + 
    theme_bw(base_size = 16) + 
    scale_colour_manual(values=SeasonParLamb[[x]])+ 
    labs(y=labsLamb[x]) + 
    theme(legend.position="none", 
          axis.text.x = (element_text(angle=45, hjust=1, colour="black")), 
          axis.title.x = element_blank())
}) -> SeasonPlotsLamb

require(patchwork); require(gridExtra)

TopPanel <- SeasonPlotsAdult[[1]] + SeasonPlotsAdult[[2]] + SeasonPlotsAdult[[3]] + 
  SeasonPlotsAdult[[4]] + SeasonPlotsAdult[[5]]  + 
  plot_layout(nrow=1) #+ plot_annotation(tag_levels = "A") 

BottomPanel <- SeasonPlotsLamb[[1]] + SeasonPlotsLamb[[2]] + SeasonPlotsLamb[[3]] + 
  SeasonPlotsLamb[[4]] + SeasonPlotsLamb[[5]]  + 
  plot_layout(nrow=1) #+ plot_annotation(tag_levels = "A") 


TopPanel/ BottomPanel + plot_layout(nrow=2) + 
  plot_annotation(tag_levels = "A") -> SeasonPanelSplit 

#TopPanelComps/ BottomPanelComps + plot_layout(nrow=2) + 
 # plot_annotation(tag_levels = "A") -> SeasonPanelComps 


SeasonPanelSplit + ggsave('Figures/SeasonalityPanelSplit.png', units="mm", height=190, width=480, dpi=400)
SeasonPanelSplit + ggsave('Figures/SeasonalityPanelSplit.eps', device=cairo_ps,
                          units="mm", height=190, width=480, fallback_resolution=400)

#SeasonPanelComps + ggsave('Figures/SeasonPanelComps.png', units="mm", height=150, width=380, dpi=400)


# repeatability  ----------------------------------------------------------

AdultRep <- lapply(1:length(IntensityModsAdult[1:2]), function(a) {
  df <- MCMCRep(IntensityModsAdult[[a]], scale="link") %>% as.data.frame() %>% 
    mutate(Response=RespsInt[a])
  return(df)
}) %>% bind_rows()

LambRep <- lapply(1:length(IntensityModsLamb), function(a) {
  df <- MCMCRep(IntensityModsLamb[[a]], scale="link") %>% as.data.frame()
  return(df)
})  #basically all zero 

MultiVar <-  clean.MCMC.GLMM(IntensityModsAdultInt[[3]]) %>% as.data.frame() %>%
  filter(effect!="fixed") %>% 
  mutate(Component=
           case_when(str_detect(variable, "traitStrongyles:traitStrongyles|traitCoccidea:traitCoccidea") ~ "Id",
                     str_detect(variable, "traitStrongyles:traitCoccidea.Id|traitCoccidea:traitStrongyles.Id") ~ "between-individual",
                     str_detect(variable, "traitStrongyles:traitCoccidea.units|traitCoccidea:traitStrongyles.units") ~ "within-individual")) %>% 
  rename(lower=`l.95..CI`, upper= `u.95..CI`, estimate=post.mean) 

MCMCRep(IntensityModsAdultInt[[3]])

Model <- IntensityModsAdultInt[[3]] 
Sol <- Model$Sol %>% as.data.frame %>% select(contains("traitStrongyles")) %>% as.mcmc 
X <- Model$X[, seq(1,32,2)] 
VCV <- Model$VCV[, c(1,5)]  
Beta0 <- sapply(1:dim(Sol)[1],
                function(z) mean(as.matrix(X) %*% as.matrix(Sol[z, 1:ncol(X)])))
mat <- matrix(NA, nrow = dim(VCV)[2], ncol = 4)

RepStrong <- { for (j in 1:dim(VCV)[2]) {
  Va <- VCV[, j]
  Ve <- rowSums(VCV)
  Repeatability1 <- Va/(Ve + log(1/exp(Beta0) + 
                                   1))
  mat[j, ] <- c(colnames(VCV)[j], (posterior.mode(Repeatability1)), 
                (HPDinterval(Repeatability1)[, 1]), (HPDinterval(Repeatability1)[, 
                                                                                 2]))
} 
  colnames(mat) <- c("Component", "Mode", "lHPD", "uHPD")
data.frame(mat)
}  

Model <- IntensityModsAdultInt[[3]] 
Sol <- Model$Sol %>% as.data.frame %>% select(contains("traitCoccidia")) %>% as.mcmc 
X <- Model$X[, seq(2,32,2)] 
VCV <- Model$VCV[, c(4,8)]  
Beta0 <- sapply(1:dim(Sol)[1],
                function(z) mean(as.matrix(X) %*% as.matrix(Sol[z, 1:ncol(X)])))
mat <- matrix(NA, nrow = dim(VCV)[2], ncol = 4)
RepCocc <- { for (j in 1:dim(VCV)[2]) {
  Va <- VCV[, j]
  Ve <- rowSums(VCV)
  Repeatability1 <- Va/(Ve + log(1/exp(Beta0) + 
                                   1))
  mat[j, ] <- c(colnames(VCV)[j], (posterior.mode(Repeatability1)), 
                (HPDinterval(Repeatability1)[, 1]), (HPDinterval(Repeatability1)[, 
                                                                                 2]))
} 
  colnames(mat) <- c("Component", "Mode", "lHPD", "uHPD")
  data.frame(mat)
}  



AdultRepPlot <- 
AdultRep %>% filter(Component=="Id") %>% 
  mutate(Response=fct_relevel(Response, "Strongyles")) %>% 
  ggplot(aes(x=Response, y=Mode, colour=Response)) + 
  geom_point(size=5) + 
  geom_errorbar(aes(ymin=lHPD, ymax=uHPD), width=0.3, size=1.5) + 
  scale_colour_manual(values=ParaCol4) + theme_bw(base_size = 16) + labs(y="Repeatability") + 
  theme(legend.position = "none")

MultiRepPlot<-
MultiVar %>% filter(Component!="Id") %>% 
  ggplot(aes(x=Component, y=estimate)) + 
  geom_point(size=5, colour=ParaCol3[[5]]) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.3, size=1.5, colour=ParaCol3[[5]]) + 
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw(base_size = 16) + 
  labs(y="Strongyle-Coccidia Covariance")


ModelSet3Panel <- AdultRepPlot + MultiRepPlot + 
  plot_annotation(tag_levels = "A") 
ModelSet3Panel + ggsave("ModelSet3Panel.png", units="mm", height=120, width=240, dpi=300)
  
# posterior comps ---------------------------------------------------------

StrongSummary <- data_summary(ParaModDF %>% filter(AgeClass=="Adult"),
                              "logStrongyles", c("Season","SexRepro"))
require(ggforce)
StrongAdultPlotRaw <-
  ParaModDF %>% filter(AgeClass=="Adult") %>% #filter(Season!="Winter20") %>% 
  mutate(Season= fct_relevel(Season, "Winter19", "Spring19", "Summer19", "Autumn19")) %>% 
  ggplot() + 
  geom_sina(aes(x=Season, y=logStrongyles, colour=SexRepro), 
            maxwidth=0.7, position=position_dodge(0.6), alpha=0.16) + 
  geom_point(data=StrongSummary, size=2, 
             aes(x=Season, y=logStrongyles, colour=SexRepro), 
             position=position_dodge(0.6)) + 
  geom_errorbar(data=StrongSummary, width=0.3, size=1.2,
                aes(x=Season, ymax=logStrongyles + se, 
                    ymin=logStrongyles - se, colour=SexRepro), 
                position=position_dodge(0.6)) + 
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(colour="black")) + 
  scale_color_manual(values=ParaCol6[c(3,2,6)]) + 
  labs(y="log(strongyle FEC+1)")

StrongAdultPlotRaw + ggsave("Figures/FigureForReport.png", units="mm", height=180, width=200, dpi=400)


CoccSummary <- data_summary(ParaModDF %>% filter(AgeClass=="Adult"),
                            "logCoccidia", c("Season","SexRepro"))

CoccAdultPlotRaw <-
  ParaModDF %>% filter(AgeClass=="Adult") %>% 
  mutate(Season= fct_relevel(Season, "Winter19", "Spring19", "Summer19", "Autumn19")) %>%
  ggplot() + 
  geom_sina(aes(x=Season, y=logCoccidia, colour=SexRepro), 
            maxwidth=0.7, position=position_dodge(0.6), alpha=0.16) + 
  geom_point(data=CoccSummary, size=2, 
             aes(x=Season, y=logCoccidia, colour=SexRepro), 
             position=position_dodge(0.6)) + 
  geom_errorbar(data=CoccSummary, width=0.3, size=1.2,
                aes(x=Season, ymax=logCoccidia + se, 
                    ymin=logCoccidia - se, colour=SexRepro), 
                position=position_dodge(0.6)) + 
  theme_bw(base_size = 18) + 
  theme(legend.position = "bottom", 
        axis.text=element_text(colour="black")) + 
  scale_color_manual(values=ParaCol6[c(3,2,6)]) + 
  labs(y="log(coccidia FOC+1)")


#Strongyles 
Int<- IntensityModsAdultInt[[1]]$Sol %>% as.data.frame() %>%  select(contains(":")) %>% colnames()

IntEffects <- MCMCFxCompInt(IntensityModsAdultInt[[1]], Int)
IntEffects[[3]] <- MCMCFxComp(IntensityModsAdultInt[[1]], "Season")[[1]]  # use the main effect of season for the female lamb values???? 

IntEffects <-lapply(1:length(IntEffects), function(a) { 
  df<-IntEffects[[a]]
  diag(df$Mean)<-NA
  diag(df$Lower)<-NA
  diag(df$Upper)<-NA
  diag(df$pMCMC)<-NA
  return(df) 
})

IntMeans <- lapply(IntEffects, function(a) {
  df<-reshape2::melt(a$Mean)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

IntProbs <- lapply(IntEffects, function(a) {
  df <- reshape2::melt(a$pMCMC)
  df$Sig <- cut(df$value,breaks=c(-0.1,0.001,0.01,0.05,1),labels=c("***","**","*","")) 
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  df$Sig[is.na(df$Sig)]<-""
  return(df)
})

IntLower <- lapply(IntEffects, function(a) {
  df<-reshape2::melt(a$Lower)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

IntUpper <- lapply(IntEffects, function(a) {
  df<-reshape2::melt(a$Upper)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

IntSE <- lapply(1:length(IntEffects), function(a) {
  df <- reshape2::melt(IntEffects[[a]]$Lower)
  df$value <- round(df$value, 2)
  df$value2 <- round(reshape2::melt(IntEffects[[a]]$Upper)$value,2)
  df$value[is.na(df$value)]<-""
  df$value2[is.na(df$value2)]<-""
  df$intervals<-ifelse(df$value=="","", paste0("(",df$value,", ",df$value2,")"))
  return(df)
})


IntVals<-lapply(1:length(IntEffects), function(a) {
  CI<-IntSE[[a]] %>% select(Var1, Var2, value, value2) %>% 
    dplyr::rename(lower=value, upper=value2)
  Mean<-IntMeans[[a]] %>% select(value) %>% rename(Mean=value) 
  Sig<-IntProbs[[a]] %>% select(Sig) 
  df<-cbind(CI, Mean, Sig) %>% 
    mutate(Comp=paste(Var1, Var2, sep="-")) 
  return(df)
})


Names <- c("Female_NoLamb", "Male", "Female_Lamb") # note CHANGE THESE BACK if using original models 
IntComps<- MCMCPostCompInt(IntensityModsAdultInt[[1]], ParaModDF %>% filter(AgeClass=="A"), Int)
IntComps[[3]] <- MCMCPostComp(IntensityModsAdultInt[[1]], ParaModDF %>% filter(AgeClass=="A"), "Season")
IntComps<-lapply(1:length(IntComps), function(a) {
  df<-merge(IntComps[[a]], IntVals[[a]], by="Comp", all.x=TRUE)
})


CompDif<-setdiff(IntVals[[1]]$Comp, IntComps[[1]]$Comp)

Terms<- Names

CompExtreme<-lapply(1:length(IntVals), function(a) {
  df<-IntComps[[a]] %>% mutate(Term=Names[a]) %>% 
    group_by(Term) %>% 
    summarise(Max=max(Value),Min=min(Value))
  return(df)
})


axis_y<-c("SeasonSpring19-Intercept", #"SeasonSummer19-Intercept", "SeasonAutumn19-Intercept", "SeasonWinter20-Intercept", 
          "SeasonSummer19-SeasonSpring19", #"SeasonAutumn19-SeasonSpring19", "SeasonWinter20-SeasonSpring19", 
          "SeasonAutumn19-SeasonSummer19", #"SeasonWinter20-SeasonSummer19", 
          "SeasonWinter20-SeasonAutumn19"
          )

labels_y<-c("Spring19-Winter19", #"Summer19-Intercept", "Autumn19-Intercept", "Winter20-Intercept", 
            "Summer19-Spring19", #"Autumn19-Spring19", "Winter20-Spring19", 
            "Autumn19-Summer19", #"Winter20-Summer19", 
            "Winter20-Autumn19"
          )


IntVals2 <-lapply(1:length(IntEffects), function(a) {
  df<-IntVals[[a]] %>% filter(!Comp %in% CompDif) %>% 
    filter(Comp %in% axis_y)  %>%  # add to only do chronological comps 
    mutate(Comp=factor(Comp, levels=axis_y)) %>% 
    mutate(#CompNum=as.numeric(rev(Comp)), 
      CompNum=row_number() %>% rev, 
      lower=as.numeric(lower), 
      upper=as.numeric(upper)) %>% 
    mutate(Term=Terms[a])
  df<-merge(df, CompExtreme[[a]], by="Term", all.x=TRUE)
  return(df)
})

IntComps2<-bind_rows(IntComps) %>% as.data.frame() %>% 
  mutate(Term=case_when(Term=="SexReproFemaleNoLamb" ~Terms[1], # change back if OG models !! 
                        Term== "SexReproMale" ~ Terms[2], 
                        is.na(Term) ~ Terms[3])) 

IntVals3<-bind_rows(IntVals2) 


IntVals2[[1]]$Comp


IntPlotsRidge <- lapply(1:length(IntEffects),
                        function(a){
                          require(ggridges)  
                          #IntPlotsRidge[[a]]<-
                          IntComps[[a]] %>% filter(Comp %in% axis_y) %>% 
                            ggplot(aes(x=Value, y=Comp, fill=..x..)) +
                            geom_density_ridges_gradient(rel_min_height = 0.005)+
                            scale_fill_gradient2(high=ParaCol3[4], 
                                                 low=ParaCol3[5],  mid = "white",na.value = "grey", 
                                                 midpoint = 0, name="Post. Mean")+
                            geom_vline(xintercept=0, linetype="dashed", size=1.2)+
                            geom_segment(data=IntVals2[[a]],aes(y = CompNum-0.05, yend = CompNum+0.5,
                                                                x =lower , xend =lower), size=0.8, colour="gray60")+
                            geom_segment(data=IntVals2[[a]],aes(y = CompNum-0.05, yend = CompNum+0.5,
                                                                x =upper , xend =upper), size=0.8, colour="gray60")+
                            geom_text(data=IntVals2[[a]], aes(label = Sig, x=Max, y = CompNum + 0.2), size=9)+
                            theme_ridges(font_size = 20, grid=TRUE)+
                            theme(axis.title.x = element_text(hjust = 0.5),
                                  axis.title.y = element_text(hjust = 0.5)) + 
                            scale_y_discrete(limits=rev(axis_y), labels=rev(labels_y))+
                            labs(x="Posterior Mean Difference", y="Comparison Levels") + 
                            ggtitle(Terms[a])
                        }) 


IntVals3_Trim <- IntVals3 %>% filter(Term!="Female_Lamb")

IntPlotRidgeAll <-
  IntComps2 %>% filter(Term!="Female_Lamb") %>% filter(Comp %in% axis_y) %>% 
  ggplot(aes(x=Value, y=Comp, fill=..x..)) +
  geom_density_ridges_gradient(rel_min_height = 0.005)+
  #scale_fill_brewer(palette = "Reds")+
  scale_fill_gradient2(high=ParaCol3[4], 
                       low=ParaCol3[5],  mid = "white",na.value = "grey", 
                       midpoint = 0, name="Post. Mean")+
  geom_vline(xintercept=0, linetype="dashed", size=1)+
  geom_segment(data=IntVals3_Trim,aes(y = CompNum-0.05, yend = CompNum+0.5,
                                      x =lower , xend =lower), size=0.8, colour="gray40")+
  geom_segment(data=IntVals3_Trim,aes(y = CompNum-0.05, yend = CompNum+0.5,
                                      x =upper , xend =upper), size=0.8, colour="gray40")+
  geom_text(data=IntVals3_Trim, aes(label = Sig, x=Max+0.15, y = CompNum + 0.2), size=4.5)+
  theme_ridges(font_size = 14, grid=TRUE)+
  scale_y_discrete(limits=rev(axis_y), labels=rev(labels_y))+
  labs(x="Posterior Mean Difference", y="Comparison Levels")+
  facet_grid(~Term, scales="free_x")+
  theme(strip.background = element_rect(fill="white"), 
        strip.text = element_text(size=16), 
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        #axis.title.y = element_blank(),
        legend.position = "bottom") 

IntPlotRidgeAll ->  IntPlotRidgeStrong

# Coccidia 
Int<- IntensityModsAdultInt[[2]]$Sol %>% as.data.frame() %>%  select(contains(":")) %>% colnames()

IntEffects <- MCMCFxCompInt(IntensityModsAdultInt[[2]], Int)
IntEffects[[3]] <- MCMCFxComp(IntensityModsAdultInt[[2]], "Season")[[1]]  # use the main effect of season for the female lamb values???? 

IntEffects <-lapply(1:length(IntEffects), function(a) { 
  df<-IntEffects[[a]]
  diag(df$Mean)<-NA
  diag(df$Lower)<-NA
  diag(df$Upper)<-NA
  diag(df$pMCMC)<-NA
  return(df) 
})

IntMeans <- lapply(IntEffects, function(a) {
  df<-reshape2::melt(a$Mean)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

IntProbs <- lapply(IntEffects, function(a) {
  df <- reshape2::melt(a$pMCMC)
  df$Sig <- cut(df$value,breaks=c(-0.1,0.001,0.01,0.05,1),labels=c("***","**","*","")) 
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  df$Sig[is.na(df$Sig)]<-""
  return(df)
})

IntLower <- lapply(IntEffects, function(a) {
  df<-reshape2::melt(a$Lower)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

IntUpper <- lapply(IntEffects, function(a) {
  df<-reshape2::melt(a$Upper)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

IntSE <- lapply(1:length(IntEffects), function(a) {
  df <- reshape2::melt(IntEffects[[a]]$Lower)
  df$value <- round(df$value, 2)
  df$value2 <- round(reshape2::melt(IntEffects[[a]]$Upper)$value,2)
  df$value[is.na(df$value)]<-""
  df$value2[is.na(df$value2)]<-""
  df$intervals<-ifelse(df$value=="","", paste0("(",df$value,", ",df$value2,")"))
  return(df)
})


IntVals<-lapply(1:length(IntEffects), function(a) {
  CI<-IntSE[[a]] %>% select(Var1, Var2, value, value2) %>% 
    rename(lower=value, upper=value2)
  Mean<-IntMeans[[a]] %>% select(value) %>% rename(Mean=value) 
  Sig<-IntProbs[[a]] %>% select(Sig) 
  df<-cbind(CI, Mean, Sig) %>% 
    mutate(Comp=paste(Var1, Var2, sep="-")) 
  return(df)
})

Names <- c("Female_Lamb", "Male", "Female_NoLamb")
IntComps<- MCMCPostCompInt(IntensityModsAdultInt[[2]], ParaModDF %>% filter(AgeClass!="A"), Int)
IntComps[[3]] <- MCMCPostComp(IntensityModsAdultInt[[2]], ParaModDF %>% filter(AgeClass!="A"), "Season")
IntComps<-lapply(1:length(IntComps), function(a) {
  df<-merge(IntComps[[a]], IntVals[[a]], by="Comp", all.x=TRUE)
})



CompDif<-setdiff(IntVals[[1]]$Comp, IntComps[[1]]$Comp)

Terms<- Names

CompExtreme<-lapply(1:length(IntVals), function(a) {
  df<-IntComps[[a]] %>% mutate(Term=Names[a]) %>% 
    group_by(Term) %>% 
    summarise(Max=max(Value),Min=min(Value))
  return(df)
})

IntVals2 <-lapply(1:length(IntEffects), function(a) {
  df<-IntVals[[a]] %>% filter(!Comp %in% CompDif) %>% 
    filter(Comp %in% axis_y)  %>%  # add to only do chronological comps 
    mutate(Comp=factor(Comp, levels=axis_y)) %>% 
    mutate(#CompNum=as.numeric(rev(Comp)), 
           CompNum=row_number() %>% rev, 
           lower=as.numeric(lower), 
           upper=as.numeric(upper)) %>% 
    mutate(Term=Terms[a])
  df<-merge(df, CompExtreme[[a]], by="Term", all.x=TRUE)
  return(df)
})

IntComps2<-bind_rows(IntComps) %>% as.data.frame() %>% 
  mutate(Term=case_when(Term=="SexReproFemaleNoLamb" ~Terms[1], 
                        Term== "SexReproMale" ~ Terms[2], 
                        is.na(Term) ~ Terms[3])) 

IntVals3<-bind_rows(IntVals2) 


IntVals2[[1]]$Comp

axis_y<-c("SeasonSpring19-Intercept", #"SeasonSummer19-Intercept", "SeasonAutumn19-Intercept", "SeasonWinter20-Intercept", 
          "SeasonSummer19-SeasonSpring19", #"SeasonAutumn19-SeasonSpring19", "SeasonWinter20-SeasonSpring19", 
          "SeasonAutumn19-SeasonSummer19", #"SeasonWinter20-SeasonSummer19", 
          "SeasonWinter20-SeasonAutumn19"
)

labels_y<-c("Spring19-Winter19", #"Summer19-Intercept", "Autumn19-Intercept", "Winter20-Intercept", 
            "Summer19-Spring19", #"Autumn19-Spring19", "Winter20-Spring19", 
            "Autumn19-Summer19", #"Winter20-Summer19", 
            "Winter20-Autumn19"
)

IntPlotsRidge <- lapply(1:length(IntEffects),
                        function(a){
                          require(ggridges)  
                          #IntPlotsRidge[[a]]<-
                          IntComps[[a]] %>% filter(Comp %in% axis_y) %>% 
                            ggplot(aes(x=Value, y=Comp, fill=..x..)) +
                            geom_density_ridges_gradient(rel_min_height = 0.005)+
                            scale_fill_gradient2(high=ParaCol3[4], 
                                                 low=ParaCol3[5],  mid = "white",na.value = "grey", 
                                                 midpoint = 0, name="Post. Mean")+
                            geom_vline(xintercept=0, linetype="dashed", size=1.2)+
                            geom_segment(data=IntVals2[[a]],aes(y = CompNum-0.05, yend = CompNum+0.5,
                                                                x =lower , xend =lower), size=0.8, colour="gray60")+
                            geom_segment(data=IntVals2[[a]],aes(y = CompNum-0.05, yend = CompNum+0.5,
                                                                x =upper , xend =upper), size=0.8, colour="gray60")+
                            geom_text(data=IntVals2[[a]], aes(label = Sig, x=Max, y = CompNum + 0.2), size=9)+
                            theme_ridges(font_size = 20, grid=TRUE)+
                            theme(axis.title.x = element_text(hjust = 0.5),
                                  axis.title.y = element_text(hjust = 0.5)) + 
                            scale_y_discrete(limits=rev(axis_y), labels=rev(labels_y))+
                            labs(x="Posterior Mean Difference", y="Comparison Levels") + 
                          ggtitle(Terms[a])
                        }) 


IntVals3_Trim <- IntVals3 %>% filter(Term!="Female_Lamb")
IntPlotRidgeAll <-
  IntComps2 %>% filter(Term!="Female_Lamb") %>% 
  filter(Comp %in% axis_y) %>% 
  ggplot(aes(x=Value, y=Comp, fill=..x..)) +
  geom_density_ridges_gradient(rel_min_height = 0.005)+
  #scale_fill_brewer(palette = "Reds")+
  scale_fill_gradient2(high=ParaCol3[4], 
                       low=ParaCol3[5],  mid = "white",na.value = "grey", 
                       midpoint = 0, name="Post. Mean")+
  geom_vline(xintercept=0, linetype="dashed", size=1)+
  geom_segment(data=IntVals3_Trim,aes(y = CompNum-0.05, yend = CompNum+0.5,
                                 x =lower , xend =lower), size=0.8, colour="gray40")+
  geom_segment(data=IntVals3_Trim,aes(y = CompNum-0.05, yend = CompNum+0.5,
                                 x =upper , xend =upper), size=0.8, colour="gray40")+
  geom_text(data=IntVals3_Trim, aes(label = Sig, x=Max+0.15, y = CompNum + 0.2), size=4.5)+
  theme_ridges(font_size = 14, grid=TRUE)+
  scale_y_discrete(limits=rev(axis_y), labels=rev(labels_y))+
  labs(x="Posterior Mean Difference", y="Comparison Levels")+
  facet_grid(~Term, scales="free_x")+
  theme(strip.background = element_rect(fill="white"), 
        strip.text = element_text(size=16), 
        axis.title.x = element_text(hjust = 0.5),
        #axis.title.y = element_text(hjust = 0.5)
        axis.title.y = element_blank(), 
        legend.position = "bottom")

IntPlotRidgeAll -> IntPlotRidgeCocc

library(patchwork); library(ggpubr)

IntP <- IntPlots[[1]] + IntPlots[[2]]
IntPRidge <- IntPlotsRidge[[1]] + IntPlotsRidge[[3]] + IntPlotsRidge[[]] #weirdly not working - some fidgeting needed 
IntPRidge <- ggarrange(plotlist = IntPlotsRidge, nrow=1, ncol=3, common.legend = TRUE, legend = "right")


IntP + ggsave("ParasiteIntPlots.tiff", units="mm", height=100, width=260, dpi=300)
IntPlotRidgeAll + ggsave("ParasiteIntPlotRidge.tiff", units="mm", height=200, width=360, dpi=300)


(StrongAdultPlotRaw + CoccAdultPlotRaw) + 
  #plot_annotation(tag_levels = "A") + 
  plot_layout(guides = 'collect') & theme(legend.position = "bottom") -> ModelSet2PlotRaw

ggarrange(StrongAdultPlot, CoccAdultPlot, 
          labels=c("A", "B"), label.y=0.99, label.x = 0.01, 
          common.legend = TRUE, legend = "bottom") -> ModelSet2Plot

ggarrange(StrongAdultPlotRaw, CoccAdultPlotRaw, 
          labels=c("A", "B"), label.y=0.99, label.x = 0.01, 
          common.legend = TRUE, legend = "bottom") -> AdultIntPlotRaw

IntPRidge <- ggarrange(IntPlotRidgeStrong, IntPlotRidgeCocc, labels=c("C", "D"), 
                       label.y=0.99, label.x = 0.01)

AdultIntPlotRaw / IntPRidge  -> InteractionModPanel


InteractionModPanel + ggsave("Figures/InteractionModPanel.png", units="mm", height=380, width=380, dpi=500)

#  Relative change plots  ------------------------------------------------------

IntComps2 %<>% 
  mutate(TermComp = paste(Term, Comp, sep="_")) %>% 
  mutate_at(vars(lower, upper, Mean), as.numeric)

ggplot() + 
  geom_sina(data=IntComps2 %>% filter(Term!="Female_Lamb") %>%  
              filter(Comp %in% axis_y),
            aes(x=Comp, y=Value, colour=Term), 
            alpha=0.35, position=position_dodge(0.5)) + 
  geom_errorbar(data=IntComps2 %>% filter(Term!="Female_Lamb") %>%  
                  filter(!duplicated(TermComp)) %>% 
                  filter(Comp %in% axis_y),
                aes(x=Comp, ymax=upper, ymin=lower, group=Term), 
                colour="black", size=0.8, width=0.4, position=position_dodge(0.5)) +
  geom_point(data=IntComps2 %>% filter(Term!="Female_Lamb") %>%  
               filter(!duplicated(TermComp)) %>% 
               filter(Comp %in% axis_y),
             aes(x=Comp, y=Mean, group=Term), # pch=21,
             colour="black", size=3, position=position_dodge(0.5)) + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  scale_x_discrete(limits=axis_y, labels=labels_y) + 
  scale_y_continuous(breaks=seq(-3, 3, by=1)) + 
  scale_color_manual(values=ParaCol3[4:5]) + 
  geom_hline(yintercept=0, linetype="dashed") + 
  labs(y='Rel. change \n(compared to baseline females w/ lamb)', x='', 
       title="Interaction Posterior Comparisons") -> relChangeStrong


