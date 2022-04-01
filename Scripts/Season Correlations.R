
rm(list=ls())

list.files("Output/", pattern="IntensityModSplitAdult")[c(3,2,1)] -> IntensityNamesAdult 
list.files("Output/", pattern="IntensityModSplitLamb")[c(3,1,2)] -> IntensityNamesLamb 
list.files("Output/", pattern="PrevModSplit")[c(2,1)] -> PrevNamesAdult 
list.files("Output/", pattern="PrevModSplit")[c(4,3)] -> PrevNamesLamb 

IntensityAdult <- list()
IntensityLamb <-  list()
PrevAdult <- list()
PrevLamb <- list()

for(x in 1:3){
  IntensityAdult[[x]] <- readRDS(paste0("Output/", IntensityNamesAdult[x]))
}

for(x in 1:3){
  IntensityLamb[[x]] <- readRDS(paste0("Output/", IntensityNamesLamb[x]))
}


for(x in 1:2){
  PrevAdult[[x]] <- readRDS(paste0("Output/", PrevNamesAdult[x]))
}

for(x in 1:2){
  PrevLamb[[x]] <- readRDS(paste0("Output/", PrevNamesLamb[x]))
}


TermsAdult <- 
  IntensityAdult[[1]]$Sol %>% as.data.frame() %>%select(contains("Season")) %>% 
  colnames()

TermsLamb <- 
  IntensityLamb[[1]]$Sol %>% as.data.frame() %>%select(contains("Season")) %>% 
  colnames()

ModelListAdult <-  do.call(c, list(IntensityAdult[1:3], PrevAdult))
ModelListLamb <-  do.call(c, list(IntensityLamb[1:3], PrevLamb))

SeasonEffectsAdult <- list()
SeasonEffectsLamb <- list()


## Adult comps - need to reload new function 
for(x in 1:length(ModelListAdult)){ 
 SeasonEffectsAdult[[x]] <- MCMCFxComp(ModelListAdult[[x]], "Season")[[1]]
} 

lapply(SeasonEffectsAdult[[3]], as.matrix) -> SeasonEffectsAdult[[3]] 
CapOverwrite <- SeasonEffectsAdult[[3]]

SeasonEffectsAdult <-lapply(1:length(SeasonEffectsAdult), function(a) { 
  df<-SeasonEffectsAdult[[a]]
  diag(df$Mean)<-NA
  diag(df$Lower)<-NA
  diag(df$Upper)<-NA
  diag(df$pMCMC)<-NA
  #df$Mean[upper.tri(df$Mean)]<-NA
  #df[upper.tri(df$Mean)]<-NA
  return(df) 
})

SeasonEffectsAdult[[3]] <- CapOverwrite


SeasonEffectsAdultLower <-lapply(1:length(SeasonEffectsAdult), function(a) { 
  df<-SeasonEffectsAdult[[a]]
  df$Mean[lower.tri(df$Mean, diag=FALSE)]<-NA
  df$Lower[lower.tri(df$Lower, diag=FALSE)]<-NA
  df$Upper[lower.tri(df$Upper, diag=FALSE)]<-NA
  df$pMCMC[lower.tri(df$pMCMC, diag=FALSE)]<-NA
  diag(df$Mean)<-NA
  diag(df$Lower)<-NA
  diag(df$Upper)<-NA
  diag(df$pMCMC)<-NA
  return(df) 
})

SeasonEffectsAdultLower[[3]] <- CapOverwrite


SeasonMeansAdult <- lapply(SeasonEffectsAdultLower, function(a) {
  df<-reshape2::melt(a$Mean)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

SeasonProbsAdult <- lapply(SeasonEffectsAdultLower, function(a) {
  df <- reshape2::melt(a$pMCMC)
  df$Sig <- cut(df$value,breaks=c(-0.1,0.001,0.01,0.05,1),labels=c("***","**","*","")) 
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  df$Sig[is.na(df$Sig)]<-""
  return(df)
})

SeasonLowerAdult <- lapply(SeasonEffectsAdultLower, function(a) {
  df<-reshape2::melt(a$Lower)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

SeasonUpperAdult <- lapply(SeasonEffectsAdultLower, function(a) {
  df<-reshape2::melt(a$Upper)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

SeasonSEAdult <- lapply(1:length(SeasonEffectsAdultLower), function(a) {
  df <- reshape2::melt(SeasonEffectsAdultLower[[a]]$Lower)
  df$value <- round(df$value, 2)
  df$value2 <- round(reshape2::melt(SeasonEffectsAdultLower[[a]]$Upper)$value,2)
  df$value[is.na(df$value)]<-""
  df$value2[is.na(df$value2)]<-""
  df$intervals<-ifelse(df$value=="","", paste0("(",df$value,", ",df$value2,")"))
  return(df)
})


SeasonValsAdult<-lapply(1:length(SeasonEffectsAdult), function(a) {
  CI<-SeasonSEAdult[[a]] %>% select(Var1, Var2, value, value2) %>% 
    rename(lower=value, upper=value2)
  Mean<-SeasonMeansAdult[[a]] %>% select(value) %>% rename(Mean=value) 
  Sig<-SeasonProbsAdult[[a]] %>% select(Sig) 
  df<-cbind(CI, Mean, Sig) %>% 
    mutate(Comp=paste(Var1, Var2, sep="-")) 
  return(df)
})



SeasonPlotsCompAdult <- lapply(1:length(SeasonEffectsAdult),
                    function(a){
                      ggplot(SeasonMeansAdult[[a]], aes(Var1, Var2)) +
                        geom_tile(aes(fill = as.numeric(value))) +
                        coord_fixed() + theme_bw(base_size = 10)+#
                        geom_text(aes(label = value), size=2.4)+
                        geom_text(data = SeasonProbsAdult[[a]],
                                  aes(label = Sig, y = as.numeric(rev(Var2)) + 0.2)) +
                        #geom_text(data = SeasonSE[[a]],
                         #         aes(label = intervals, y = as.numeric(rev(Var2))- 0.3),
                        #          size = 1) +
                        scale_fill_gradient2(mid="white",   low=SeasonParAdult[[a]][4], 
                                            high=SeasonParAdult[[a]][4], na.value = "white")+  # midpoint = 0)+
                        xlab(NULL) + ylab(NULL) + theme(legend.position="none") +
                        #ggtitle(paste0(names(FullModelMeansXS[a]), ", cross-sectional")) +
                        theme(axis.text.x = element_text(angle=45, hjust=1, colour="black"),
                              axis.text.y= element_text(colour="black"))+
                        #theme(axis.text.y = element_text(size=8,angle=90,hjust=0.5)) +
                        scale_y_discrete(limits = rev(levels(SeasonMeansAdult[[a]][,1])))
                      #scale_x_discrete(limits = levels(FullModelXSMeans[[a]][,1]))
                      
                             })





## lamb  comps - need to reload new function 
for(x in 1:length(ModelListLamb)){ 
  SeasonEffectsLamb[[x]] <- MCMCFxComp(ModelListLamb[[x]], "Season")[[1]]
} 

lapply(SeasonEffectsLamb[[4]], as.matrix) -> SeasonEffectsLamb[[4]] 
StrongyloidesOverwrite <- SeasonEffectsLamb[[4]]

SeasonEffectsLamb <-lapply(1:length(SeasonEffectsLamb), function(a) { 
  df<-SeasonEffectsLamb[[a]]
  diag(df$Mean)<-NA
  diag(df$Lower)<-NA
  diag(df$Upper)<-NA
  diag(df$pMCMC)<-NA
  #df$Mean[upper.tri(df$Mean)]<-NA
  #df[upper.tri(df$Mean)]<-NA
  return(df) 
})

SeasonEffectsLamb[[4]] <- StrongyloidesOverwrite


SeasonEffectsLambLower <-lapply(1:length(SeasonEffectsLamb), function(a) { 
  df<-SeasonEffectsLamb[[a]]
  df$Mean[lower.tri(df$Mean, diag=FALSE)]<-NA
  df$Lower[lower.tri(df$Lower, diag=FALSE)]<-NA
  df$Upper[lower.tri(df$Upper, diag=FALSE)]<-NA
  df$pMCMC[lower.tri(df$pMCMC, diag=FALSE)]<-NA
  diag(df$Mean)<-NA
  diag(df$Lower)<-NA
  diag(df$Upper)<-NA
  diag(df$pMCMC)<-NA
  return(df) 
})

SeasonEffectsLambLower[[4]] <- StrongyloidesOverwrite


SeasonMeansLamb <- lapply(SeasonEffectsLambLower, function(a) {
  df<-reshape2::melt(a$Mean)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

SeasonProbsLamb <- lapply(SeasonEffectsLambLower, function(a) {
  df <- reshape2::melt(a$pMCMC)
  df$Sig <- cut(df$value,breaks=c(-0.1,0.001,0.01,0.05,1),labels=c("***","**","*","")) 
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  df$Sig[is.na(df$Sig)]<-""
  return(df)
})

SeasonLowerLamb <- lapply(SeasonEffectsLambLower, function(a) {
  df<-reshape2::melt(a$Lower)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

SeasonUpperLamb <- lapply(SeasonEffectsLambLower, function(a) {
  df<-reshape2::melt(a$Upper)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

SeasonSELamb <- lapply(1:length(SeasonEffectsLambLower), function(a) {
  df <- reshape2::melt(SeasonEffectsLambLower[[a]]$Lower)
  df$value <- round(df$value, 2)
  df$value2 <- round(reshape2::melt(SeasonEffectsLambLower[[a]]$Upper)$value,2)
  df$value[is.na(df$value)]<-""
  df$value2[is.na(df$value2)]<-""
  df$intervals<-ifelse(df$value=="","", paste0("(",df$value,", ",df$value2,")"))
  return(df)
})


SeasonValsLamb<-lapply(1:length(SeasonEffectsLamb), function(a) {
  CI<-SeasonSELamb[[a]] %>% select(Var1, Var2, value, value2) %>% 
    rename(lower=value, upper=value2)
  Mean<-SeasonMeansLamb[[a]] %>% select(value) %>% rename(Mean=value) 
  Sig<-SeasonProbsLamb[[a]] %>% select(Sig) 
  df<-cbind(CI, Mean, Sig) %>% 
    mutate(Comp=paste(Var1, Var2, sep="-")) 
  return(df)
})



SeasonPlotsCompLamb <- lapply(1:length(SeasonEffectsLamb),
                               function(a){
                                 ggplot(SeasonMeansLamb[[a]], aes(Var1, Var2)) +
                                   geom_tile(aes(fill = as.numeric(value))) +
                                   coord_fixed() + theme_bw(base_size = 10)+#
                                   geom_text(aes(label = value), size=2.4)+
                                   geom_text(data = SeasonProbsLamb[[a]],
                                             aes(label = Sig, y = as.numeric(rev(Var2)) + 0.2)) +
                                   #geom_text(data = SeasonSE[[a]],
                                   #         aes(label = intervals, y = as.numeric(rev(Var2))- 0.3),
                                   #          size = 1) +
                                   scale_fill_gradient2(mid="white",   low=SeasonParLamb[[a]][4], 
                                                        high=SeasonParLamb[[a]][4], na.value = "white")+  # midpoint = 0)+
                                   xlab(NULL) + ylab(NULL) + theme(legend.position="none") +
                                   #ggtitle(paste0(names(FullModelMeansXS[a]), ", cross-sectional")) +
                                   theme(axis.text.x = element_text(angle=45, hjust=1, colour="black"),
                                         axis.text.y= element_text(colour="black"))+
                                   #theme(axis.text.y = element_text(size=8,angle=90,hjust=0.5)) +
                                   scale_y_discrete(limits = rev(levels(SeasonMeansLamb[[a]][,1])))
                                 #scale_x_discrete(limits = levels(FullModelXSMeans[[a]][,1]))
                                 
                               })


require(patchwork)

SeasonPlotsCompAdult[[1]] + SeasonPlotsCompAdult[[2]] + 
  SeasonPlotsCompAdult[[4]] + SeasonPlotsCompAdult[[5]]  


SeasonPlotsCompLamb[[1]] + SeasonPlotsCompLamb[[2]] + 
  SeasonPlotsCompLamb[[3]] + SeasonPlotsCompLamb[[5]]  
