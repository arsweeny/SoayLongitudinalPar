
library(MCMCglmm)

list.files("Output/", pattern="IntensityModAdultInt")[c(6,3,1)] -> modFile

IntensityModsAdultInt <- list()

for(i in 1:length(modFile)){
  IntensityModsAdultInt[[i]] <- readRDS(paste0("Output/", modFile[i]))
}



# Grand Means  ------------------------------------------------------------

Resp <- c("Strongyles", "Coccidia")
Covar <- c("Season", "SexRepro", "AgeScale", "Id")

TestDF <- 
  lapply(1:length(Resp), function(i){
    adults %>% 
      select(Resp[i], Covar) %>% 
      mutate(Id=as.integer(Id)) %>% 
      mutate(SexRepro = fct_relevel(SexRepro, "FemaleLamb", "FemaleNoLamb", "Male"))
  })

Model <- c(IntensityModsAdultInt[1:2])

Estimates <- 
  lapply(1:length(Model), function(x){
    Model[[x]]$Sol %>% as.matrix
  })

XMatrix <- lapply(1:length(TestDF), function(x){
  model.matrix(~ AgeScale, TestDF[[x]]) # Season + SexRepro + Season:SexRepro  - take out so it doesnt mess up adding these effects in 
})

map(XMatrix, colnames) -> xcol
map(Estimates, colnames) -> ecol

setdiff(xcol[[1]], ecol[[1]]) #fine 

PredictedIntercept <-list()

PredictedIntercept<-
  lapply(1:length(Estimates), function(x){
    est <- Estimates[[x]]
    xmat <- XMatrix[[x]]
  lapply(1:nrow(Estimates[[x]]), function(i){
    xmat %*% est[i,colnames(xmat)]
  })
})

PredictedInterceptList <- 
  lapply(1:length(PredictedIntercept), function(i){
    lapply(PredictedIntercept[[i]], mean) %>% unlist()
  })

#for(i in 1:nrow(Estimates)){
 # PredictedIntercept[[i]] <- XMatrix %*% Estimates[i,colnames(XMatrix)]
#} 

#PredictedIntercepts<- lapply(PredictedIntercept, mean) %>% unlist()



# predicted values  -------------------------------------------------------
newdat<- list()
newdat[[1]] <- expand.grid(Strongyles = mean(adults$Strongyles),
                      Season = as.factor(c("Winter19", "Spring19",
                                           "Summer19", "Autumn19", "Winter20")),
                      SexRepro = as.factor(c("FemaleLamb", "FemaleNoLamb", "Male")),
                      AgeScale =  mean(adults$AgeScale)) %>% as.data.frame() %>% 
  mutate(Season=fct_relevel(Season, "Winter19", "Spring19", "Summer19", "Autumn19", "Winter20"))

newdat[[2]] <- expand.grid(Coccidia = mean(adults$Coccidia),
                           Season = as.factor(c("Winter19", "Spring19",
                                                "Summer19", "Autumn19", "Winter20")),
                           SexRepro = as.factor(c("FemaleLamb", "FemaleNoLamb", "Male")),
                           AgeScale =  mean(adults$AgeScale)) %>% as.data.frame() %>% 
  mutate(Season=fct_relevel(Season, "Winter19", "Spring19", "Summer19", "Autumn19", "Winter20"))



XMatrix<- lapply(1:length(newdat), function(i){
  model.matrix(~ Season + SexRepro + Season:SexRepro+ AgeScale, newdat[[i]])
})

#XMatrix <- model.matrix(~ Season + SexRepro + Season:SexRepro+ AgeScale, newdat)

FocalColumns <- colnames(XMatrix[[1]]) %>% setdiff("(Intercept)") %>% setdiff("AgeScale")

PredictedVals<-list()
OldNames <- c(1:15) %>% as.character()
NewNames <- c("Intercept", FocalColumns)


PredictedVals<-
  lapply(1:length(Estimates), function(x){
    est <- Estimates[[x]]
    xmat <- XMatrix[[x]]
    lapply(1:nrow(Estimates[[x]]), function(i){
      xmat[,FocalColumns] %*% est[i,FocalColumns]
    })
  })


#PredictedVals<- 
#  lapply(1:nrow(Estimates), function(i){
#  XMatrix[,FocalColumns] %*% Estimates[i,FocalColumns]
#})  old code for one 

WideDF <- lapply(1:length(PredictedVals), function(i){
  PredictedVals[[i]] %>% bind_cols 
})

# WideDF <- PredictedVals %>% bind_cols - old code for one 

for(x in 1:length(WideDF)){
  colnames(WideDF[[x]]) <- paste0("I.", 1:1000)
}

#WideDF %<>% bind_cols(newdat, .)

WideDF <-
  lapply(1:length(WideDF), function(i){
  bind_cols(newdat[[i]], WideDF[[i]])
})

LongDF <-
  lapply(1:length(WideDF), function(i){
    WideDF[[i]] %>% gather("Iteration", "Value", -c(1:4))
})

#LongDF <- WideDF %>% gather("Iteration", "Value", -c(1:4))


PredictedValuesDF <- 
  lapply(1:length(LongDF), function(i){
    data.frame(Iteration = paste0("I.", 1:1000), PredictedInterceptList[[i]]) %>% 
      rename(PredictedIntercept=2) %>% 
      left_join(LongDF[[i]]) %>% mutate(Fit = PredictedIntercept + Value)
  })

#data.frame(Iteration = paste0("I.", 1:1000), PredictedIntercepts) %>% 
 # left_join(LongDF) %>% mutate(Fit = PredictedIntercepts + Value) -> 
#  PredictedValuesDF

PredictedValuesDF[[2]] %>% dim()

library(ggforce)

PredictedValuesDF[[2]] %>% 
  ggplot(aes(Season, (Fit), colour = SexRepro)) + 
  geom_sina(alpha = 0.3)

for(i in 1:length(PredictedValuesDF)){
  write.csv(PredictedValuesDF[[i]], paste0("Output/", Resp[i], "PredictedValuesDF.csv"))
}

# predicted comparisons  --------------------------------------------------

comp_trim<-c("SeasonSpring19-SeasonWinter19", #"SeasonSummer19-Intercept", "SeasonAutumn19-Intercept", "SeasonWinter20-Intercept", 
          "SeasonSummer19-SeasonSpring19", #"SeasonAutumn19-SeasonSpring19", "SeasonWinter20-SeasonSpring19", 
          "SeasonAutumn19-SeasonSummer19", #"SeasonWinter20-SeasonSummer19", 
          "SeasonWinter20-SeasonAutumn19"
) %>% str_remove("Season") %>% str_remove("Season")

labels_trim<-c("Winter19 to Spring19", #"Summer19-Intercept", "Autumn19-Intercept", "Winter20-Intercept", 
            "Spring19 to Summer19", #"Autumn19-Spring19", "Winter20-Spring19", 
            "Summer19 to Autumn19", #"Winter20-Summer19", 
            "Autumn19 to Winter20"
)

# functions ! 
PredictComp <-function(df, var, group){
  require(tidybayes)
  InteractionComps<-list()  
  df <- df 
  groupLevels <- levels(df[,group])
  df$group <- df[,group]
  df$var <- df[,var]
  enddf<-list()
  if(!is.null(group)){
    for(x in 1:length(groupLevels)){

      df1<-
        df %>% filter(group==groupLevels[[x]]) %>% 
        select(Iteration, var, Fit) %>% 
        spread(value="Fit", key="var") %>% 
        select(-Iteration)
      
      Names <-colnames(df1) 
      LevelCombs <- combn(ncol(df1),2) %>% t()
      
      df2 = apply(LevelCombs, 1, function(a){
        
        df3 = data.frame(Value = df1[,a[2]] - df1[,a[1]] ) %>% 
          mutate(Comp = paste(colnames(df1)[a[2]],colnames(df1)[a[1]], sep = "-")) 
        
      }) %>% bind_rows() 
      
      #df4 = 1:ncol(df1) %>% 
       # lapply(function(a) data.frame(Value = c(df1[,a]), Comp = paste(colnames(df1)[a],"Intercept", sep = "-"))) %>% 
       #  bind_rows()
      
      enddf[[x]] = df2 %>% 
        mutate(SexRepro = groupLevels[x]) %>% 
        mutate(Comp=as.factor(as.character(Comp)))
    }}
  
  enddf
}

CompComp <-function(df, var, group){

  df <- df 
  groupLevels <- levels(df[,group])
  df$group <- df[,group]
  df$var <- df[,var]
  iter <- 1:1000
  iter <- c(iter, iter, iter)
  FactorCompList<-list()
  
  enddf<-list()
  if(!is.null(group)){
    for(x in 1:length(groupLevels)){
      
      Columns <- df %>% 
        filter(group==groupLevels[[x]]) %>% 
        mutate(Iteration=iter) %>% 
        select(Iteration, Value, var) %>% 
        spread(value="Value", key="var") %>% 
        select(-Iteration)
      
      df1<-apply(Columns,2,function(b)
        apply(Columns,2,function(a) mean(a - b)))
      
      df1l<-apply(Columns,2,function(b)
        apply(Columns,2,function(a) HPDinterval(as.mcmc(a - b))[1]))
      
      df1u<-apply(Columns,2,function(b)
        apply(Columns,2,function(a) HPDinterval(as.mcmc(a - b))[2]))
      
      df1p<-apply(Columns,2,function(b) 
        apply(Columns,2,function(a) (dim(Columns)[1]-max(table(a>b)))/(dim(Columns)[1]/2)))
      
      FactorCompList[[x]] <- list(Mean=df1, Lower=df1l, Upper=df1u, pMCMC=df1p)
      names(FactorCompList)[x]<-groupLevels[x]

    }
  }
  FactorCompList
}

CompsSeason <- lapply(1:length(PredictedValuesDF), function(i){
  PredictComp(PredictedValuesDF[[i]], "Season", "SexRepro")
})


# do the stuff 
CompsSeasonDF <- map(CompsSeason, bind_rows)

CompsSeasonSex <- list()

Effects <- lapply(1:length(CompsSeasonDF), function(i){
  CompComp(CompsSeasonDF[[i]], "SexRepro", "Comp")
})

Effects <- Effects[[1]]

Effects <-lapply(1:length(Effects), function(a) { 
  df<-Effects[[a]]
  diag(df$Mean)<-NA
  diag(df$Lower)<-NA
  diag(df$Upper)<-NA
  diag(df$pMCMC)<-NA
  return(df) 
})

Means <- lapply(Effects, function(a) {
  df<-reshape2::melt(a$Mean)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

Probs <- lapply(Effects, function(a) {
  df <- reshape2::melt(a$pMCMC)
  df$Sig <- cut(df$value,breaks=c(-0.1,0.001,0.01,0.05,1),labels=c("***","**","*","")) 
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  df$Sig[is.na(df$Sig)]<-""
  return(df)
})

Lower <- lapply(Effects, function(a) {
  df<-reshape2::melt(a$Lower)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

Upper <- lapply(Effects, function(a) {
  df<-reshape2::melt(a$Upper)
  df$value <- round(df$value, 2)
  df$value[is.na(df$value)]<-""
  return(df)
})

SE <- lapply(1:length(Effects), function(a) {
  df <- reshape2::melt(Effects[[a]]$Lower)
  df$value <- round(df$value, 2)
  df$value2 <- round(reshape2::melt(Effects[[a]]$Upper)$value,2)
  df$value[is.na(df$value)]<-""
  df$value2[is.na(df$value2)]<-""
  df$ervals<-ifelse(df$value=="","", paste0("(",df$value,", ",df$value2,")"))
  return(df)
})

groupLevels <- levels(CompsSeasonDF[[1]][,"Comp"])

Vals<-lapply(1:length(Effects), function(a) {
  CI<-SE[[a]] %>% select(Var1, Var2, value, value2) %>% 
    rename(lower=value, upper=value2)
  Mean<-Means[[a]] %>% select(value) %>% rename(Mean=value) 
  pVal = Probs[[a]] %>% select(value)  
  Sig<-Probs[[a]] %>% select(Sig) 
  df<-cbind(CI, Mean, pVal, Sig) %>% 
    mutate(SexComp=paste(Var1, Var2, sep="-"),
           SeasonComp=groupLevels[a]) 
  return(df)
})

names(Vals) <- groupLevels

Vals2 <- bind_rows(Vals) %>% 
  filter(lower!='') %>% 
  mutate(value=ifelse(value==0, "<0.001", value))

Vals2 %>% filter(SeasonComp %in% comp_trim) %>% View()

write.csv(Vals2, "Output/StrongylesRelChangeComps.csv", row.names = F)
write.csv(Vals2, "Output/CoccidiaRelChangeComps.csv", row.names = F)



## plot the comps ------------ 
CompSum <- lapply(1:2, function(i){
  CompsSeasonDF[[i]] %>% 
  filter(Comp %in% comp_trim) %>% 
  group_by(Comp, SexRepro) %>% 
  mutate(Comp=factor(Comp)) %>% 
  summarise(mean=mean(Value),
         lower=HPDinterval(as.mcmc(Value))[1],
         upper=HPDinterval(as.mcmc(Value))[2]) %>% as.data.frame() %>% 
  mutate(Comp=fct_relevel(Comp, comp_trim)) %>% 
  mutate(CompNum=case_when(Comp==comp_trim[1] ~ 4, 
                           Comp==comp_trim[2] ~ 3,
                           Comp==comp_trim[3] ~ 2,
                           Comp==comp_trim[4] ~ 1)) 
})

PredSum <-lapply(1:2, function(i){
  PredictedValuesDF[[i]] %>% 
  group_by(Season, SexRepro) %>% 
  summarise(mean=mean(Fit),
            lower=HPDinterval(as.mcmc(Fit))[1],
            upper=HPDinterval(as.mcmc(Fit))[2]) 
}) 


ggplot() + 
  geom_sina(data=CompsSeasonDF[[1]] %>% 
              filter(Comp %in% comp_trim),
            aes(x=Comp, y=Value, colour=SexRepro), 
            alpha=0.1, position=position_dodge(0.5)) + 
  geom_errorbar(data=CompSum[[1]] %>% 
                  filter(Comp %in% comp_trim),
                aes(x=Comp, ymax=upper, ymin=lower, group=SexRepro), 
                colour="black", size=0.8, width=0.4, position=position_dodge(0.5)) +
  geom_point(data=CompSum[[1]] %>% 
               filter(Comp %in% comp_trim),
             aes(x=Comp, y=mean, group=SexRepro), # pch=21,
             colour="black", size=3, position=position_dodge(0.5)) + 
  theme_bw(base_size = 18) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, colour="black"),
      legend.position = "none") + 
  scale_x_discrete(limits=comp_trim, labels=labels_trim) + 
  scale_y_continuous(breaks=seq(-3, 3, by=1)) + 
  scale_color_manual(values=ParaCol6[c(3,2,6)]) + 
  geom_hline(yintercept=0, linetype="dashed") + 
  labs(y='Relative change in log strongyles', x='') -> StrongylePredChange 

ggplot() + 
  geom_sina(data=CompsSeasonDF[[2]] %>% 
              filter(Comp %in% comp_trim),
            aes(x=Comp, y=Value, colour=SexRepro), 
            alpha=0.1, position=position_dodge(0.5)) + 
  geom_errorbar(data=CompSum[[2]] %>% 
                  filter(Comp %in% comp_trim),
                aes(x=Comp, ymax=upper, ymin=lower, group=SexRepro), 
                colour="black", size=0.8, width=0.4, position=position_dodge(0.5)) +
  geom_point(data=CompSum[[2]] %>% 
               filter(Comp %in% comp_trim),
             aes(x=Comp, y=mean, group=SexRepro), # pch=21,
             colour="black", size=3, position=position_dodge(0.5)) + 
  theme_bw(base_size = 18) + 
  theme(axis.text.x = element_text(angle=45, hjust=1, colour="black"), 
        legend.position = "none") + 
  scale_x_discrete(limits=comp_trim, labels=labels_trim) + 
  scale_y_continuous(breaks=seq(-3, 3, by=1)) + 
  scale_color_manual(values=ParaCol6[c(3,2,6)]) + 
  geom_hline(yintercept=0, linetype="dashed") + 
  labs(y='Relative change in log coccidia', x='') -> CoccidiaPredChange 


library(ggpubr); library(patchwork)
ggarrange(StrongAdultPlotRaw, CoccAdultPlotRaw, 
          labels=c("A", "B"), label.y=0.99, label.x = 0.01, 
          common.legend = TRUE, legend = "bottom") -> AdultIntPlotRaw

IntPred <- ggarrange(StrongylePredChange, CoccidiaPredChange, labels=c("C", "D"), 
                       label.y=0.99, label.x = 0.01)

AdultIntPlotRaw / IntPred  -> InteractionModPanel


InteractionModPanel + ggsave("Figures/InteractionModPanelPred.png", units="mm", height=350, width=350, dpi=500)

