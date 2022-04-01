
## ggplot2 theme for effects plot with gridded background 
THEMEGRID<-theme(axis.text.x=element_text(size=12,colour="black"),
                 axis.text.y=element_text(size=12,colour="black"))+
  theme(axis.title.x=element_text(vjust=-0.35),axis.title.y=element_text(vjust=1.2))+theme_bw()



## functions to construct forest plot of effect sizes 
## applicable to glmmTMB and lme4 and similar classes of models 

# Vertical Forest Plot 
EfxplotTidy<-function(modellist, sig=TRUE, ModelNames = NULL, relabels=NULL, submodelNames=NULL){
  require(ggplot2)
  require(dplyr)
  require(broom.mixed)
  require(dotwhisker)
  graphlist<-list()
  for(i in 1:length(modellist)){
    model<-modellist[[i]]
    graph<-tidy(model, conf.int=TRUE, effects="fixed") %>% 
      #relabel_predictors(relabels) %>% 
      mutate(Sig=ifelse(p.value<0.0011, "***", ifelse(p.value>0.001&p.value<0.011, "**", 
                                                 ifelse( p.value<0.05&p.value>0.01, "*", ""))))
    values<-c("estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
    graph[,values]<-
      mapply(round, digits=2, graph[, values])
    
    graph$Model<-i
    
    if(!is.null(submodelNames)){     #useful if aesthetics according to a submodel group desired later 
      graph$submodel<-submodelNames[[i]]
    }
    
    graphlist[[i]]<-graph
    
    
  }
  
  graph<-bind_rows(graphlist)
  
  graph$starloc<-NA   
  
  min<-min(graph$conf.low,na.rm=T)
  max<-max(graph$conf.high,na.rm=T)
  
  if(sig==TRUE){
    graph$starloc <- max+(max-min)/10
  }
  
  
  graph$Model<-as.factor(graph$Model)
  
  if(!is.null(ModelNames)){
    levels(graph$Model)<-ModelNames
  }

  
  
  ggplot(as.data.frame(graph),aes(x=term,y=estimate,colour=Model))+
    geom_point(position=position_dodge(w=0.5), size=2)+
    geom_errorbar(position=position_dodge(w=0.5),aes(ymin=conf.low,ymax=conf.high),size=0.3,width=0.2)+
    geom_hline(aes(yintercept=0),lty=2)+THEMEGRID+labs(x=NULL)+coord_flip()+
    geom_text(aes(label=Sig,y=starloc),size=5, position=position_dodge(w=0.5), 
              show.legend = FALSE)+
    theme(strip.text.x = element_text(size = 12), 
          strip.background =element_rect(fill="white"))
  
}

# horizontal effect sizes 
EfxplotWide<-function(modellist, sig=TRUE, ModelNames = NULL, submodelNames=NULL, relabels=NULL){
  require(ggplot2)
  require(broom.mixed)
  graphlist<-list()
  for(i in 1:length(modellist)){
    model<-modellist[[i]]
    graph<-tidy(model, conf.int=TRUE) %>% 
      relabel_predictors(relabels) %>% 
      mutate(Sig=ifelse(p.value<0.0011, "***", ifelse(p.value>0.001&p.value<0.011, "**", 
                                                      ifelse( p.value<0.05&p.value>0.01, "*", ""))))
    values<-c("estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
    graph[,values]<-
      mapply(round, digits=2, graph[, values])
    
    graph$Model<-i
    
    if(!is.null(submodelNames)){     #useful if aesthetics according to a submodel group desired later 
      graph$submodel<-submodelNames[[i]]
    }
    
    
    graphlist[[i]]<-graph
    
    
  }
  
  graph<-bind_rows(graphlist)
  
  graph$starloc<-NA
  
  min<-min(graph$conf.low,na.rm=T)
  max<-max(graph$conf.high,na.rm=T)
  
  if(sig==TRUE){
    graph$starloc <- min-(max-min)/10
  }
  
  
  graph$Model<-as.factor(graph$Model)
  
  if(!is.null(ModelNames)){
    levels(graph$Model)<-ModelNames
  }
  
  
  ggplot(as.data.frame(graph),aes(x=term,y=estimate,colour=Model, shape=Model))+
    geom_point(position=position_dodge(w=0.5), size=2)+
    geom_errorbar(position=position_dodge(w=0.5),aes(ymin=conf.low,ymax=conf.high),size=0.3,width=0.2)+
    geom_hline(aes(yintercept=0),lty=2)+THEMEGRID+labs(x=NULL)+#coord_flip()+
    geom_text(aes(label=Sig,y=starloc),size=5, position=position_dodge(w=0.5), 
              show.legend = FALSE)+
    theme(axis.text.x = element_text(angle=90, vjust=0.6), strip.text.x = element_text(size = 12), 
          strip.background =element_rect(fill="white"))
  
}
