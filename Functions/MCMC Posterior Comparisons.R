# Liverpool data posterior comparisons 

MCMCPostCompInt<-function(model, df, vars){
  require(tidybayes)
  InteractionComps<-list()  
  
  vars<-vars[c(1, length(vars))]
  enddf<-list()
  if(!is.null(vars)){
    for(x in 1:length(vars)){
      #Factors <- vars[x]
      #Columns<-Solutions[,substr(colnames(Solutions),1,nchar(Factors))==Factors] %>% colnames()
      selected<-unlist(str_split(vars[x], "[:]"))[2]
      
      df1<-
        model %>% 
        recover_types(df) %>% 
        tidy_draws() %>% 
        select(contains(selected)) %>% 
        select(-selected) %>%  as.data.frame()
      
      Names <-colnames(df1) %>% as.data.frame()
      NamesNew <- Names %>% separate(1, into=c("Level", "Term"), sep="[:]") %>% pull(Level)
      colnames(df1)<-NamesNew
      LevelCombs <- combn(ncol(df1),2) %>% t()
      
      df2 = apply(LevelCombs, 1, function(a){
        
        df3 = data.frame(Value = df1[,a[2]] - df1[,a[1]] ) %>% 
          mutate(Comp = paste(colnames(df1)[a[2]],colnames(df1)[a[1]], sep = "-")) 
        
      }) %>% bind_rows() 
      
      df4 = 1:ncol(df1) %>% 
        lapply(function(a) data.frame(Value = c(df1[,a]), Comp = paste(colnames(df1)[a],"Intercept", sep = "-"))) %>% 
        bind_rows()
      

      enddf[[x]] = bind_rows(df2, df4) %>% 
        mutate(Term=selected, Comp=as.factor(as.character(Comp)))
    }}
  
  enddf
}

MCMCPostComp<-function(model, df, vars){
  require(tidybayes)
  InteractionComps<-list()  
  
  enddf<-list()
  if(!is.null(vars)){
    for(x in 1:length(vars)){
      #Factors <- vars[x]
      #Columns<-Solutions[,substr(colnames(Solutions),1,nchar(Factors))==Factors] %>% colnames()
      selected<-unlist(str_split(vars[x], "[:]"))[2]
      
      df1<-
        model %>% 
        recover_types(df) %>% 
        tidy_draws() %>% 
        select(contains(vars)) %>% 
        select(-contains(":")) %>%  as.data.frame()
      
      Names <-colnames(df1) 
      LevelCombs <- combn(ncol(df1),2) %>% t()
      
      df2 = apply(LevelCombs, 1, function(a){
        
        df3 = data.frame(Value = df1[,a[2]] - df1[,a[1]] ) %>% 
          mutate(Comp = paste(colnames(df1)[a[2]],colnames(df1)[a[1]], sep = "-")) 
        
      }) %>% bind_rows() 
      
      df4 = 1:ncol(df1) %>% 
        lapply(function(a) data.frame(Value = c(df1[,a]), Comp = paste(colnames(df1)[a],"Intercept", sep = "-"))) %>% 
        bind_rows()
      
      
      enddf[[x]] = bind_rows(df2, df4) %>% 
        mutate(Comp=as.factor(as.character(Comp)))
    }}
  
  enddf
}
