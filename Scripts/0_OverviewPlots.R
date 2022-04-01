library(tidyverse); library(ggregplot); library(RColorBrewer); library(nationalparkcolors)

ParaColOverview <- park_palette("Arches", n=5)

ParaModDF <- read_csv("Data/ParaModDF.csv") %>% 
  mutate(Season=as.factor(Season),
         Season=fct_relevel(Season, "Winter19", "Spring19", "Summer19","Autumn19", "Winter20"))

ParaBin <- ParaModDF %>% select(ends_with("Bin")) %>% colnames()

ParaModLong <- ParaModDF %>% 
  select(ParaBin, "Season", "AgeClass", "Id", "Sex") %>% 
  gather(key="Parasite", value="InfStatus", 1:7) %>% 
  mutate(Id=as.factor(Id)) %>% 
  mutate(Parasite=str_remove(Parasite, "Bin")) %>% 
  #filter(!(Season=="Winter"&Parasite=="Nematodirus")) %>% 
  #filter(!(Season=="Spring"&Parasite=="Nematodirus")) %>% 
  mutate(IdParasite = paste(Id, Parasite, sep=":"), 
         Unique = paste(Id, Parasite,Season, sep=":")) 

IdAge <- ParaModLong %>% select(Id, AgeClass) %>% 
  filter(!duplicated(Id))

ParaModLong2 <- ParaModLong %>% 
  group_by(Id, Season, Parasite) %>% 
  summarise(InfStatus=round(mean(InfStatus))) %>% as.data.frame() %>% 
  #filter(InfStatus==max(InfStatus)) %>% 
  #filter(!duplicated(Unique))
  mutate(IdParasite = paste(Id, Parasite, sep=":"), 
         Unique = paste(Id, Parasite,Season, sep=":")) 

ParaModWide <- ParaModLong2 %>% 
  pivot_wider(
              id_cols = "IdParasite", 
              values_from = "InfStatus", 
              names_from = "Season")

ParaModLong3 <- pivot_longer(ParaModWide, 
                            cols = levels(ParaModDF$Season),
                            values_to = "InfStatus", 
                            names_to = "Season") %>% as.data.frame() %>% 
  mutate(Parasite=str_split(IdParasite, "[[:]]") %>% map_chr(.,2),
         Id=str_split(IdParasite, "[[:]]") %>% map_chr(.,1),
         Unique = paste(Id, Parasite, Season, sep=":")) %>% 
  left_join(., IdAge, by="Id", all.x=TRUE)


ParaPlotLong <- 
  ParaModLong3 %>%  
  mutate(SeasonId = paste(Id, Season, sep="_"), 
         ParasiteId = paste(Id, Parasite, sep= "_")) %>% 
  select(SeasonId, ParasiteId, Season, AgeClass, Id, Parasite, InfStatus) 


Labels <- ParaPlotLong %>% 
      pull(SeasonId) %>% unique %>% sort %>% as.character

FillLabels <-
  lapply(1:length(Labels), function(a) { 
    rep(" ", length(Labels[[a]])*1.5) -> FillLabels
    return(FillLabels)
  })   

for(i in 1:length(FillLabels)){
  FillLabels[[i]][1:length(Labels[[i]]) + 
                    rep(c(0:length(Labels[[i]])), each = 2)[1:length(Labels[[i]])]] <-
    Labels[[i]]
} 


ParaPlotLong %>% 
  mutate(Parasite = fct_relevel(Parasite, "Coccidia", "Strongyles", "Capillaria", 
                                "Nematodirus", "Strongyloides", "Moniezia", "Trichuris"),
         Season= fct_relevel(Season, "Winter19", "Spring19", "Summer19", "Autumn19", "Winter20")) %>% 
  mutate(InfStatus = as.factor(as.character(InfStatus)), 
         InfStatusAge = case_when(InfStatus == "0" ~ "0", 
                                  InfStatus == "1" & AgeClass == "Lamb" ~ "1 (Lamb)", 
                                  InfStatus == "1" & AgeClass == "Adult" ~ "1 (Adult)")) %>% 
  #mutate(InfStatusAge = fct_relevel(InfStatusAge, "0", "1 (Lamb)"), 
  #       AgeClass = fct_relevel(AgeClass, "Lamb")) %>% 
  ggplot(aes(x=Season, y=Id, fill=InfStatusAge)) + 
  geom_tile(colour="grey") + 
  scale_fill_manual(values=c("white", ParaColOverview[4], ParaColOverview[5]), na.value="grey89") + 
  facet_grid(AgeClass ~ Parasite, scales="free_y") + 
  theme_bw(base_size = 18) + 
  labs(y='Sheep ID') +
  theme(axis.text.y = element_blank(), 
        axis.text.x = element_text(angle=45, hjust=1, colour="black"), 
        strip.background = element_rect(fill="white")) + 
  #scale_x_discrete(limits=c("Coccidea", "Strongyles", "Capillaria", 
   #                         "Nematodirus", "Strongyloides", "Monezia", "Trichuris")) + 
  ggsave('Figures/InfStatusHeatMapFin.png', units='mm', height=250, width=350, dpi=500)


# supp tables  -------------------------------------------------------------

Para <- c("Strongyles", "Coccidia", "Capillaria", "Nematodirus")

ParaLongFEC <- ParaModDF %>% select(Para, covar, "Id", "Sex") %>% 
  gather(key="Parasite", value="FEC", 1:4) %>% 
  mutate(Id=as.factor(Id)) %>% 
  mutate(Parasite=str_remove(Parasite, "Bin")) %>% 
  filter(!(Season=="Winter"&Parasite=="Nematodirus")) %>% 
  filter(!(Season=="Spring"&Parasite=="Nematodirus"))

ParaLongFEC %>% 
  #mutate(Parasite = fct_relevel(Parasite, "Coccidia", "Strongyles", "Capillaria", 
   #                             "Nematodirus", "Strongyloides", "Moniezia", "Trichuris")) %>% 
  group_by(Parasite, Season, AgeClass) %>% 
  summarise(Mean = mean(FEC), 
            Min = min(FEC), 
            Max = max(FEC)) %>% as.data.frame() %>% 
  mutate_if(is.numeric, round, digits=2) %>%  
  mutate(Range= paste(Min, Max, sep=" - ")) %>% 
  arrange(Parasite, AgeClass, match(Season, c("Winter19", "Spring19", "Summer19", "Autumn19", "Winter20")))-> suppTable 

DemoNos <- ParaModDF %>% 
  mutate(Season=fct_relevel(Season, "Winter19", "Spring19", "Summer19", "Autumn19")) %>% 
  group_by(AgeClass, Season, SexRepro) %>% 
  count() 

write.csv(suppTable, "Output/rawMeansTable.csv", row.names = F)
write.csv(DemoNos, "Output/DemoSampleNos.csv", row.names = F)
