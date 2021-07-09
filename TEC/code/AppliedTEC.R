#Part 2 Applied TEC Framework
library(vegan)
library(tidyverse)
library(MuMIn)
library(ape)
library(adespatial)
library(betapart)
library(cowplot)
library(ggbiplot)

setwd("~/Dropbox/Users/matthewdouglasgreen/Dropbox/Manuscipts/L-S Biodviersity Streams_RCC_SDH")



species<-read.csv(file= "Analysis/TEC/sp.density.update.12.28.19_2.csv", row.name=1)
#species <- mutate_all(species, function(x) as.numeric(as.character(x)))
env <-read.csv(file= "Analysis/TEC/dave.matt.env.full.12.29.19.csv", row.name=1)
spa<-read.csv(file= "Analysis/TEC/spa.csv", row.name=1)

env<-env%>%dplyr::select(-c(WOODY_WETLANDS))
summary(env)

species<-species%>%dplyr::select(-c(Arachnida,Chironomidae,Nematomorpha,Nemerata,Oligochaeta,Ostracoda,Turbellaria,Euhirudinea))
#env<-env%>%mutate(Euc.dist.lake=log(1+Euc.dist.lake),River.dist.lake=log(1+River.dist.lake),Elevation=log(1+Elevation),Head.river.dist=log(1+Head.river.dist),Size.net.dist=Head.river.dist*Up.Lake.area,Size.river.dist=River.dist.lake*Up.Lake.area,Elev.dist=River.dist.lake/Elevation)

#Calcualte diversity and bind with envrioemtnal data, remove network if necessary

#aa%>%column_to_rownames("site")
Sites<-as.data.frame(rownames(all))

diversity<-species%>%
  #group_by(Site,Network)%>%
  transmute(N0=rowSums(species > 0),H= diversity(species),N1 =exp(H),N2 =diversity(species, "inv"),J= H/log(N0),E10= (N1/N0),E20= (N2/N0),
            Com.Size=rowSums(species),betas.LCBD=beta.div(species, method="hellinger",sqrt.D=TRUE)$LCBD ,betas.LCBD.p=beta.div(species, method="chord",sqrt.D=TRUE)$p.LCBD)

rownames(diversity)<-rownames(species)
diversity<-diversity%>%rownames_to_column("site")
env<-env%>%rownames_to_column("site")
species<-species%>%rownames_to_column("site")
#aa<-aa%>%rownames_to_column("site")

all<-left_join(species,env, by="site")
all<-left_join(all,diversity, by="site")
#all<-left_join(all,aa, by="site")

all<-all%>%filter(Network != "YOUNG")%>%column_to_rownames("site")

#Subset by network
Casc<-all%>%filter(O.NET=="CASCADE")
Evo<-all%>%filter(O.NET=="EVO")
Bubb<-all%>%filter(O.NET=="BUBBS")
Bubb<-Bubb%>%rownames_to_column("Site")%>%
  filter(Site !=	"Outlet.11007.fishless.2003")%>%
  filter(Site !=	"Outlet.11007.fishless.2004")%>%
  filter(Site !=	"Outlet.10487.trt.2003")%>%
  filter(Site !=	"Outlet.10487.trt.2004")%>%
  filter(Site !=	"Outlet.10477.trt.2003")%>%
  filter(Site !=	"Outlet.10477.trt.2004")%>%
  filter(Site !=	"Outlet.10494.trt.2012")%>%
  filter(Site !=	"	Outlet.Vidette.below.2003")%>%
  filter(Site !=	"	Outlet.Vidette.below.2004")%>%
  filter(Site !=	"	Outlet.Vidette.below.20012")%>%
  column_to_rownames("Site")


#Cascade
spatials<-Casc%>%dplyr::select(c(Head.river.dist, River.dist.lake, Up.Lake.area, Elevation))
envs<-Casc%>%dplyr::select(c(Temp,Chlorophyll.mean,Conductivity,DO,pH,Discharge.Mean,SHRUB_SCRUB,EVERGREEN_FOREST))

Epca = prcomp(envs, scale.=TRUE)
## and plot it
biplot(Epca)
##Examine and extract the axes
summary(Epca)
Epca$rotation
ggbiplot(Epca)
ggbiplot(Epca, labels=rownames(all$O.NET), groups=interaction(all$O.NET), ellipse=TRUE)

Epca$x
env_pc_scores <- data.frame(Epca$x[,1:6])
colnames(env_pc_scores)<-c("E_PC1","E_PC2","E_PC3","E_PC4","E_PC5","E_PC6")
rownames(env_pc_scores)<-rownames(Casc)#$site

#Space
Epca = prcomp(spatials, scale.=TRUE)
## and plot it
biplot(Epca)
##Examine and extract the axes
summary(Epca)
Epca$rotation


Epca$x
spatial_pc_scores <- data.frame(Epca$x[,1:4])
colnames(spatial_pc_scores)<-c("S_PC1","S_PC2","S_PC3","S_PC4")
rownames(spatial_pc_scores)<-rownames(Casc)

################
spatial_pc_scores<-spatial_pc_scores%>%rownames_to_column("Site")
env_pc_scores<-env_pc_scores%>%rownames_to_column("Site")
Casc<-Casc%>%rownames_to_column("Site")

all_dat<-left_join(Casc,spatial_pc_scores, by="Site")
all_dat<-left_join(all_dat,env_pc_scores, by="Site")
casc_all_dat<-all_dat%>%add_column(Com.Size.Gradient=log(Casc$Com.Size))

casc_all_dat

#Patterns of diversity acrosds environmental, spatial and communtiy size gradients

################################################################################################################
#EVO
spatials<-Evo%>%dplyr::select(c(Head.river.dist, River.dist.lake, Up.Lake.area, Elevation))
envs<-Evo%>%dplyr::select(c(Temp,Chlorophyll.mean,Conductivity,DO,pH,Discharge.Mean,SHRUB_SCRUB))

Epca = prcomp(envs, scale.=TRUE)
## and plot it
biplot(Epca)
##Examine and extract the axes
summary(Epca)
Epca$rotation
ggbiplot(Epca)
ggbiplot(Epca, labels=rownames(all$O.NET), groups=interaction(all$O.NET), ellipse=TRUE)

Epca$x
env_pc_scores <- data.frame(Epca$x[,1:6])
colnames(env_pc_scores)<-c("E_PC1","E_PC2","E_PC3","E_PC4","E_PC5","E_PC6")
rownames(env_pc_scores)<-rownames(Evo)#$site

#Space
Epca = prcomp(spatials, scale.=TRUE)
## and plot it
biplot(Epca)
##Examine and extract the axes
summary(Epca)
Epca$rotation


Epca$x
spatial_pc_scores <- data.frame(Epca$x[,1:4])
colnames(spatial_pc_scores)<-c("S_PC1","S_PC2","S_PC3","S_PC4")
rownames(spatial_pc_scores)<-rownames(Evo)

################
spatial_pc_scores<-spatial_pc_scores%>%rownames_to_column("Site")
env_pc_scores<-env_pc_scores%>%rownames_to_column("Site")
Evo<-Evo%>%rownames_to_column("Site")

all_dat<-left_join(Evo,spatial_pc_scores, by="Site")
all_dat<-left_join(all_dat,env_pc_scores, by="Site")
evo_all_dat<-all_dat%>%add_column(Com.Size.Gradient=log(Evo$Com.Size))

evo_all_dat


################################################################################################################
#Bubb
spatials<-Bubb%>%dplyr::select(c(Head.river.dist, River.dist.lake, Up.Lake.area, Elevation))
envs<-Bubb%>%dplyr::select(c(Temp,Chlorophyll.mean,Conductivity,DO,pH,Discharge.Mean,SHRUB_SCRUB,EVERGREEN_FOREST))

Epca = prcomp(envs, scale.=TRUE)
## and plot it
biplot(Epca)
##Examine and extract the axes
summary(Epca)
Epca$rotation
ggbiplot(Epca)
ggbiplot(Epca, labels=rownames(Bubb$Site) ,ellipse=F)

Epca$x
env_pc_scores <- data.frame(Epca$x[,1:6])
colnames(env_pc_scores)<-c("E_PC1","E_PC2","E_PC3","E_PC4","E_PC5","E_PC6")
rownames(env_pc_scores)<-rownames(Bubb)#$site

#Space
Epca = prcomp(spatials, scale.=TRUE)
## and plot it
biplot(Epca)
##Examine and extract the axes
summary(Epca)
Epca$rotation


Epca$x
spatial_pc_scores <- data.frame(Epca$x[,1:4])
colnames(spatial_pc_scores)<-c("S_PC1","S_PC2","S_PC3","S_PC4")
rownames(spatial_pc_scores)<-rownames(Bubb)

################
spatial_pc_scores<-spatial_pc_scores%>%rownames_to_column("Site")
env_pc_scores<-env_pc_scores%>%rownames_to_column("Site")
Bubb<-Bubb%>%rownames_to_column("Site")

all_dat<-left_join(Bubb,spatial_pc_scores, by="Site")
all_dat<-left_join(all_dat,env_pc_scores, by="Site")
bubb_all_dat<-all_dat%>%add_column(Com.Size.Gradient=log(Bubb$Com.Size))

bubb_all_dat


###################################################################################################################################################################################################
all_big_dat<-rbind(bubb_all_dat,evo_all_dat,casc_all_dat)


a<-all_big_dat%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = betas.LCBD, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  facet_grid(O.NET~ var, scales="free_x") +
  theme_bw()+
  ylab("Beta-diversity (LCBD)")+
  xlab("Log Community Size           Env. Gradeint (E_PC1)        Spa. Gradient (S_PC1)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

b<-all_big_dat%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = E10, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  facet_grid(O.NET~ var, scales = "free") +
  theme_bw()+
  ylab("Species Evenness")+
  xlab("Log Community Size           Env. Gradeint (E_PC1)        Spa. Gradient (S_PC1)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

c<-all_big_dat%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = N0, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  facet_grid(O.NET~ var, scales="free_x") +
  ylab("Species Richness")+
  xlab("Log Community Size           Env. Gradeint (E_PC1)        Spa. Gradient (S_PC1)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

plot_grid(a,b,c, ncol=3)


all_big_dat%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = River.dist.lake, y = value, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  facet_grid(O.NET~ var, scales="free_x") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")


new_labels <- c( "Com.Size.Gradient" = "Community Size Gradient","S_PC1" = "Spatial Gradient", "E_PC1" = "Environmental Gradient")

e2<-all_big_dat%>%
  pivot_longer(c(Com.Size.Gradient,S_PC1,E_PC1) , names_to = "key", values_to = "value")
summary(e2)

g1<-e2%>% 
  ggplot(aes(x = value, y =betas.LCBD, colour=O.NET )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="S_PC1"), shape=19)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="S_PC1"), shape=19)+
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="E_PC1"), shape=1)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="E_PC1"), shape=19)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="E_PC1"), shape=19)+
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="Com.Size.Gradient"), shape=1)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="Com.Size.Gradient"), shape=19)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="Com.Size.Gradient"), shape=1)+
  geom_smooth(data = filter(e2, O.NET =="BUBBS" & key=="S_PC1"),method="lm")+
  #geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="S_PC1"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="EVO"& key=="S_PC1"), method = "lm")+
 # geom_smooth(data=filter(e2,O.NET=="BUBBS"& key=="E_PC1"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="EVO"& key=="E_PC1"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="E_PC1"), method = "lm")+
  #geom_smooth(data = filter(e2, O.NET =="BUBBS" & key=="Com.Size.Gradient"),method="lm")+
  geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="Com.Size.Gradient"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="EVO"& key=="Com.Size.Gradient"), method = "lm")+
  facet_grid(O.NET~key,scales="free_x", labeller =labeller(key= new_labels)) + #axis.title.x+
  ylab("Local Contribution to Beta Diversity (LCBD)")+
  xlab("   Log Com. Size         Env. Gradeint (E_PC1)     Spa. Gradient (S_PC1)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

g2<-e2%>% 
  ggplot(aes(x = value, y =N0, colour=O.NET )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="E_PC1"), shape=19)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="E_PC1"), shape=19)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="E_PC1"), shape=1)+
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="Com.Size.Gradient"), shape=1)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="Com.Size.Gradient"), shape=19)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="Com.Size.Gradient"), shape=1)+
  #geom_smooth(data = filter(e2, O.NET =="BUBBS" & key=="S_PC1"),method="lm")+
  #geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="S_PC1"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="EVO"& key=="S_PC1"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="BUBBS"& key=="E_PC1"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="EVO"& key=="E_PC1"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="E_PC1"), method = "lm")+
  #geom_smooth(data = filter(e2, O.NET =="BUBBS" & key=="Com.Size.Gradient"),method="lm")+
  geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="Com.Size.Gradient"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="EVO"& key=="Com.Size.Gradient"), method = "lm")+
  facet_grid(O.NET~key,scales="free_x", labeller =labeller(key= new_labels)) + #axis.title.x+
  ylab("Species Richness")+
  xlab("   Log Com. Size         Env. Gradeint (E_PC1)     Spa. Gradient (S_PC1)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")


g3<-e2%>% 
  ggplot(aes(x = value, y =E10, colour=O.NET )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="S_PC1"), shape=19)+
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="E_PC1"), shape=1)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="E_PC1"), shape=1)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="E_PC1"), shape=19)+
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="Com.Size.Gradient"), shape=19)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="Com.Size.Gradient"), shape=19)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="Com.Size.Gradient"), shape=19)+
  #geom_smooth(data = filter(e2, O.NET =="BUBBS" & key=="S_PC1"),method="lm")+
  #geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="S_PC1"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="EVO"& key=="S_PC1"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="BUBBS"& key=="E_PC1"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="EVO"& key=="E_PC1"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="E_PC1"), method = "lm")+
  geom_smooth(data = filter(e2, O.NET =="BUBBS" & key=="Com.Size.Gradient"),method="lm")+
  geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="Com.Size.Gradient"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="EVO"& key=="Com.Size.Gradient"), method = "lm")+
  facet_grid(O.NET~key,scales="free_x", labeller =labeller(key= new_labels)) + #axis.title.x+
  ylab("Species Evenness")+
  xlab("   Log Com. Size         Env. Gradeint (E_PC1)     Spa. Gradient (S_PC1)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

plot_grid(g1,g2,g3, ncol=3)


#Option 1
#GLMS
cascade<-glm(betas.LCBD~S_PC1+E_PC1+Com.Size.Gradient,data=casc_all_dat)
cascade1<-glm(betas.LCBD~E_PC1+Com.Size.Gradient,data=casc_all_dat)
cascade2<-glm(betas.LCBD~S_PC1+Com.Size.Gradient,data=casc_all_dat)
cascade3<-glm(betas.LCBD~S_PC1+E_PC1,data=casc_all_dat)
cascade4<-glm(betas.LCBD~S_PC1,data=casc_all_dat)
cascade5<-glm(betas.LCBD~E_PC1,data=casc_all_dat)
cascade6<-glm(betas.LCBD~Com.Size.Gradient,data=casc_all_dat)
cascade7<-glm(betas.LCBD~1,data=casc_all_dat)
reported.table0 <- bbmle::AICtab(cascade,cascade1,cascade2,cascade3,cascade4,cascade5,cascade6,cascade7, weights = TRUE, sort = TRUE)
reported.table0 <- bbmle::AICtab(cascade,cascade1,cascade2,cascade3,cascade4,cascade5,cascade6,cascade7, weights = TRUE, sort = TRUE)

cascade<-glm(N0~S_PC1+E_PC1+Com.Size.Gradient,data=casc_all_dat)
cascade1<-glm(N0~E_PC1+Com.Size.Gradient,data=casc_all_dat)
cascade2<-glm(N0~S_PC1+Com.Size.Gradient,data=casc_all_dat)
cascade3<-glm(N0~S_PC1+E_PC1,data=casc_all_dat)
cascade4<-glm(N0~S_PC1,data=casc_all_dat)
cascade5<-glm(N0~E_PC1,data=casc_all_dat)
cascade6<-glm(N0~Com.Size.Gradient,data=casc_all_dat)
cascade7<-glm(N0~1,data=casc_all_dat)
reported.table0 <- bbmle::AICtab(cascade,cascade1,cascade2,cascade3,cascade4,cascade5,cascade6,cascade7, weights = TRUE, sort = TRUE)

evo<-glm(betas.LCBD~S_PC1+E_PC1+Com.Size.Gradient,data=evo_all_dat)
evo1<-glm(betas.LCBD~E_PC1+Com.Size.Gradient,data=evo_all_dat)
evo2<-glm(betas.LCBD~S_PC1+Com.Size.Gradient,data=evo_all_dat)
evo3<-glm(betas.LCBD~S_PC1+E_PC1,data=evo_all_dat)
evo4<-glm(betas.LCBD~S_PC1,data=evo_all_dat)
evo5<-glm(betas.LCBD~E_PC1,data=evo_all_dat)
evo6<-glm(betas.LCBD~Com.Size.Gradient,data=evo_all_dat)
evo7<-glm(betas.LCBD~1,data=evo_all_dat)
reported.table0 <- bbmle::AICtab(evo,evo1,evo2,evo3,evo4,evo5,evo6,evo7, weights = TRUE, sort = TRUE)
reported.table0 <- bbmle::AICtab(evo4,evo5,evo6,evo7, weights = TRUE, sort = TRUE)

bubb<-glm(betas.LCBD~S_PC1+E_PC1+Com.Size.Gradient,data=bubb_all_dat)
bubb1<-glm(betas.LCBD~E_PC1+Com.Size.Gradient,data=bubb_all_dat)
bubb2<-glm(betas.LCBD~S_PC1+Com.Size.Gradient,data=bubb_all_dat)
bubb3<-glm(betas.LCBD~S_PC1+E_PC1,data=bubb_all_dat)
bubb4<-glm(betas.LCBD~S_PC1,data=bubb_all_dat)
bubb5<-glm(betas.LCBD~E_PC1,data=bubb_all_dat)
bubb6<-glm(betas.LCBD~Com.Size.Gradient,data=bubb_all_dat)
bubb7<-glm(betas.LCBD~1,data=bubb_all_dat)
reported.table0 <- bbmle::AICtab(bubb,bubb1,bubb2,bubb3,bubb4,bubb5,bubb6,bubb7, weights = TRUE, sort = TRUE)
reported.table0 <- bbmle::AICtab(bubb4,bubb5,bubb6,bubb7, weights = TRUE, sort = TRUE)

casc4<-lm(betas.LCBD~S_PC1,data=casc_all_dat)
summary(casc4)
casc5<-lm(betas.LCBD~E_PC1,data=casc_all_dat)
summary(casc5)
casc6<-lm(betas.LCBD~Com.Size.Gradient,data=casc_all_dat)
summary(casc6)
casc7<-lm(S_PC1~E_PC1,data=casc_all_dat)
summary(casc7)

evo4<-lm(betas.LCBD~S_PC1,data=evo_all_dat)
summary(evo4)
evo5<-lm(betas.LCBD~E_PC1,data=evo_all_dat)
summary(evo5)
evo6<-lm(betas.LCBD~Com.Size.Gradient,data=evo_all_dat)
summary(evo6)
evo7<-lm(S_PC1~E_PC1,data=evo_all_dat)
summary(evo5)

bubb4<-lm(betas.LCBD~S_PC1,data=bubb_all_dat)
summary(bubb4)
bubb5<-lm(betas.LCBD~E_PC1,data=bubb_all_dat)
summary(bubb5)
bubb6<-lm(betas.LCBD~Com.Size.Gradient,data=bubb_all_dat)
summary(bubb6)
bubb7<-lm(S_PC1~E_PC1,data=bubb_all_dat)
summary(bubb7)

evo<-glm(N0~S_PC1+E_PC1+Com.Size.Gradient,data=evo_all_dat)
evo1<-glm(N0~E_PC1+Com.Size.Gradient,data=evo_all_dat)
evo2<-glm(N0~S_PC1+Com.Size.Gradient,data=evo_all_dat)
evo3<-glm(N0~S_PC1+E_PC1,data=evo_all_dat)
evo4<-glm(N0~S_PC1,data=evo_all_dat)
evo5<-glm(N0~E_PC1,data=evo_all_dat)
evo6<-glm(N0~Com.Size.Gradient,data=evo_all_dat)
evo7<-glm(N0~1,data=evo_all_dat)
reported.table0 <- bbmle::AICtab(evo,evo1,evo2,evo3,evo4,evo5,evo6,evo7, weights = TRUE, sort = TRUE)

casc4<-lm(N0~S_PC1,data=casc_all_dat)
summary(casc4)
casc5<-lm(N0~E_PC1,data=casc_all_dat)
summary(casc5)
casc6<-lm(N0~Com.Size.Gradient,data=casc_all_dat)
summary(casc6)
casc7<-lm(S_PC1~E_PC1,data=casc_all_dat)
summary(casc7)

evo4<-lm(N0~S_PC1,data=evo_all_dat)
summary(evo4)
evo5<-lm(N0~E_PC1,data=evo_all_dat)
summary(evo5)
evo6<-lm(N0~Com.Size.Gradient,data=evo_all_dat)
summary(evo6)
evo7<-lm(S_PC1~E_PC1,data=evo_all_dat)
summary(evo5)

bubb4<-lm(N0~S_PC1,data=bubb_all_dat)
summary(bubb4)
bubb5<-lm(N0~E_PC1,data=bubb_all_dat)
summary(bubb5)
bubb6<-lm(N0~Com.Size.Gradient,data=bubb_all_dat)
summary(bubb6)
bubb7<-lm(S_PC1~E_PC1,data=bubb_all_dat)
summary(bubb7)

casc4<-lm(E10~S_PC1,data=casc_all_dat)
summary(casc4)
casc5<-lm(E10~E_PC1,data=casc_all_dat)
summary(casc5)
casc6<-lm(E10~Com.Size.Gradient,data=casc_all_dat)
summary(casc6)
casc7<-lm(S_PC1~E_PC1,data=casc_all_dat)
summary(casc7)

evo4<-lm(E10~S_PC1,data=evo_all_dat)
summary(evo4)
evo5<-lm(E10~E_PC1,data=evo_all_dat)
summary(evo5)
evo6<-lm(E10~Com.Size.Gradient,data=evo_all_dat)
summary(evo6)
evo7<-lm(S_PC1~E_PC1,data=evo_all_dat)
summary(evo5)

bubb4<-lm(E10~S_PC1,data=bubb_all_dat)
summary(bubb4)
bubb5<-lm(E10~E_PC1,data=bubb_all_dat)
summary(bubb5)
bubb6<-lm(E10~Com.Size.Gradient,data=bubb_all_dat)
summary(bubb6)
bubb7<-lm(S_PC1~E_PC1,data=bubb_all_dat)
summary(bubb7)

bubb_all_dat_2<-bubb_all_dat%>%
  filter(Site !=	"Outlet.11007.fishless.2003")%>%
  filter(Site !=	"Outlet.11007.fishless.2004")%>%
  filter(Site !=	"Outlet.10487.trt.2003")%>%
  filter(Site !=	"Outlet.10487.trt.2004")%>%
  filter(Site !=	"Outlet.10477.trt.2003")%>%
  filter(Site !=	"Outlet.10477.trt.2004")%>%
  filter(Site !=	"Outlet.10494.trt.2012")%>%
  filter(Site !=	"	Outlet.Vidette.below.2003")%>%
  filter(Site !=	"	Outlet.Vidette.below.2004")%>%
  filter(Site !=	"	Outlet.Vidette.below.20012")

bubb4<-lm(N0~S_PC1,data=bubb_all_dat_2)
summary(bubb4)
bubb5<-lm(N0~E_PC1,data=bubb_all_dat_2)
summary(bubb5)
bubb6<-lm(N0~Com.Size.Gradient,data=bubb_all_dat_2)
summary(bubb6)
bubb7<-lm(S_PC1~E_PC1,data=bubb_all_dat_2)
summary(bubb7) 



#Option 2
#SEM's
##### Lavaan modeling
library(lavaan)
library(semPlot)

smod1 = ' betas.LCBD~Space+Com.Size.Gradient+Env
          Space=~S_PC1+S_PC2+S_PC3+S_PC4
          Env=~E_PC1+E_PC2+E_PC3+E_PC4
         '

#space<-'pred.persistence =~connect.per + meta.size +total.vol '
smod1.fit <- sem(smod1,data=all_dat)
summary(smod1.fit,standardized=TRUE,rsq=T)
fitMeasures(smod1.fit)
modindices(smod1.fit)

#quick plot of path analysis
semPaths(smod1.fit, what='std', layout = "tree3", intercepts = FALSE, residuals = FALSE, 
         edge.label.cex=1.25, curvePivot = FALSE,  fade=FALSE, rotation = 2)

##########


a.casc<-bubb_all_dat_2%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = betas.LCBD, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  facet_grid(O.NET~ var, scales="free_x") +
  theme_bw()+
  ylab("Beta-diversity (LCBD)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

b.casc<-casc_all_dat%>%
  gather(S_PC2,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = E10, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  facet_grid(O.NET~ var, scales = "free") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

c.casc<-bubb_all_dat_2%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = N0, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  facet_grid(O.NET~ var, scales="free_x") +
  ylab("Species Richness")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

d.casc<-casc_all_dat%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = N1, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  facet_grid(O.NET~ var, scales="free_x") +
  ylab("Species Diversity")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

plot_grid(a,b,c,d, ncol=4)
plot_grid(a,c,d, ncol=3)


#Patterns of diversity acrosds environmental, spatial and communtiy size gradients


a.evo<-evo_all_dat%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = betas.LCBD, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  facet_grid(O.NET~ var, scales="free_x") +
  theme_bw()+
  ylab("Beta-diversity (LCBD)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

b.evo<-evo_all_dat%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = E10, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  facet_grid(O.NET~ var, scales = "free") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

c.evo<-evo_all_dat%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = N0, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  facet_grid(O.NET~ var, scales="free_x") +
  ylab("Species Richness")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

d.evo<-evo_all_dat%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = N1, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=T)+
  facet_grid(O.NET~ var, scales="free_x") +
  ylab("Species Diversity")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

plot_grid(a,b,c,d, ncol=4)
plot_grid(a,c,d, ncol=3)

#Patterns of diversity acrosds environmental, spatial and communtiy size gradients

a.bubb<-bubb_all_dat%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = betas.LCBD, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  facet_grid(O.NET~ var, scales="free_x") +
  theme_bw()+
  ylab("Beta-diversity (LCBD)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

b.bubb<-bubb_all_dat%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = E10, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  facet_grid(O.NET~ var, scales = "free") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

c.bubb<-bubb_all_dat%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = N0, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  facet_grid(O.NET~ var, scales="free_x") +
  ylab("Species Richness")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

d.bubb<-bubb_all_dat%>%
  gather(S_PC1,Com.Size.Gradient,E_PC1,key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = N1, colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there tredns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  facet_grid(O.NET~ var, scales="free_x") +
  ylab("Species Diversity")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

plot_grid(a,b,c,d, ncol=4)
plot_grid(a,c,d, ncol=3)
