#paper title: Rethinking Biodiversity in Stream Ecology Frameworks
#Part 2 Applied TEC Framework
#Author: Matthew Douglas Green
#Date: Sept 14,2021

########################################################################################################################
#load packages
Packages <- c("olsrr","semPlot","lavaan","MuMIn","lme4","vegan", "ggplot2", "tidyverse", "ape","MuMIn","adespatial", "betapart", "cowplot")
lapply(Packages, library, character.only = TRUE)
library(ggbiplot)
library(insight)
library(performance)
library(glmmTMB)

########################################################################################################################
#Open Data and Clean
setwd("~/Dropbox/Users/matthewdouglasgreen/Dropbox/Manuscipts/L-S Biodviersity Streams_RCC_SDH")
getwd()
species<-read.csv(file= "TEC/data/sp.density.update.12.28.19.csv", row.name=1)
env <-read.csv(file= "TEC/data/dave.matt.env.full.12.29.19.csv", row.name=1)

env<-env%>%dplyr::select(-c(WOODY_WETLANDS))
summary(env)

species<-species%>%dplyr::select(-c(Arachnida,Chironomidae,Nematomorpha,Oligochaeta,Ostracoda,Turbellaria,Euhirudinea))
#env<-env%>%mutate(Euc.dist.lake=log(1+Euc.dist.lake),River.dist.lake=log(1+River.dist.lake),Elevation=log(1+Elevation),Head.river.dist=log(1+Head.river.dist),Size.net.dist=Head.river.dist*Up.Lake.area,Size.river.dist=River.dist.lake*Up.Lake.area,Elev.dist=River.dist.lake/Elevation)

########################################################################################################################
#Calculate diversity and bind with environmental data, remove network if necessary

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

all<-left_join(species,env, by="site")
all<-left_join(all,diversity, by="site")
#all<-left_join(all,aa, by="site")
all<-all%>%column_to_rownames("site")

#all<-all%>%filter(Network != "YOUNG")%>%column_to_rownames("site")

#Subset by network
Casc<-all%>%filter(O.NET=="CASCADE")
Evo<-all%>%filter(O.NET=="EVO")
Bubb<-all%>%filter(O.NET=="BUBBS")
Rock<-all%>%filter(O.NET=="ROCK")

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

########################################################################################################################
#PCA's Gradients of Spatial and Environment

#1)Cascade
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
casc_all_dat<-casc_all_dat%>%
  mutate(Spatial=-1*S_PC1)
#Patterns of diversity acrosds environmental, spatial and communtiy size gradients

################################################################################################################
#2)EVO
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

evo_all_dat<-evo_all_dat%>%
  mutate(Spatial=-1*S_PC1)
################################################################################################################
#3) Bubb
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

bubb_all_dat<-bubb_all_dat%>%
  mutate(Spatial=S_PC1)
################################################################################################################
#4) Rock
spatials<-Rock%>%dplyr::select(c(Head.river.dist, River.dist.lake, Up.Lake.area, Elevation))
envs<-Rock%>%dplyr::select(c(Temp,Chlorophyll.mean,Conductivity,DO,pH,Discharge.Mean,SHRUB_SCRUB,EVERGREEN_FOREST))

Epca = prcomp(envs, scale.=TRUE)
## and plot it
biplot(Epca)
##Examine and extract the axes
summary(Epca)
Epca$rotation
ggbiplot(Epca)
ggbiplot(Epca, labels=rownames(Rock$Site) ,ellipse=F)

Epca$x
env_pc_scores <- data.frame(Epca$x[,1:6])
colnames(env_pc_scores)<-c("E_PC1","E_PC2","E_PC3","E_PC4","E_PC5","E_PC6")
rownames(env_pc_scores)<-rownames(Rock)#$site

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
rownames(spatial_pc_scores)<-rownames(Rock)

################
spatial_pc_scores<-spatial_pc_scores%>%rownames_to_column("Site")
env_pc_scores<-env_pc_scores%>%rownames_to_column("Site")
Rock<-Rock%>%rownames_to_column("Site")

all_dat<-left_join(Rock,spatial_pc_scores, by="Site")
all_dat<-left_join(all_dat,env_pc_scores, by="Site")

Rock_all_dat<-all_dat%>%add_column(Com.Size.Gradient=log(Rock$Com.Size))

Rock_all_dat<-Rock_all_dat%>%
  mutate(Spatial=-1*S_PC1)



###################################################################################################################################################################################################
#Combine PCA's with Diversity data and Species pool
all_big_dat<-rbind(bubb_all_dat,evo_all_dat,casc_all_dat,Rock_all_dat)
all_big_dat$Reg.Pool<-if_else(all_big_dat$O.NET =="BUBBS","88", 
                              if_else(all_big_dat$O.NET =="EVO","39",
                                      if_else(all_big_dat$O.NET =="ROCK","67",
                                              if_else(all_big_dat$O.NET =="KERN","56","47"))))

all_big_dat%>%
  ggplot(aes(x=Reg.Pool,y=N1, fill=O.NET))+
  geom_boxplot()

all_big_dat%>%
  ggplot(aes(x=Reg.Pool,y=betas.LCBD, fill=O.NET))+
  geom_boxplot()

####################################################################################################################################################################################
#Plotting
p1<-all_big_dat%>%
  #filter(Head.river.dist>3)%>%
  ggplot(aes(x = E_PC1, y =betas.LCBD )) + 
  ggtitle("d)") +#remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point()+
  stat_smooth(method = glm,method.args = list(family = gaussian(link = "identity")))+
  #geom_smooth(method = "lm")+
  #facet_wrap(~Network,scales="free") + 
  theme_bw()+
  xlab("Environmental Gradient") +labs(y=(("Local Contribution to \u03B2-diversity (LCBD) ")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

p2<-all_big_dat%>%
  #filter(Head.river.dist>3)%>%
  ggplot(aes(x = Spatial, y =betas.LCBD )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  ggtitle("b)") +
  geom_point()+
  stat_smooth(method = glm,method.args = list(family = gaussian(link = "identity")))+
  #geom_smooth(method = "lm")+
  xlab("Spatial Gradient") +labs(y=(("Local Contribution to \u03B2-diversity (LCBD) ")))+
  #facet_wrap(~Network,scales="free") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

p3<-all_big_dat%>%
  #filter(Head.river.dist>3)%>%
  ggplot(aes(x = log(Com.Size+1), y =betas.LCBD )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  ggtitle("f)") +
  geom_point()+
  stat_smooth(method = glm,method.args = list(family = gaussian(link = "identity")))+
  #geom_smooth(method = "lm")+
  xlab("Log Community Size") +labs(y=(("Local Contribution to \u03B2-diversity (LCBD) ")))+
  #facet_wrap(~Network,scales="free") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

p4<-all_big_dat%>%
  #filter(Head.river.dist>3)%>%
  ggplot(aes(x = E_PC1, y =N1 )) +
  ggtitle("d)") +#remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point()+
  #stat_smooth(method = glm,method.args = list(family = gaussian(link = "identity")))+
  #geom_smooth(method = "lm")+
  xlab(" Environmental Gradient") +labs(y=(("Shannon Diversity")))+
  #facet_wrap(~Network,scales="free") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

p5<-all_big_dat%>%
  # filter(Head.river.dist>3)%>%
  ggplot(aes(x = Spatial, y =N1 )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  ggtitle("b)") +
  geom_point()+
  stat_smooth(method = glm,method.args = list(family = gaussian(link = "identity")))+
  #geom_smooth(method = "lm")+
  #facet_wrap(~Network,scales="free") + 
  xlab("Spatial Gradient") +labs(y=(("Shannon Diversity")))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

p6<-all_big_dat%>%
  #filter(Head.river.dist>3)%>%
  ggplot(aes(x = log(Com.Size+1), y =N1 )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  ggtitle("f)") +
  geom_point()+
  stat_smooth(method = glm,method.args = list(family = gaussian(link = "identity")))+
  #geom_smooth(method = "lm")+
  #facet_wrap(~Network,scales="free") + 
  xlab(" Log Community Size") +labs(y=(("Shannon Diversity")))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

plot_grid(p2,p5,p1,p4,p3,p6,ncol=2)
############################################################################################################################################
#GLM's

#LCBD
mod1<-betareg(betas.LCBD~Spatial,data=all_big_dat)
mod2<-betareg(betas.LCBD~E_PC1,data=all_big_dat)
mod3<-betareg(betas.LCBD~Com.Size.Gradient,data=all_big_dat)
mod4<-betareg(betas.LCBD~Spatial*E_PC1,data=all_big_dat)
mod5<-betareg(betas.LCBD~Spatial*Com.Size.Gradient,data=all_big_dat)
mod6<-betareg(betas.LCBD~E_PC1*Com.Size.Gradient,data=all_big_dat)
mod7<-betareg(betas.LCBD~E_PC1*Com.Size.Gradient*Spatial,data=all_big_dat)
null<-betareg(betas.LCBD~1,data=all_big_dat)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod3,mod4,mod5,mod6, null,weights = TRUE, sort = FALSE)
r2(mod1)
r2(mod2)
r2(mod3)
r2(mod4)
r2(mod5)
r2(mod6)
r2(null)
anova(mod1,mod2,mod5,mod6, null)

pseudoR1 <- ((mod1$null.deviance-mod1$deviance)/mod1$null.deviance)
pseudoR2 <- ((mod2$null.deviance-mod2$deviance)/mod2$null.deviance)
pseudoR3 <- ((mod3$null.deviance-mod3$deviance)/mod3$null.deviance)
pseudoR4 <- ((mod4$null.deviance-mod4$deviance)/mod4$null.deviance)
pseudoR5 <- ((mod5$null.deviance-mod5$deviance)/mod5$null.deviance)
pseudoR6 <- ((mod6$null.deviance-mod6$deviance)/mod6$null.deviance)
pseudoRnull <- ((null$null.deviance-null$deviance)/null$null.deviance)

#N1
mod1<-glm(N1~Spatial,family=gaussian(link = "identity"),data=all_big_dat)
mod2<-glm(N1~E_PC1,family=gaussian(link = "identity"),data=all_big_dat)
mod3<-glm(N1~Com.Size.Gradient,family=gaussian(link = "identity"),data=all_big_dat)
mod4<-glm(N1~Spatial*E_PC1,family=gaussian(link = "identity"),data=all_big_dat)
mod5<-glm(N1~Spatial*Com.Size.Gradient,family=gaussian(link = "identity"),data=all_big_dat)
mod6<-glm(N1~E_PC1*Com.Size.Gradient,family=gaussian(link = "identity"),data=all_big_dat)
mod7<-glm(N1~E_PC1*Com.Size.Gradient*Spatial,family=gaussian(link = "identity"),data=all_big_dat)
null<-glm(N1~1,family=gaussian(link = "identity"),data=all_big_dat)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod3,mod4,mod5,mod6, null,weights = TRUE, sort = FALSE)
r2(mod1)
r2(mod2)
r2(mod3)
r2(mod4)
r2(mod5)
r2(mod6)
r2(null)
anova(mod1,mod2,mod5,mod6, mod7null)

pseudoR1 <- ((mod1$null.deviance-mod1$deviance)/mod1$null.deviance)
pseudoR2 <- ((mod2$null.deviance-mod2$deviance)/mod2$null.deviance)
pseudoR3 <- ((mod3$null.deviance-mod3$deviance)/mod3$null.deviance)
pseudoR4 <- ((mod4$null.deviance-mod4$deviance)/mod4$null.deviance)
pseudoR5 <- ((mod5$null.deviance-mod5$deviance)/mod5$null.deviance)
pseudoR6 <- ((mod6$null.deviance-mod6$deviance)/mod6$null.deviance)
pseudoRnull <- ((null$null.deviance-null$deviance)/null$null.deviance)

########################################################################################################################################################################
#Nww  Stuff

#Gaussian dist
qqnorm(all_big_dat$N1)
model <- lm(N1 ~ Spatial, data = all_big_dat)
ols_plot_resid_fit(model)

qqnorm(all_big_dat$betas.LCBD)
model <- lm(betas.LCBD ~ Spatial, data = all_big_dat)
ols_plot_resid_fit(model)

###Autocorrelation among variables
cor(all_big_dat$E_PC1,all_big_dat$Spatial, method = "pearson")
cor(all_big_dat$Com.Size.Gradient,all_big_dat$Spatial, method = "pearson")
cor(all_big_dat$E_PC1,all_big_dat$Com.Size.Gradient, method = "pearson")

cor(all_big_dat$Head.river.dist,all_big_dat$Com.Size.Gradient, method = "pearson")
cor(all_big_dat$Head.river.dist,all_big_dat$Spatial, method = "pearson")
cor(all_big_dat$Head.river.dist,all_big_dat$E_PC1, method = "pearson")
cor(all_big_dat$River.dist.lake,all_big_dat$Com.Size.Gradient, method = "pearson")
cor(all_big_dat$River.dist.lake,all_big_dat$Spatial, method = "pearson")
cor(all_big_dat$River.dist.lake,all_big_dat$E_PC1, method = "pearson")

#Supplemtary Figure 1:
  
s1<-all_big_dat%>%
  filter(log(River.dist.lake+1) >0)%>%
  ggplot(aes(x=log(Head.river.dist+1),y=E_PC1))+
  geom_point()+
  ggtitle("a)") +
  geom_smooth(method = "lm",se=F)+ xlab("Distance from Headwaters (m)")+ylab("Environmental Gradient")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.border = element_rect(colour = "black"))
s2<-all_big_dat%>%
  filter(log(River.dist.lake+1) >0)%>%
  ggplot(aes(x=log(Head.river.dist+1),y=Spatial))+
  geom_point()+
  ggtitle("b)") +
  geom_smooth(method = "lm",se=F)+ xlab("Distance from Headwaters (m)")+ylab("Spatial Gradient")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.border = element_rect(colour = "black"))

s3<-all_big_dat%>%
  filter(log(River.dist.lake+1) >0)%>%
  ggplot(aes(x=log(Head.river.dist+1),y=Com.Size.Gradient))+
  geom_point()+
  ggtitle("c)") +
  #geom_smooth(method = "lm",se=F)+ 
  xlab("Distance from Headwaters (m)")+ylab("Community Size Gradient")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.border = element_rect(colour = "black"))

s4<-all_big_dat%>%
  filter(log(River.dist.lake+1) >0)%>%
  ggplot(aes(x=log(River.dist.lake+1),y=E_PC1))+
  geom_point()+
  ggtitle("d)") +
  #geom_smooth(method = "lm",se=F)+
  xlab("Distance from Upstream Lakes (m)")+ylab("Environmental Gradient")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black"))
s5<-all_big_dat%>%
  filter(log(River.dist.lake+1) >0)%>%
  ggplot(aes(x=log(River.dist.lake+1),y=Spatial))+
  geom_point()+
  ggtitle("e)") +
  geom_smooth(method = "lm",se=F)+xlab("Distance from Upstream Lakes (m)")+ylab("Spatial Gradient")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black"))

s6<-all_big_dat%>%
  filter(log(River.dist.lake+1) >0)%>%
  ggplot(aes(x=log(River.dist.lake+1),y=Com.Size.Gradient))+
  ggtitle("f)") +
  geom_point()+
  geom_smooth(method = "lm",se=F)+xlab("Distance from Upstream Lakes (m)")+ylab("Community Size Gradient")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(colour = "black"))
plot_grid(s1,s2,s3,s4,s5,s6, nrow=2)

all_big_dat%>%
  filter(log(River.dist.lake+1) >0)%>%
  gather(E_PC1,Spatial,Com.Size.Gradient, key = "var", value = "value") %>%
  ggplot(aes(x=log(River.dist.lake+1),y=value))+
  geom_point()+
  geom_smooth(method = "lm",se=F)+xlab("Distance from Upstream Lakes (m)")+
  facet_wrap(~var,scales="free")+theme_bw()+theme(panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       panel.border = element_rect(colour = "black"))

all_big_dat%>%
  filter(log(River.dist.lake+1) >0)%>%
  gather(SHRUB_SCRUB,Chlorophyll.mean,Temp, Conductivity,  DO,   pH, Discharge.Mean, key = "var", value = "value") %>%
  ggplot(aes(x=log(Head.river.dist+1),y=value))+
  geom_point()+
  geom_smooth(method = "lm",se=F)+xlab("Distance from Headwaters (m)")+
  facet_wrap(~var,scales="free")+theme_bw()+theme(panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.border = element_rect(colour = "black"))

all_big_dat%>%
  filter(log(River.dist.lake+1) >0)%>%
  gather(SHRUB_SCRUB,Chlorophyll.mean,Temp, Conductivity,  DO,   pH, Discharge.Mean, key = "var", value = "value") %>%
  ggplot(aes(x=log(River.dist.lake+1),y=value))+
  geom_point()+
  geom_smooth(method = "lm",se=F)+xlab("Distance from Upstream Lakes (m)")+
  facet_wrap(~var,scales="free")+theme_bw()+theme(panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       panel.border = element_rect(colour = "black"))

reg.pool<-c(88,47,39,56,67)
beta<-c(0.6313422,0.5546429,0.548144,0.7042006,0.5172451)

df<-cbind(reg.pool,beta)
df<-as.data.frame(df)
df%>%
  ggplot(aes(x=(reg.pool),y=beta))+
  geom_point()+
  geom_smooth(method = "lm")

all_big_dat%>%
  ggplot(aes(x=E_PC1,y=Spatial))+
  geom_point(aes(x=E_PC1,y=Spatial, colour=N1, size=N1))+
  geom_smooth(method = "lm",se=F)+
  xlab("Env")+ theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  facet_grid(~O.NET)

all_big_dat%>%
  gather(E_PC1,Spatial, Com.Size.Gradient,key = "var", value = "value") %>% 
  ggplot(aes(x = log(Head.river.dist+1), y = value)) + geom_point() +geom_smooth(method = "lm")+
  facet_wrap(~ var, scales = "free") + theme_bw()

all_big_dat%>%
  gather(E_PC1,Spatial, Com.Size.Gradient,key = "var", value = "value") %>% 
  ggplot(aes(x = log(River.dist.lake+1), y = value)) + geom_point() +geom_smooth(method = "lm")+
  facet_wrap(~ var, scales = "free") + theme_bw()

########################################################################################################################################################################
#SEMS
summary(all_big_dat$Spatial)
summary(all_big_dat$E_PC1)
summary(all_big_dat$Com.Size.Gradient)
summary(all_big_dat$betas.LCBD)
summary(all_big_dat$N1)



all_big_datz<-all_big_dat%>%
  mutate(N1=N1*0.001+0.05,Spatial=Spatial*.01+0.05, E_PC1=E_PC1*.01+0.05,Com.Size.Gradient=Com.Size.Gradient*.01+0.05,betas.LCBD=betas.LCBD*10+0.05)

summary(all_big_datz$Spatial)
summary(all_big_datz$E_PC1)
summary(all_big_datz$Com.Size.Gradient)
summary(all_big_datz$betas.LCBD)
summary(all_big_datz$N1)



smod1 = ' N1 ~ Spatial +  E_PC1 + Com.Size.Gradient
          betas.LCBD ~ Spatial +  E_PC1 + Com.Size.Gradient
          '
smod1 = ' betas.LCBD ~ Spatial +  E_PC1 + Com.Size.Gradient
          '
smod1 = ' N1 ~ Spatial +  E_PC1 + Com.Size.Gradient
          '

smod1.fit <- sem(smod1,data=all_big_datz)
summary(smod1.fit,standardized=TRUE,rsq=T)
fitMeasures(smod1.fit)
modindices(smod1.fit)

#quick plot of path analysis
semPaths(smod1.fit, what='std', layout = "tree3", intercepts = FALSE, residuals = FALSE, 
         edge.label.cex=1.25, curvePivot = FALSE,  fade=FALSE, rotation = 2)


smod1 = ' diversity ~ Spatial +  E_PC1 + Com.Size.Gradient
          diversity =~ betas.LCBD +  N1 
          '

#space<-'pred.persistence =~connect.per + meta.size +total.vol '
smod1.fit <- sem(smod1,data=all_big_datz)
summary(smod1.fit,standardized=TRUE,rsq=T)
fitMeasures(smod1.fit)
modindices(smod1.fit)

#quick plot of path analysis
semPaths(smod1.fit, what='std', layout = "tree3", intercepts = FALSE, residuals = FALSE, 
         edge.label.cex=1.25, curvePivot = FALSE,  fade=FALSE, rotation = 2)


smod1 = ' N1 ~ Head.river.dist +  River.dist.lake
          betas.LCBD ~  Head.river.dist +  River.dist.lake 
          '

#space<-'pred.persistence =~connect.per + meta.size +total.vol '
smod1.fit <- sem(smod1,data=dd_specie)
summary(smod1.fit,standardized=TRUE,rsq=T)
fitMeasures(smod1.fit)
modindices(smod1.fit)

#quick plot of path analysis
semPaths(smod1.fit, what='std', layout = "tree3", intercepts = FALSE, residuals = FALSE, 
         edge.label.cex=1.25, curvePivot = FALSE,  fade=FALSE, rotation = 2)

smod1 = ' N1 ~ Head.river.dist +  River.dist.lake
          '
smod1 = 'betas.LCBD ~  Head.river.dist +  River.dist.lake '
smod1.fit <- sem(smod1,data=dd_specie)
summary(smod1.fit,standardized=TRUE,rsq=T)
fitMeasures(smod1.fit)
modindices(smod1.fit)

#quick plot of path analysis
semPaths(smod1.fit, what='std', layout = "tree3", intercepts = F, residuals = FALSE, 
         edge.label.cex=1.25, curvePivot = FALSE,  fade=FALSE, rotation = 2)


smod1 = ' diversity ~ Head.river.dist +  River.dist.lake
           diversity =~ betas.LCBD +  N1 '
smod1.fit <- sem(smod1,data=dd_specie)
summary(smod1.fit,standardized=TRUE,rsq=T)
fitMeasures(smod1.fit)
modindices(smod1.fit)

#quick plot of path analysis
semPaths(smod1.fit, what='std', layout = "tree3", intercepts = FALSE, residuals = FALSE, 
         edge.label.cex=1.25, curvePivot = FALSE,  fade=FALSE, rotation = 2)



########################################################################################################################################################################
#Old figs
new_labels <- c( "Com.Size.Gradient" = "Community Size Gradient","S_PC1" = "Spatial Gradient", "E_PC1" = "Environmental Gradient")

e6<-all_big_dat%>%
  pivot_longer(c(Com.Size.Gradient,S_PC1,E_PC1) , names_to = "key", values_to = "value")
summary(e2)

g1<-e6%>% 
  ggplot(aes(x = value, y =betas.LCBD, colour=O.NET )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point(data = filter(e6, O.NET =="BUBBS" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e6, O.NET =="CASCADE" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e6, O.NET =="EVO" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e6, O.NET =="ROCK" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e6, O.NET =="BUBBS" & key=="E_PC1"), shape=19)+
  geom_point(data = filter(e6, O.NET =="CASCADE" & key=="E_PC1"), shape=19)+
  geom_point(data = filter(e6, O.NET =="EVO" & key=="E_PC1"), shape=1)+
  geom_point(data = filter(e6, O.NET =="ROCK" & key=="E_PC1"), shape=1)+
  geom_point(data = filter(e6, O.NET =="BUBBS" & key=="Com.Size.Gradient"), shape=1)+
  geom_point(data = filter(e6, O.NET =="CASCADE" & key=="Com.Size.Gradient"), shape=19)+
  geom_point(data = filter(e6, O.NET =="EVO" & key=="Com.Size.Gradient"), shape=1)+
  geom_point(data = filter(e6, O.NET =="ROCK" & key=="Com.Size.Gradient"), shape=1)+
  #geom_smooth(data = filter(e6, O.NET =="BUBBS" & key=="S_PC1"),method="lm")+
  #geom_smooth(data=filter(e6,O.NET=="CASCADE"& key=="S_PC1"), method = "lm")+
  #geom_smooth(data=filter(e6,O.NET=="EVO"& key=="S_PC1"), method = "lm")+
  #geom_smooth(data=filter(e6,O.NET=="ROCK"& key=="S_PC1"), method = "lm")+
  geom_smooth(data=filter(e6,O.NET=="BUBBS"& key=="E_PC1"), method = "lm")+
  #geom_smooth(data=filter(e6,O.NET=="EVO"& key=="E_PC1"), method = "lm")+
  geom_smooth(data=filter(e6,O.NET=="CASCADE"& key=="E_PC1"), method = "lm")+
  geom_smooth(data=filter(e6,O.NET=="ROCK"& key=="E_PC1"), method = "lm")+
  #geom_smooth(data = filter(e6, O.NET =="BUBBS" & key=="Com.Size.Gradient"),method="lm")+
  geom_smooth(data=filter(e6,O.NET=="CASCADE"& key=="Com.Size.Gradient"), method = "lm")+
  #geom_smooth(data=filter(e6,O.NET=="EVO"& key=="Com.Size.Gradient"), method = "lm")+
  #geom_smooth(data=filter(e6,O.NET=="ROCK"& key=="Com.Size.Gradient"), method = "lm")+
  facet_grid(O.NET~key,scales="free_x", labeller =labeller(key= new_labels)) + #axis.title.x+
  ylab("Local Contribution to Beta Diversity (LCBD)")+
  xlab("   Log Com. Size         Env. Gradeint (E_PC1)     Spa. Gradient (S_PC1)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

g4<-e6%>% 
  ggplot(aes(x = value, y =N1, colour=O.NET )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point(data = filter(e6, O.NET =="BUBBS" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e6, O.NET =="CASCADE" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e6, O.NET =="EVO" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e6, O.NET =="ROCK" & key=="S_PC1"), shape=1)+
  geom_point(data = filter(e6, O.NET =="BUBBS" & key=="E_PC1"), shape=19)+
  geom_point(data = filter(e6, O.NET =="CASCADE" & key=="E_PC1"), shape=19)+
  geom_point(data = filter(e6, O.NET =="EVO" & key=="E_PC1"), shape=1)+
  geom_point(data = filter(e6, O.NET =="ROCK" & key=="E_PC1"), shape=1)+
  geom_point(data = filter(e6, O.NET =="BUBBS" & key=="Com.Size.Gradient"), shape=1)+
  geom_point(data = filter(e6, O.NET =="CASCADE" & key=="Com.Size.Gradient"), shape=19)+
  geom_point(data = filter(e6, O.NET =="EVO" & key=="Com.Size.Gradient"), shape=1)+
  geom_point(data = filter(e6, O.NET =="ROCK" & key=="Com.Size.Gradient"), shape=1)+
  #geom_smooth(data = filter(e6, O.NET =="BUBBS" & key=="S_PC1"),method="lm")+
  #geom_smooth(data=filter(e6,O.NET=="CASCADE"& key=="S_PC1"), method = "lm")+
  #geom_smooth(data=filter(e6,O.NET=="EVO"& key=="S_PC1"), method = "lm")+
  #geom_smooth(data=filter(e6,O.NET=="ROCK"& key=="S_PC1"), method = "lm")+
  geom_smooth(data=filter(e6,O.NET=="BUBBS"& key=="E_PC1"), method = "lm")+
  #geom_smooth(data=filter(e6,O.NET=="EVO"& key=="E_PC1"), method = "lm")+
  geom_smooth(data=filter(e6,O.NET=="CASCADE"& key=="E_PC1"), method = "lm")+
  geom_smooth(data=filter(e6,O.NET=="ROCK"& key=="E_PC1"), method = "lm")+
  #geom_smooth(data = filter(e6, O.NET =="BUBBS" & key=="Com.Size.Gradient"),method="lm")+
  geom_smooth(data=filter(e6,O.NET=="CASCADE"& key=="Com.Size.Gradient"), method = "lm")+
  #geom_smooth(data=filter(e6,O.NET=="EVO"& key=="Com.Size.Gradient"), method = "lm")+
  #geom_smooth(data=filter(e6,O.NET=="ROCK"& key=="Com.Size.Gradient"), method = "lm")+
  facet_grid(O.NET~key,scales="free_x", labeller =labeller(key= new_labels)) + #axis.title.x+
  ylab("Species Diversity")+
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
plot_grid(g1,g4, ncol=2)
