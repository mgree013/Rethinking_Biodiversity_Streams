setwd("~/Users/matthewdouglasgreen/Dropbox/Manuscipts/L-S Biodviersity Streams_RCC_SDH") 
getwd()
library(vegan)
library(tidyverse)
library(dplyr)
library(adespatial)
library(MuMIn)
library(gridExtra)
library(gtable)
library(grid)

species<-read.csv(file = "Analysis/Stream_frame/sp.density.update.12.28.19.csv", row.name=1)
env <-read.csv(file = "Analysis/Stream_frame/dave.matt.env.all.csv", row.name=1)
summary(env)


#resolved SPecies list now
species<-species%>%dplyr::select(-c(Chironomidae,Arachnida,Nematomorpha,Oligochaeta,Ostracoda,Turbellaria,Euhirudinea))


#Transform Env variables
#plot(density(specie$Head.river.dist))
#shapiro.test(env$River.dist.lake)

env<-env%>%mutate(Euc.dist.lake=log(1+Euc.dist.lake),River.dist.lake=log(1+River.dist.lake),Elevation=log(1+Elevation),Head.river.dist=log(1+Head.river.dist))


#Calcualte diversity and bind with envrioemtnal data, remove network if necessary

diversity<-species%>%
  #group_by(Site,Network)%>%
  transmute(N0=rowSums(species > 0),H= diversity(species),N1 =exp(H),N2 =diversity(species, "inv"),J= H/log(N0),E10= (N1/N0),E20= (N2/N0),Com.Size=rowSums(species)) #,betas.LCBD=beta.div(species, method="hellinger",sqrt.D=TRUE)$LCBD ,betas.LCBD.p=beta.div(species, method="chord",sqrt.D=TRUE)$p.LCBD )



specie%>%
  gather(N0, N1, E10, betas.LCBD, key = "var", value = "value") %>% 
  ggplot(aes(x =log(Head.river.dist+1), y = value)) + #remove , fill=Network and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap( ~var  , scales = "free") +
  theme_bw()
combine<-cbind(species,env$O.NET)

kern<-subset(combine, `env$O.NET`=="KERN")
kern.var<-kern[, colSums(kern != 0) > 0]
kern.beta<-beta.div(kern.var[1:61], method="hellinger",sqrt.D=TRUE)
kern.beta.comp<-beta.div.comp(kern.var[1:61], coef = "S", quant=T)

casc<-subset(combine, `env$O.NET`=="CASCADE")
casc.var<-casc[, colSums(casc != 0) > 0]
casc.beta<-beta.div(casc.var[1:44], method="hellinger",sqrt.D=TRUE)
casc.beta.comp<-beta.div.comp(casc.var[1:44], coef = "S", quant=T)

evo<-subset(combine, `env$O.NET`=="EVO")
evo.var<-evo[, colSums(evo != 0) > 0]
evo.beta<-beta.div(evo.var[1:36], method="hellinger",sqrt.D=TRUE)
evo.beta.comp<-beta.div.comp(evo.var[1:36], coef = "S", quant=T)

bubb<-subset(combine, `env$O.NET`=="BUBBS")
bubb.var<-bubb[, colSums(bubb != 0) > 0]
bubb.beta<-beta.div(bubb.var[1:83], method="hellinger",nperm=999,sqrt.D=TRUE)
bubb.beta.comp<-beta.div.comp(bubb.var[1:83], coef = "S", quant=T)

young<-subset(combine, `env$O.NET`=="YOUNG")
young.var<-young[, colSums(young != 0) > 0]
young.beta<-beta.div(young.var[1:29], nperm=999,method="hellinger",sqrt.D=TRUE)
young.beta.comp<-beta.div.comp(young.var[1:29], coef = "S", quant=T)
betas.LCBD<-c(kern.beta$LCBD,casc.beta$LCBD,evo.beta$LCBD,bubb.beta$LCBD,young.beta$LCBD)

all<-cbind(diversity,betas.LCBD, env)
specie<-all%>%filter( O.NET != "KERN" & O.NET != "YOUNG")
###########################################################################################

#NEW Stuff
ke<-as.data.frame(kern.beta$SCBD)
ke<-ke%>%rename(SCBD=`kern.beta$SCBD`)%>%rownames_to_column(var="species")%>%add_column(Net="KERN")

ca<-as.data.frame(casc.beta$SCBD)
ca<-ca%>%rename(SCBD=`casc.beta$SCBD`)%>%rownames_to_column(var="species")%>%add_column(Net="CASCADE")

e<-as.data.frame(evo.beta$SCBD)
e<-e%>%rename(SCBD=`evo.beta$SCBD`)%>%rownames_to_column(var="species")%>%add_column(Net="EVO")

b<-as.data.frame(bubb.beta$SCBD)
b<-b%>%rename(SCBD=`bubb.beta$SCBD`)%>%rownames_to_column(var="species")%>%add_column(Net="BUBBS")
  
all.scbd<-rbind(ca,ke,e,b)
all.scbd%>%
  ggplot(aes(x=species,y=SCBD))+
  geom_boxplot()+ 
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_bw()

all.scbd%>%
  ggplot(aes(x=species,y=SCBD, colour=Net))+
  geom_point()+ 
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_bw()
  
  
env.evo<-env%>%filter(Network=="EVO")
distance<-matrix(c(env.evo$Head.river.dist,env.evo$Elevation), nrow=21, ncol=2)
dist<-dist(distance, method="euclidean")

env2 = env.evo%>%dplyr::select(-c(Catch.Area.m,MIXED_FOREST,WOODY_WETLANDS,New.Net,HERBACEUOUS,
                                  O.NET ,Network,Lat,Lon,Local.Dist,Euc.dist.lake,River.dist.lake,
                                  Head.river.dist,BARREN_LAND, EVERGREEN_FOREST,EMERGENT_HERBACEUOUS_WETLANDS))
head(env2)
Epca = prcomp(env2, scale.=T, tol = 0.1)

biplot(Epca)
summary(Epca)
Epca$rotation
Epca$x
pc_scores <- data.frame(Epca$x[,1:4])
head(pc_scores)

envdist<-vegdist(pc_scores, method = "euclidean")

plot(dist,evo.beta.comp$repl)
plot(dist,evo.beta.comp$rich)
plot(dist,evo.beta.comp$D)

plot(envdist,evo.beta.comp$repl)
plot(envdist,evo.beta.comp$rich)
plot(envdist,evo.beta.comp$D)

env.bubbs<-env%>%filter(O.NET=="BUBBS")

distance<-matrix(c(env.bubbs$Lat,env.bubbs$Lon), nrow=26, ncol=2)
dist<-dist(distance, method="euclidean")
plot(dist,bubb.beta.comp$repl)
plot(dist,bubb.beta.comp$rich)
plot(dist,bubb.beta.comp$D)
##########################################################################################################
specie%>%
  gather(N0, N1, E10, betas.LCBD, key = "var", value = "value") %>% 
  ggplot(aes(x =log(Head.river.dist+1), y = value)) + #remove , fill=Network and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap( ~var  , scales = "free") +
  theme_bw()

specie%>%
  gather(N0, N1, E10, betas.LCBD, key = "var", value = "value") %>% 
  ggplot(aes(x = log(Size.net.dist), y = value,colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm")+
  facet_grid( var ~ O.NET, scales = "free") +
  theme_bw()

specie%>%
  gather(N0,key = "var", value = "value") %>% 
  summarise(N0=sum)%>%
  ggplot(aes(x = Network, y = value,colour=Network)) + #remove , fill=Network and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_boxplot()
  facet_grid( var ~ Network, scales = "free") +
  theme_bw()
  
specieaa<-cbind(env$O.NET,env$Network,species)

d<-species%>%
  add_column(O.NET=env$O.NET, .before="Acentrella")%>%
  add_column(Network=env$Network, .before="O.NET")%>%
  filter( Network != "UP.KERN" & Network != "YOUNG")%>%
  select(-c(Network))%>%
  group_by(O.NET)%>%
  summarise_all(sum) %>%
  ungroup %>% 
  transmute(O.NET,Richness = apply(.[4:(ncol(.)-1)] > 0, 1, sum))


d<-species%>%
  add_column(O.NET=env$O.NET, .before="Acentrella")%>%
  add_column(Network=env$Network, .before="O.NET")%>%
  filter( Network != "UP.KERN" & Network != "YOUNG")%>%
  select(-c(Network))

bubb<-d%>%filter( O.NET == "BUBBS")%>%select(-c(O.NET))
cas<-d%>%filter( O.NET == "CASCADE")%>%select(-c(O.NET))
evo<-d%>%filter( O.NET == "EVO")%>%select(-c(O.NET))
ker<-d%>%filter( O.NET == "KERN")%>%select(-c(O.NET))


b.betas.LCBD=beta.div(bubb, method="hellinger",sqrt.D=TRUE)
c.betas.LCBD=beta.div(cas, method="hellinger",sqrt.D=TRUE)
e.betas.LCBD=beta.div(evo, method="hellinger",sqrt.D=TRUE)
k.betas.LCBD=beta.div(ker, method="hellinger",sqrt.D=TRUE)

b<-subset(d[,2:149], O.NET=="BUBBS")
########################################################################################################################################
##Final Figure 1
new_labels <- c( "Head.river.dist" = "Headwater River Distance", "River.dist.lake" = "River Distance from Lakes")

e2<-specie%>%
  pivot_longer(c(Head.river.dist,River.dist.lake) , names_to = "key", values_to = "value")
summary(e2)

#PCBD
g1<-e2%>% 
  ggplot(aes(x = value, y =betas.LCBD, colour=O.NET )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="Head.river.dist"), shape=19)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="Head.river.dist"), shape=19)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="Head.river.dist"), shape=19)+
  geom_point(data = filter(e2, O.NET =="KERN" & key=="Head.river.dist"), shape=19)+
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="River.dist.lake"), shape=1)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="River.dist.lake"), shape=19)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="River.dist.lake"), shape=1)+
  geom_point(data = filter(e2, O.NET =="KERN" & key=="River.dist.lake"), shape=1)+
  geom_smooth(data = filter(e2, O.NET =="BUBBS" & key=="Head.river.dist"),method="lm")+
  geom_smooth(data=filter(e2,O.NET=="BUBBS"& key=="Head.river.dist"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="Head.river.dist"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="EVO"& key=="Head.river.dist"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="KERN"& key=="Head.river.dist"),  method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="BUBBS"& key=="River.dist.lake"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="River.dist.lake"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="EVO"& key=="River.dist.lake"),  method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="KERN"& key=="River.dist.lake"),  method = "lm")+
  facet_grid(O.NET~key,scales="free_x", labeller =labeller(key= new_labels)) + #axis.title.x+
  ylab("LCBD")+
  xlab("Log Headwater Distance (m)           Log Distance from Lake (m)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

#Richness
g2<-e2%>% 
  ggplot(aes(x = value, y =N0, colour=O.NET )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="Head.river.dist"), shape=1)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="Head.river.dist"), shape=19)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="Head.river.dist"), shape=19)+
  geom_point(data = filter(e2, O.NET =="KERN" & key=="Head.river.dist"), shape=19)+
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="River.dist.lake"), shape=19)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="River.dist.lake"), shape=1)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="River.dist.lake"), shape=1)+
  geom_point(data = filter(e2, O.NET =="KERN" & key=="River.dist.lake"), shape=1)+
  #geom_smooth(data = filter(e2, O.NET =="BUBBS" & key=="Head.river.dist"),method="lm")+
  geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="Head.river.dist"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="EVO"& key=="Head.river.dist"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="KERN"& key=="Head.river.dist"),  method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="BUBBS"& key=="River.dist.lake"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="River.dist.lake"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="EVO"& key=="River.dist.lake"),  method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="KERN"& key=="River.dist.lake"),  method = "lm")+
  facet_grid(  O.NET~key,scales="free_x", labeller =labeller(key= new_labels)) + #axis.title.x+
  ylab("Species Richness")+
  xlab("Log Headwater Distance (m)           Log Distance from Lake (m)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

#EVENNess
e3<-specie%>%
  pivot_longer(c(Head.river.dist,River.dist.lake) , names_to = "key", values_to = "value")
summary(e2)

g3<-e3%>% 
  ggplot(aes(x = value, y =E10, colour=O.NET )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="Head.river.dist"), shape=19)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="Head.river.dist"), shape=1)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="Head.river.dist"), shape=1)+
  geom_point(data = filter(e2, O.NET =="KERN" & key=="Head.river.dist"), shape=1)+
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="River.dist.lake"), shape=19)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="River.dist.lake"), shape=19)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="River.dist.lake"), shape=19)+
  geom_point(data = filter(e2, O.NET =="KERN" & key=="River.dist.lake"), shape=1)+
  geom_smooth(data=filter(e3,O.NET=="BUBBS"& key=="Head.river.dist"), method = "lm")+
  #geom_smooth(data=filter(e3,O.NET=="CASCADE"& key=="Head.river.dist"), method = "lm")+
  #geom_smooth(data=filter(e3,O.NET=="EVO"& key=="Head.river.dist"), method = "lm")+
  #geom_smooth(data=filter(e3,O.NET=="KERN"& key=="Head.river.dist"),  method = "lm")+
  geom_smooth(data=filter(e3,O.NET=="BUBBS"& key=="River.dist.lake"), method = "lm")+
  geom_smooth(data=filter(e3,O.NET=="CASCADE"& key=="River.dist.lake"), method = "lm")+
  geom_smooth(data=filter(e3,O.NET=="EVO"& key=="River.dist.lake"),  method = "lm")+
  #geom_smooth(data=filter(e3,O.NET=="KERN"& key=="River.dist.lake"),  method = "lm")+
  facet_grid(  O.NET~key,scales="free_x", labeller =labeller(key= new_labels)) + #axis.title.x+
  ylab("Species Evenness")+
  xlab("Log Headwater Distance (m)           Log Distance from Lake (m)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

plot_grid(g1,g2,g3, ncol=3)
############################################################################################################################################################################################
#GLMS
mod1<-glm(betas.LCBD~(River.dist.lake),family=gaussian(link = "identity"),data=specie)
mod2<-glm(betas.LCBD~(Head.river.dist),family=gaussian(link = "identity"),data=specie)
mod3<-glm(betas.LCBD~(Head.river.dist)*O.NET,family=gaussian(link = "identity"),data=specie)
mod4<-glm(betas.LCBD~(River.dist.lake)*O.NET,family=gaussian(link = "identity"),data=specie)
mod5<-glm(betas.LCBD~(River.dist.lake)*Head.river.dist,family=gaussian(link = "identity"),data=specie)
mod6<-glm(betas.LCBD~(River.dist.lake)*Head.river.dist*O.NET,family=gaussian(link = "identity"),data=specie)
nullmod<-glm(betas.LCBD~1+ O.NET,family=gaussian(link = "identity"),data=specie)
null<-glm(betas.LCBD~1,family=gaussian(link = "identity"),data=specie)

summary(mod4)

ggplot(specie, aes(x = Head.river.dist, y = betas.LCBD)) +
  geom_point(size = 4, aes(fill = O.NET), shape = 21) +
  geom_line(aes(y = predict(mod2, specie),
                colour = "log-log", linetype = "log-log")) +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  theme_bw()

ggplot(specie, aes(x = River.dist.lake, y = betas.LCBD)) +
  geom_point(size = 4, aes(fill = O.NET), shape = 21) +
  geom_line(aes(y = predict(mod1, specie),
                colour = "log-log", linetype = "log-log")) +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  theme_bw()

mod1<-glm(N0~(River.dist.lake),family=poisson(link = "log"),data=specie)
mod2<-glm(N0~(Head.river.dist),family=poisson(link = "log"),data=specie)
mod3<-glm(N0~(Head.river.dist)*O.NET,family=poisson(link = "log"),data=specie)
mod4<-glm(N0~(River.dist.lake)*O.NET,family=poisson(link = "log"),data=specie)
mod5<-glm(N0~(River.dist.lake)*Head.river.dist,family=poisson(link = "log"),data=specie)
mod6<-glm(N0~(River.dist.lake)*Head.river.dist*O.NET,family=poisson(link = "log"),data=specie)
nullmod<-glm(N0~1+ O.NET,family=poisson(link = "log"),data=specie)
null<-glm(N0~1,family=poisson(link = "log"),data=specie)

ggplot(specie, aes(x = Head.river.dist, y = N0)) +
  geom_point(size = 4, aes(fill = O.NET), shape = 21) +
  geom_line(aes(y = predict(mod2, specie),
                colour = "log-log", linetype = "log-log")) +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  theme_bw()

ggplot(specie, aes(x = River.dist.lake, y = N0)) +
  geom_point(size = 4, aes(fill = O.NET), shape = 21) +
  geom_line(aes(y = predict(mod1, specie),
                colour = "log-log")) +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  theme_bw()


mod1<-glm(E10~(River.dist.lake),family=gaussian(link = "identity"),data=specie)
mod2<-glm(E10~(Head.river.dist),family=gaussian(link = "identity"),data=specie)
mod3<-glm(E10~(Head.river.dist)*O.NET,family=gaussian(link = "identity"),data=specie)
mod4<-glm(E10~(River.dist.lake)*O.NET,family=gaussian(link = "identity"),data=specie)
mod5<-glm(E10~(River.dist.lake)*Head.river.dist,family=gaussian(link = "identity"),data=specie)
mod6<-glm(E10~(River.dist.lake)*Head.river.dist*O.NET,family=gaussian(link = "identity"),data=specie)
nullmod<-glm(E10~1+ O.NET,family=gaussian(link = "identity"),data=specie)
null<-glm(E10~1,family=gaussian(link = "identity"),data=specie)


reported.table2 <- bbmle::AICtab(mod1,mod2,mod3,mod4,mod5,mod6,nullmod, null,weights = TRUE, sort = TRUE)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod5, null,weights = TRUE, sort = TRUE)

form<-c(mod1$formula[[1]],mod2$formula,mod3$formula,mod4$formula,mod5$formula,mod6$formula,nullmod$formula,null$formula)
export<-data.frame((reported.table2))

pseudoR1 <- ((mod1$null.deviance-mod1$deviance)/mod1$null.deviance)
pseudoR2 <- ((mod2$null.deviance-mod2$deviance)/mod2$null.deviance)
pseudoR3 <- ((mod3$null.deviance-mod3$deviance)/mod3$null.deviance)
pseudoR4 <- ((mod4$null.deviance-mod4$deviance)/mod4$null.deviance)
pseudoR5 <- ((mod5$null.deviance-mod5$deviance)/mod5$null.deviance)
pseudoR6 <- ((mod6$null.deviance-mod6$deviance)/mod6$null.deviance)
pseudoRnullmod <- ((nullmod$null.deviance-nullmod$deviance)/nullmod$null.deviance)
pseudoRnull <- ((null$null.deviance-null$deviance)/null$null.deviance)

psuedoR2<-rbind(pseudoR1,pseudoR2,pseudoR3,pseudoR4,pseudoR5,pseudoR6,pseudoRnullmod,pseudoRnull)


####TABle2 1b
mod0<-glm(betas.LCBD~(Head.river.dist)*(River.dist.lake),family=gaussian(link = "identity"),data=specie)
mod1<-glm(betas.LCBD~(River.dist.lake),family=gaussian(link = "identity"),data=specie)
mod2<-glm(betas.LCBD~(Head.river.dist),family=gaussian(link = "identity"),data=specie)
null<-glm(betas.LCBD~1,family=gaussian(link = "identity"),data=specie)

mod0<-glm(N0~(River.dist.lake)*(Head.river.dist),family=poisson(link = "identity"),data=specie)
mod1<-glm(N0~(River.dist.lake),family=poisson(link = "log"),data=specie)
mod2<-glm(N0~(Head.river.dist),family=poisson(link = "log"),data=specie)
null<-glm(N0~1,family=poisson(link = "log"),data=specie)

mod0<-glm(E10~(River.dist.lake)*(Head.river.dist),family=gaussian(link = "identity"),data=specie)
mod1<-glm(E10~(River.dist.lake),family=gaussian(link = "identity"),data=specie)
mod2<-glm(E10~(Head.river.dist),family=gaussian(link = "identity"),data=specie)
null<-glm(E10~1,family=gaussian(link = "identity"),data=specie)

summary(mod3)
reported.table2 <- bbmle::AICtab(mod0,mod1,mod2, null,weights = TRUE, sort = FALSE)

pseudoR0 <- ((mod0$null.deviance-mod0$deviance)/mod0$null.deviance)
pseudoR1 <- ((mod1$null.deviance-mod1$deviance)/mod1$null.deviance)
pseudoR2 <- ((mod2$null.deviance-mod2$deviance)/mod2$null.deviance)
pseudoRnull <- ((null$null.deviance-null$deviance)/null$null.deviance)

psuedoR2<-rbind(pseudoR0,pseudoR1,pseudoR2,pseudoRnull)
summary(mod0)

library(interplot)
library(sjPlot)
devtools::install_github("strengejacke/sjmisc")
library(sjmisc)
library(sjlabelled)
library(ggiraphExtra)
interplot(m=mod0, var1="Head.river.dist",var2="River.dist.lake")
#LMER

mod0<-glmer(betas.LCBD~(River.dist.lake)*(Head.river.dist) +(1|O.NET),family=gaussian(link = "identity"),data=specie)
mod1<-glmer(betas.LCBD~(River.dist.lake)+(1|O.NET),family=gaussian(link = "identity"),data=specie)
mod2<-glmer(betas.LCBD~(Head.river.dist)+(1|O.NET),family=gaussian(link = "identity"),data=specie)
null<-glmer(betas.LCBD~1+(1|O.NET),family=gaussian(link = "identity"),data=specie)

summary(mod0)
sjPlot::sjp.glmer(mod0, type = "fe.pc")
plot_models(mod0,mod1,mod2, show.values = FALSE, show.p = FALSE, p.shape = TRUE)
mod0<-glm(betas.LCBD~(Head.river.dist)+O.NET,family=gaussian(link = "identity"),data=specie)

ggPredict(mod0,se=F,interactive=TRUE)
#Table 3:


#Figure 2 Opt 2

d.b1<-ggplot(specie,aes(x = River.dist.lake, y = betas.LCBD)) + geom_point()+
  #geom_smooth(method = "lm", se=F)+ 
  xlab(" Log Distance from Upstream Lakes (m)") +labs(y=(("Local Contribution to \u03B2-diversity (LCBD) ")))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

d.b2<-ggplot(specie,aes(x = (Head.river.dist), y = betas.LCBD)) + geom_point()+
  geom_smooth(method = "lm", se=F)+ 
  xlab(" Log Distance from Headwaters (m)") +labs(y=(("Local Contribution to \u03B2-diversity (LCBD) ")))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

d.r1<-ggplot(specie,aes(x = River.dist.lake, y = N0)) + geom_point()+
  #geom_smooth(method = "lm", se=F)+ 
  xlab(" Log Distance from Upstream Lakes (m)") +labs(y=(("Species Richness")))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

d.r2<-ggplot(specie,aes(x = (Head.river.dist), y = N0)) + geom_point()+
  #geom_smooth(method = "lm", se=F)+ 
  xlab(" Log Distance from Headwaters (m)") +labs(y=(("Species Richness")))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

d.e1<-ggplot(specie,aes(x = River.dist.lake, y = E10)) + geom_point()+
  geom_smooth(method = "lm", se=F)+ xlab(" Log Distance from Upstream Lakes (m)") +labs(y=(("Species Evenness")))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

d.e2<-ggplot(specie,aes(x = Head.river.dist, y = E10)) + geom_point()+
 # geom_smooth(method = "lm", se=F)+ 
  xlab(" Log Distance from Headwaters (m)") +labs(y=(("Species Evenness")))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


plot_grid(d.b1,d.r1,d.e1,d.b2,d.r2,d.e2)

dd<-lm(betas.LCBD~Head.river.dist, data=specie)
summary(dd)
dd<-lm(betas.LCBD~River.dist.lake, data=specie)
summary(dd)
dd<-lm(N0~Head.river.dist, data=specie)
summary(dd)
dd<-lm(N0~River.dist.lake, data=specie)
summary(dd)
dd<-lm(E10~Head.river.dist, data=specie)
summary(dd)
dd<-lm(E10~River.dist.lake, data=specie)
summary(dd)

########################################################################################################################
##################################################################################################################################################################################################################################################################################################################
#beta diversity

ggplot(specie,aes(x = (Head.river.dist), y = N0)) + geom_point()+
  geom_smooth(method = "lm", se=F)+ xlab(" log Distance from Headwaters (m)") +labs(y=(("Local Contribution to \u03B2-diversity (LCBD)")))

specie%>%
  gather(E10, key = "var", value = "value") %>% 
  ggplot(aes(x = Head.river.dist, y = value,colour=O.NET)) + #remove , fill=Network and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm")+
  facet_grid(  ~ O.NET, scales = "free") +
  theme_bw()

mod0.1<-glm(betas.LCBD~(Head.river.dist),family=gaussian (link = "identity"),data=subset(specie, O.NET=="BUBBS"))
mod0.2<-glm(betas.LCBD~(River.dist.lake),family=gaussian (link = "identity"),data=subset(specie, O.NET=="BUBBS"))
mod0.3<-glm(betas.LCBD~1,family=gaussian (link = "identity"),data=subset(specie, O.NET=="BUBBS"))

mod1.1<-glm(betas.LCBD~(Head.river.dist),family=gaussian (link = "identity"),data=subset(specie, O.NET=="CASCADE"))
mod1.2<-glm(betas.LCBD~(River.dist.lake),family=gaussian (link = "identity"),data=subset(specie, O.NET=="CASCADE"))
mod1.3<-glm(betas.LCBD~1,family=gaussian (link = "identity"),data=subset(specie, O.NET=="CASCADE"))

mod2.1<-glm(betas.LCBD~(Head.river.dist),family=gaussian (link = "identity"),data=subset(specie, O.NET=="EVO"))
mod2.2<-glm(betas.LCBD~(River.dist.lake),family=gaussian (link = "identity"),data=subset(specie, O.NET=="EVO"))
mod2.3<-glm(betas.LCBD~(1),family=gaussian (link = "identity"),data=subset(specie, O.NET=="EVO"))
anova (mod0.1, mod0.2,mod0.3,test="Chisq")

mod3.1<-glm(betas.LCBD~(Head.river.dist),family=gaussian (link = "identity"),data=subset(specie, O.NET=="KERN"))
mod3.2<-glm(betas.LCBD~(River.dist.lake),family=gaussian (link = "identity"),data=subset(specie, O.NET=="KERN"))
mod3.3<-glm(betas.LCBD~(1),family=gaussian (link = "identity"),data=subset(specie, O.NET=="KERN"))

mod4.1<-glm(betas.LCBD~(Head.river.dist),family=gaussian (link = "identity"),data=subset(specie))
mod4.2<-glm(betas.LCBD~(River.dist.lake),family=gaussian (link = "identity"),data=subset(specie))
mod4.3<-glm(betas.LCBD~1,family=gaussian (link = "identity"),data=subset(specie))

anova (mod0.1,mod0.3,test="Chisq")
anova (mod0.2,mod0.3,test="Chisq")

anova (mod1.1, mod1.3,test="Chisq")
anova (mod1.2, mod1.3,test="Chisq")

anova (mod2.1, mod2.3,test="Chisq")
anova (mod2.2, mod2.3,test="Chisq")

anova (mod3.1, mod3.3,test="Chisq")
anova (mod3.2, mod3.3,test="Chisq")

anova (mod4.1, mod4.3,test="Chisq")
anova (mod4.2, mod4.3,test="Chisq")

pseudoR0.1 <- ((mod0.1$null.deviance-mod0.1$deviance)/mod0.1$null.deviance)
pseudoR0.2 <- ((mod0.2$null.deviance-mod0.2$deviance)/mod0.2$null.deviance)
pseudoR0.3 <- ((mod0.3$null.deviance-mod0.3$deviance)/mod0.3$null.deviance)

pseudoR1.1 <- ((mod1.1$null.deviance-mod1.1$deviance)/mod1.1$null.deviance)
pseudoR1.2 <- ((mod1.2$null.deviance-mod1.2$deviance)/mod1.2$null.deviance)
pseudoR1.3 <- ((mod1.3$null.deviance-mod1.3$deviance)/mod1.3$null.deviance)

pseudoR2.1 <- ((mod2.1$null.deviance-mod2.1$deviance)/mod2.1$null.deviance)
pseudoR2.2 <- ((mod2.2$null.deviance-mod2.2$deviance)/mod2.2$null.deviance)
pseudoR2.3 <- ((mod2.2$null.deviance-mod2.3$deviance)/mod2.3$null.deviance)

pseudoR3.1 <- ((mod3.1$null.deviance-mod3.1$deviance)/mod3.1$null.deviance)
pseudoR3.2 <- ((mod3.2$null.deviance-mod3.2$deviance)/mod3.2$null.deviance)
pseudoR3.3 <- ((mod3.3$null.deviance-mod3.3$deviance)/mod3.3$null.deviance)

pseudoR4.1 <- ((mod4.1$null.deviance-mod4.1$deviance)/mod4.1$null.deviance)
pseudoR4.2 <- ((mod4.2$null.deviance-mod4.2$deviance)/mod4.2$null.deviance)

reported.table0 <- bbmle::AICtab(mod0.1,mod0.2,mod0.3, weights = TRUE, sort = TRUE)
reported.table1 <- bbmle::AICtab(mod1.1,mod1.2,mod1.3, weights = TRUE, sort = TRUE)
reported.table2 <- bbmle::AICtab(mod2.1,mod2.2,mod2.3, weights = TRUE, sort = TRUE)
reported.table3 <- bbmle::AICtab(mod3.1,mod3.2,mod3.3, weights = TRUE, sort = TRUE)

#########################################################################################################################################################
#Species Richness

ggplot(specie,aes(x = (new), y = N0)) + geom_point()+
  geom_smooth(method = "lm", se=F)+ xlab("Distance from Headwaters (m)") +labs(y=(("Species Richness")))

specie%>%
  gather(N0, key = "var", value = "value") %>% 
  ggplot(aes(x = Head.river.dist, y = value,colour=O.NET)) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm")+
  facet_grid(  ~ O.NET, scales = "free") +
  theme_bw()

specie%>%
  gather(N0, key = "var", value = "value") %>% 
  ggplot(aes(x = River.dist.lake, y = value,colour=O.NET)) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm")+
  facet_grid(  ~ O.NET, scales = "free") +
  theme_bw()

mod0.1<-glm(N0~(Head.river.dist),family=poisson (link = "log"),data=subset(specie, O.NET=="BUBBS"))
mod0.2<-glm(N0~(River.dist.lake),family=poisson (link = "log"),data=subset(specie, O.NET=="BUBBS"))
mod0.3<-glm(N0~1,family=poisson (link = "log"),data=subset(specie, O.NET=="BUBBS"))

mod1.1<-glm(N0~(Head.river.dist),family=poisson (link = "log"),data=subset(specie, O.NET=="CASCADE"))
mod1.2<-glm(N0~(River.dist.lake),family=poisson (link = "log"),data=subset(specie, O.NET=="CASCADE"))
mod1.3<-glm(N0~1,family=poisson (link = "log"),data=subset(specie, O.NET=="CASCADE"))

mod2.1<-glm(N0~(Head.river.dist),family=poisson (link = "log"),data=subset(specie, O.NET=="EVO"))
mod2.2<-glm(N0~(River.dist.lake),family=poisson (link = "log"),data=subset(specie, O.NET=="EVO"))
mod2.3<-glm(N0~(1),family=poisson (link = "log"),data=subset(specie, O.NET=="EVO"))

mod3.1<-glm(N0~(Head.river.dist),family=poisson (link = "log"),data=subset(specie, O.NET=="KERN"))
mod3.2<-glm(N0~(River.dist.lake),family=poisson (link = "log"),data=subset(specie, O.NET=="KERN"))
mod3.3<-glm(N0~(1),family=poisson (link = "log"),data=subset(specie, O.NET=="KERN"))

mod4.1<-glm(N0~(Head.river.dist),family=poisson (link = "log"),data=subset(specie))
mod4.2<-glm(N0~(River.dist.lake),family=poisson (link = "log"),data=subset(specie))
mod4.3<-glm(N0~(1),family=poisson (link = "log"),data=subset(specie))

#ANovas compare models to null and could be ggod for p value r2 aaporoach
anova (mod0.1,mod0.3,test="Chisq")
anova (mod0.2,mod0.3,test="Chisq")

anova (mod1.1, mod1.3,test="Chisq")
anova (mod1.2, mod1.3,test="Chisq")

anova (mod2.1, mod2.3,test="Chisq")
anova (mod2.2, mod2.3,test="Chisq")

anova (mod3.1, mod3.3,test="Chisq")
anova (mod3.2, mod3.3,test="Chisq")

anova (mod4.1, mod4.3,test="Chisq")
anova (mod4.2, mod4.3,test="Chisq")
anova.glm(mod4.1, mod4.3)
pseudoR0.1 <- ((mod0.1$null.deviance-mod0.1$deviance)/mod0.1$null.deviance)
pseudoR0.2 <- ((mod0.2$null.deviance-mod0.2$deviance)/mod0.2$null.deviance)
pseudoR0.3 <- ((mod0.3$null.deviance-mod0.3$deviance)/mod0.3$null.deviance)

pseudoR1.1 <- ((mod1.1$null.deviance-mod1.1$deviance)/mod1.1$null.deviance)
pseudoR1.2 <- ((mod1.2$null.deviance-mod1.2$deviance)/mod1.2$null.deviance)
pseudoR1.3 <- ((mod1.3$null.deviance-mod1.3$deviance)/mod1.3$null.deviance)

pseudoR2.1 <- ((mod2.1$null.deviance-mod2.1$deviance)/mod2.1$null.deviance)
pseudoR2.2 <- ((mod2.2$null.deviance-mod2.2$deviance)/mod2.2$null.deviance)
pseudoR2.3 <- ((mod2.2$null.deviance-mod2.3$deviance)/mod2.3$null.deviance)

pseudoR3.1 <- ((mod3.1$null.deviance-mod3.1$deviance)/mod3.1$null.deviance)
pseudoR3.2 <- ((mod3.2$null.deviance-mod3.2$deviance)/mod3.2$null.deviance)
pseudoR3.3 <- ((mod3.3$null.deviance-mod3.3$deviance)/mod3.3$null.deviance)

pseudoR4.1 <- ((mod4.1$null.deviance-mod4.1$deviance)/mod4.1$null.deviance)
pseudoR4.2 <- ((mod4.2$null.deviance-mod4.2$deviance)/mod4.2$null.deviance)

reported.table0 <- bbmle::AICtab(mod0.1,mod0.2,mod0.3, weights = TRUE, sort = TRUE)
reported.table1 <- bbmle::AICtab(mod1.1,mod1.2,mod1.3, weights = TRUE, sort = TRUE)
reported.table2 <- bbmle::AICtab(mod2.1,mod2.2,mod2.3, weights = TRUE, sort = TRUE)
reported.table3 <- bbmle::AICtab(mod3.1,mod3.2,mod3.3, weights = TRUE, sort = TRUE)

################################################################################################################################################################################################
#Species Evenneess

#IMP:Need to use resolved data for eveness

ggplot(specie,aes(x = (Euc.dist.lake), y = E10)) + geom_point()+
  geom_smooth(method = "lm", se=F)+ xlab("Distance from Upstream Lakes (m)") +labs(y=(("Species Evenness")))
summary(h)

specie%>%
  gather(E10, key = "var", value = "value") %>% 
  ggplot(aes(x = new, y = value,colour=O.NET)) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point()+
  geom_smooth(method = "lm")+
  facet_grid(  ~O.NET ,scales = "free") +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

mod0.1<-glm(E10~(Head.river.dist),family=gaussian (link = "identity"),data=subset(specie, O.NET=="BUBBS"))
mod0.2<-glm(E10~(Euc.dist.lake),family=gaussian (link = "identity"),data=subset(specie, O.NET=="BUBBS"))
mod0.3<-glm(E10~1,family=gaussian (link = "identity"),data=subset(specie, O.NET=="BUBBS"))

mod1.1<-glm(E10~(Head.river.dist),family=gaussian (link = "identity"),data=subset(specie, O.NET=="CASCADE"))
mod1.2<-glm(E10~(Euc.dist.lake),family=gaussian (link = "identity"),data=subset(specie, O.NET=="CASCADE"))
mod1.3<-glm(E10~1,family=gaussian (link = "identity"),data=subset(specie, O.NET=="CASCADE"))

mod2.1<-glm(E10~(Head.river.dist),family=gaussian (link = "identity"),data=subset(specie, O.NET=="EVO"))
mod2.2<-glm(E10~(Euc.dist.lake),family=gaussian (link = "identity"),data=subset(specie, O.NET=="EVO"))
mod2.3<-glm(E10~(1),family=gaussian (link = "identity"),data=subset(specie, O.NET=="EVO"))

mod3.1<-glm(E10~(Head.river.dist),family=gaussian (link = "identity"),data=subset(specie, O.NET=="KERN"))
mod3.2<-glm(E10~(Euc.dist.lake),family=gaussian (link = "identity"),data=subset(specie, O.NET=="KERN"))
mod3.3<-glm(E10~(1),family=gaussian (link = "identity"),data=subset(specie, O.NET=="KERN"))

mod4.1<-glm(E10~(Head.river.dist),family=gaussian (link = "identity"),data=subset(specie))
mod4.2<-glm(E10~(Euc.dist.lake),family=gaussian (link = "identity"),data=subset(specie))
mod4.3<-glm(E10~(1),family=gaussian (link = "identity"),data=subset(specie))


anova (mod0.1,mod0.3,test="Chisq")
anova (mod0.2,mod0.3,test="Chisq")

anova (mod1.1, mod1.3,test="Chisq")
anova (mod1.2, mod1.3,test="Chisq")

anova (mod2.1, mod2.3,test="Chisq")
anova (mod2.2, mod2.3,test="Chisq")

anova (mod3.1, mod3.3,test="Chisq")
anova (mod3.2, mod3.3,test="Chisq")

anova (mod4.1, mod4.3,test="Chisq")
anova (mod4.2, mod4.3,test="Chisq")

pseudoR0.1 <- ((mod0.1$null.deviance-mod0.1$deviance)/mod0.1$null.deviance)
pseudoR0.2 <- ((mod0.2$null.deviance-mod0.2$deviance)/mod0.2$null.deviance)
pseudoR0.3 <- ((mod0.3$null.deviance-mod0.3$deviance)/mod0.3$null.deviance)

pseudoR1.1 <- ((mod1.1$null.deviance-mod1.1$deviance)/mod1.1$null.deviance)
pseudoR1.2 <- ((mod1.2$null.deviance-mod1.2$deviance)/mod1.2$null.deviance)
pseudoR1.3 <- ((mod1.3$null.deviance-mod1.3$deviance)/mod1.3$null.deviance)

pseudoR2.1 <- ((mod2.1$null.deviance-mod2.1$deviance)/mod2.1$null.deviance)
pseudoR2.2 <- ((mod2.2$null.deviance-mod2.2$deviance)/mod2.2$null.deviance)
pseudoR2.3 <- ((mod2.2$null.deviance-mod2.3$deviance)/mod2.3$null.deviance)

pseudoR3.1 <- ((mod3.1$null.deviance-mod3.1$deviance)/mod3.1$null.deviance)
pseudoR3.2 <- ((mod3.2$null.deviance-mod3.2$deviance)/mod3.2$null.deviance)
pseudoR3.3 <- ((mod3.3$null.deviance-mod3.3$deviance)/mod3.3$null.deviance)

pseudoR4.1 <- ((mod4.1$null.deviance-mod4.1$deviance)/mod4.1$null.deviance)
pseudoR4.2 <- ((mod4.2$null.deviance-mod4.2$deviance)/mod4.2$null.deviance)

reported.table0 <- bbmle::AICtab(mod0.1,mod0.2,mod0.3, weights = TRUE, sort = TRUE)
reported.table1 <- bbmle::AICtab(mod1.1,mod1.2,mod1.3, weights = TRUE, sort = TRUE)
reported.table2 <- bbmle::AICtab(mod2.1,mod2.2,mod2.3, weights = TRUE, sort = TRUE)
reported.table3 <- bbmle::AICtab(mod3.1,mod3.2,mod3.3, weights = TRUE, sort = TRUE)

##################################################################################################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################################################
#Community Structure
#NMDS (Non-Metric Multideminsional Scaling)
#https://en.wikipedia.org/wiki/Multidimensional_scaling
alls<-cbind(species, env)
alles<-alls%>%filter( Network != "UP.KERN" & Network != "YOUNG")%>%drop_na()
dune.rel<-decostand(alles[,1:139],"hellinger")
dune.bray<-vegdist(dune.rel)
set.seed(134)
dune.nmds=metaMDS(dune.rel, k=2, try=100)
dune.nmds
stressplot(dune.nmds)

plot(dune.nmds,typ= "n", xlab = "NMDS Axis 1", ylab = "NMDS Axis 2")
text(dune.nmds$species[,1], dune.nmds$species[,2], rownames(dune.nmds$species), cex=0.6, col ="black")
points(dune.nmds$points[,1], dune.nmds$points[,2],  pch = 1) #I think it looks better with black dots and cant invert elevation so maybe you can run it on your end
ordihull(dune.nmds, groups=specie$Network, draw="polygon", label=T)
ordisurf(dune.nmds, alles$Head.river.dist, prioirty=,labcex=0.9, add = T,col="forestgreen")

adonis2(dune.bray ~ alles$Head.river.dist+alles$Elevation+alles$Network+alles$River.dist.lake, permutations = 999, method = "bray")
###############################################################################################################################################

#Beta.reg
#WEcould use this as an option...
library(betareg)
s<-betareg(betas.LCBD~(Head.river.dist),data=subset(specie))
s1<-betareg(betas.LCBD~(River.dist.lake),data=subset(specie))
s2<-betareg(betas.LCBD~(River.dist.lake)+(Head.river.dist),data=subset(specie))
s3<-betareg(betas.LCBD~1,data=subset(specie))

ggplot(specie, aes(x = Head.river.dist, y = betas.LCBD)) +
  geom_point(size = 4, aes(fill = O.NET), shape = 21) +
  geom_line(aes(y = predict(s, specie),
                colour = "log-log", linetype = "log-log")) +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  theme_bw()

summary(s)

s<-betareg(betas.LCBD~(Head.river.dist),data=subset(specie, O.NET=="BUBBS"))
s1<-betareg(betas.LCBD~(Head.river.dist),data=subset(specie, O.NET=="CASCADE"))
s2<-betareg(betas.LCBD~(Head.river.dist),data=subset(specie, O.NET=="EVO"))
s3<-betareg(betas.LCBD~(Head.river.dist),data=subset(specie, O.NET=="KERN"))

s<-betareg(betas.LCBD~(River.dist.lake),data=subset(specie, O.NET=="BUBBS"))
s1<-betareg(betas.LCBD~(River.dist.lake),data=subset(specie, O.NET=="CASCADE"))
s2<-betareg(betas.LCBD~(River.dist.lake),data=subset(specie, O.NET=="EVO"))
s3<-betareg(betas.LCBD~(River.dist.lake),data=subset(specie, O.NET=="KERN"))

s.1<-betareg(E10~(Head.river.dist),data=subset(specie))
s.11<-betareg(E10~(River.dist.lake),data=subset(specie))

s<-betareg(E10~(Head.river.dist),data=subset(specie, O.NET=="BUBBS"))
s1<-betareg(E10~(Head.river.dist),data=subset(specie, O.NET=="CASCADE"))
s2<-betareg(E10~(Head.river.dist),data=subset(specie, O.NET=="EVO"))
s3<-betareg(E10~(Head.river.dist),data=subset(specie, O.NET=="KERN"))

s4<-betareg(E10~(River.dist.lake),data=subset(specie, O.NET=="BUBBS"))
s5<-betareg(E10~(River.dist.lake),data=subset(specie, O.NET=="CASCADE"))
s6<-betareg(E10~(River.dist.lake),data=subset(specie, O.NET=="EVO"))
s7<-betareg(E10~(River.dist.lake),data=subset(specie, O.NET=="KERN"))

df.residual.betareg(s)
summary.betareg(s)
summary(s1)
summary(s2)
summary(s3)

ev.df <- data.frame=(df= c())

cat<-c(s.1$pseudo.r.squared,s.11$pseudo.r.squared,s$pseudo.r.squared,s4$pseudo.r.squared,s1$pseudo.r.squared,s5$pseudo.r.squared,s2$pseudo.r.squared,s6$pseudo.r.squared,s3$pseudo.r.squared,s7$pseudo.r.squared)
dog<-c(s.1$weights,s.11$pseudo.r.squared,s$pseudo.r.squared,s4$pseudo.r.squared,s1$pseudo.r.squared,s5$pseudo.r.squared,s2$pseudo.r.squared,s6$pseudo.r.squared,s3$pseudo.r.squared,s7$pseudo.r.squared)
s$df.null
s$df.residual 
s$loglik

library(AICcmodavg)
models<-list()
models[[1]]<-betareg(betas.LCBD~(Head.river.dist),data=subset(specie))
models[[2]]<-betareg(betas.LCBD~(River.dist.lake),data=subset(specie))

bubb.model<-list()
bubb.model[[1]]<-betareg(betas.LCBD~(Head.river.dist),data=subset(specie, O.NET=="BUBBS"))
bubb.model[[2]]<-betareg(betas.LCBD~(River.dist.lake),data=subset(specie, O.NET=="BUBBS"))

casc.model<-list()
casc.model[[1]]<-betareg(betas.LCBD~(Head.river.dist),data=subset(specie, O.NET=="CASCADE"))
casc.model[[2]]<-betareg(betas.LCBD~(River.dist.lake),data=subset(specie, O.NET=="CASCADE"))

evo.model<-list()
evo.model[[1]]<-betareg(betas.LCBD~(Head.river.dist),data=subset(specie, O.NET=="EVO"))
evo.model[[2]]<-betareg(betas.LCBD~(River.dist.lake),data=subset(specie, O.NET=="EVO"))

kern.model<-list()
kern.model[[1]]<-betareg(betas.LCBD~(Head.river.dist),data=subset(specie, O.NET=="KERN"))
kern.model[[2]]<-betareg(betas.LCBD~(River.dist.lake),data=subset(specie, O.NET=="KERN"))


mod<-aictab(cand.set = models, sort = FALSE)
b.mod<-aictab(cand.set = bubb.model, sort = FALSE)
c.mod<-aictab(cand.set = casc.model, sort = FALSE)
e.mod<-aictab(cand.set = evo.model, sort = FALSE)
k.mod<-aictab(cand.set = kern.model, sort = FALSE)

all<-rbind(mod,b.mod,c.mod,e.mod,k.mod)
summary(models[[1]])
parameters::p_value(models[[1]])
cat<-c(models[[1]]$pseudo.r.squared,s.11$pseudo.r.squared,s$pseudo.r.squared,s4$pseudo.r.squared,s1$pseudo.r.squared,s5$pseudo.r.squared,s2$pseudo.r.squared,s6$pseudo.r.squared,s3$pseudo.r.squared,s7$pseudo.r.squared)

models<-list()
models[[1]]<-betareg(E10~(Head.river.dist),data=subset(specie))
models[[2]]<-betareg(E10~(River.dist.lake),data=subset(specie))

bubb.model<-list()
bubb.model[[1]]<-betareg(E10~(Head.river.dist),data=subset(specie, O.NET=="BUBBS"))
bubb.model[[2]]<-betareg(E10~(River.dist.lake),data=subset(specie, O.NET=="BUBBS"))

casc.model<-list()
casc.model[[1]]<-betareg(E10~(Head.river.dist),data=subset(specie, O.NET=="CASCADE"))
casc.model[[2]]<-betareg(E10~(River.dist.lake),data=subset(specie, O.NET=="CASCADE"))

evo.model<-list()
evo.model[[1]]<-betareg(E10~(Head.river.dist),data=subset(specie, O.NET=="EVO"))
evo.model[[2]]<-betareg(E10~(River.dist.lake),data=subset(specie, O.NET=="EVO"))

kern.model<-list()
kern.model[[1]]<-betareg(E10~(Head.river.dist),data=subset(specie, O.NET=="KERN"))
kern.model[[2]]<-betareg(E10~(River.dist.lake),data=subset(specie, O.NET=="KERN"))


mod<-aictab(cand.set = models, sort = FALSE)
b.mod<-aictab(cand.set = bubb.model, sort = FALSE)
c.mod<-aictab(cand.set = casc.model, sort = FALSE)
e.mod<-aictab(cand.set = evo.model, sort = FALSE)
k.mod<-aictab(cand.set = kern.model, sort = FALSE)

all<-rbind(mod,b.mod,c.mod,e.mod,k.mod)


