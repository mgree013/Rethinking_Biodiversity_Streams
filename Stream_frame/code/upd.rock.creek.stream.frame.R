#paper title: Rethinking Biodiversity in Stream Ecology Frameworks
#Part 1 of Analysis: Stream Biodiversity Frameworks
#Author: Matthew Douglas Green
#Date: Sept 14,2021
#Data: https://datadryad.org/stash/dataset/doi:10.5061/dryad.2fqz612qw
########################################################################################################################
#load Libraries
Packages <- c("grid","gtable","gridExtra","dplyr","insight","performance","piecewiseSEM","betareg","adespatial","tidyverse","vegan","olsrr","semPlot","lavaan","lme4","vegan", "ggplot2", "tidyverse", "ape","MuMIn","adespatial", "betapart", "cowplot","glmmTMB")
lapply(Packages, library, character.only = TRUE)
library(ggbiplot)

##############################################################################################################
#Read data and clean
species<-read.csv(file = "Stream_frame/data/sp.density.update.12.28.19.csv", row.name=1)
env <-read.csv(file = "Stream_frame/data/dave.matt.env.all.csv", row.name=1)
summary(env)

#resolved Species list now
species<-species%>%dplyr::select(-c(Chironomidae,Arachnida,Nematomorpha,Oligochaeta,Ostracoda,Turbellaria,Euhirudinea))

#Re-scale environmental parameters
env<-env%>%mutate(Euc.dist.lake=log(1+Euc.dist.lake),River.dist.lake=log(1+River.dist.lake),Elevation=log(1+Elevation),Head.river.dist=log(1+Head.river.dist))

#########################################################################################################################
#Calculate diversity and bind with environmental data, remove network if necessary
diversity<-species%>%
  #group_by(Site,Network)%>%
  transmute(N0=rowSums(species > 0),H= diversity(species),N1 =exp(H),N2 =diversity(species, "inv"),J= H/log(N0),E10= (N1/N0),E20= (N2/N0),Com.Size=rowSums(species)) #,betas.LCBD=beta.div(species, method="hellinger",sqrt.D=TRUE)$LCBD ,betas.LCBD.p=beta.div(species, method="chord",sqrt.D=TRUE)$p.LCBD )


#Calculate Beta Diversity (LCBD) for each network
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

rock<-subset(combine, `env$O.NET`=="ROCK")
rock.var<-rock[, colSums(rock != 0) > 0]
rock.beta<-beta.div(rock.var[1:60], nperm=999,method="hellinger",sqrt.D=TRUE)
rock.beta.comp<-beta.div.comp(rock.var[1:60], coef = "S", quant=T)

betas.LCBD<-c(kern.beta$LCBD,casc.beta$LCBD,evo.beta$LCBD,bubb.beta$LCBD,young.beta$LCBD,rock.beta$LCBD)

all<-cbind(diversity,betas.LCBD, env)
specie<-all
###########################################################################################
#species Richness for each Network
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

############################################################################################################################################################################################
#1)GLMS
dd_specie<-specie%>%filter(River.dist.lake>0)%>%filter(Head.river.dist>2.5)#%>%filter(O.NET !="YOUNG")

qqnorm(log(specie$betas.LCBD))
qqline(log(specie$betas.LCBD))
wilcox.test(log(specie$betas.LCBD))

qqnorm((specie$N1))
qqline((specie$N1))
wilcox.test((specie$N1))

mod1<-glm(N1~(River.dist.lake),family=gaussian(link = "identity"),data=dd_specie)
mod2<-glm(N1~(Head.river.dist),family=gaussian(link = "identity"),data=dd_specie)
mod3<-glm(N1~(River.dist.lake)*Head.river.dist,family=gaussian(link = "identity"),data=dd_specie)
null<-glm(N1~1,family=gaussian(link = "identity"),data=dd_specie)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod5, null,weights = TRUE, sort = FALSE)
r2(mod1)
r2(mod2)
r2(mod3)
r2(null)
anova(mod1,mod2,mod5,mod6, null)

pseudoR1 <- ((mod1$null.deviance-mod1$deviance)/mod1$null.deviance)
pseudoR2 <- ((mod2$null.deviance-mod2$deviance)/mod2$null.deviance)
pseudoR0 <- ((mod5$null.deviance-mod5$deviance)/mod5$null.deviance)
pseudoRnull <- ((null$null.deviance-null$deviance)/null$null.deviance)



#2)Betareg Models
mod0<-betareg(betas.LCBD~River.dist.lake+O.NET, data=dd_specie)
mod01<-betareg(betas.LCBD~River.dist.lake, data=dd_specie)
reported.table2 <- bbmle::AICtab(mod01,mod0,weights = TRUE, sort = FALSE)

mod1<-betareg(betas.LCBD~River.dist.lake, data=dd_specie)
mod2<-betareg(betas.LCBD~Head.river.dist,data=dd_specie)
mod3<-betareg(betas.LCBD~Head.river.dist*River.dist.lake,data=dd_specie)
null<-betareg(betas.LCBD~1,data=dd_specie)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod3,null,weights = TRUE, sort = FALSE)
summary(mod01)
summary(mod0)
r2(mod1)
r2(mod2)
r2(mod3)
r2(null)

summary(mod1)
anova(mod1,mod2,mod5,mod6, null)

pseudoR1 <- ((mod1$null.deviance-mod1$deviance)/mod1$null.deviance)
pseudoR2 <- ((mod2$null.deviance-mod2$deviance)/mod2$null.deviance)
pseudoR0 <- ((mod5$null.deviance-mod5$deviance)/mod5$null.deviance)
pseudoRnull <- ((null$null.deviance-null$deviance)/null$null.deviance)
################################################################################################################################################################
#FIGURES

#Figure 4

d.b1<-dd_specie%>%
  #filter(Head.river.dist >3)%>%
  #filter(Network != "BUBBS" &O.NET != "KERN")%>%
  ggplot(aes(x = River.dist.lake, y = betas.LCBD)) + geom_point()+
  ggtitle("b)") +
  #stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)+
  #geom_smooth(method = "lm", se=F)+ 
  xlab(" Log Distance below Upstream Lakes (m)") +labs(y=(("Local Contribution to \u03B2-diversity (LCBD) ")))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

d.b2<-dd_specie%>%
  #filter(Head.river.dist >3)%>%
  #filter(Network != "BUBBS" &O.NET != "KERN")%>%
  ggplot(aes(x = (Head.river.dist), y = betas.LCBD)) + geom_point()+
  ggtitle("d)") +
  stat_smooth(method = glm,method.args = list(family = gaussian(link = "identity")))+
  # geom_smooth(method = "lm", se=F)+ 
  xlab(" Log Distance from Headwaters (m)") +labs(y=(("Local Contribution to \u03B2-diversity (LCBD) ")))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


d.e1<-dd_specie%>%
  #filter(Head.river.dist >3)%>%
  #filter(Network != "BUBBS" &O.NET != "KERN")%>%
  ggplot(aes(x = River.dist.lake, y = N1)) + geom_point()+
  ggtitle("a)") +
  stat_smooth(method = glm,method.args = list(family = gaussian(link = "identity")))+
  #geom_smooth(method = "lm", se=F)+ 
  xlab(" Log Distance below Upstream Lakes (m)") +labs(y=(("Shannon Diversity")))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

d.e2<-dd_specie%>%
  #filter(Head.river.dist >3)%>%
  #filter(Network != "BUBBS" &O.NET != "KERN")%>%
  ggplot(aes(x = Head.river.dist, y = N1)) + 
  ggtitle("c)") +
  geom_point()+
  stat_smooth(method = glm,method.args = list(family = gaussian(link = "identity")))+
  #geom_smooth(method = "lm", se=F)+ 
  xlab(" Log Distance from Headwaters (m)") +labs(y=(("Shannon Diversity")))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


plot_grid(d.b1,d.r1,d.e1,d.b2,d.r2,d.e2)

plot_grid(d.b1,d.e1,d.b2,d.e2)


########################################################################################################################
#########################################################################################################################
#########################################################################################################################################################################################
##Final Figure Appendix
new_labels <- c( "Head.river.dist" = "Headwater River Distance", "River.dist.lake" = "River Distance from Lakes")

e2<-specie%>%
  filter(Head.river.dist >4)%>%
  filter(New.Net != "U.KERN" & New.Net != "M.KERN" & New.Net != "BUBBS" & New.Net != "YOUNG")%>%
  pivot_longer(c(Head.river.dist,River.dist.lake) , names_to = "key", values_to = "value")
summary(e2)

#PCBD
g1<-e2%>% 
  ggplot(aes(x = value, y =betas.LCBD, colour=O.NET )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="Head.river.dist"), shape=19)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="Head.river.dist"), shape=19)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="Head.river.dist"), shape=19)+
  geom_point(data = filter(e2, O.NET =="KERN" & key=="Head.river.dist"), shape=19)+
  geom_point(data = filter(e2, O.NET =="ROCK" & key=="Head.river.dist"), shape=1)+
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="River.dist.lake"), shape=1)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="River.dist.lake"), shape=1)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="River.dist.lake"), shape=1)+
  geom_point(data = filter(e2, O.NET =="KERN" & key=="River.dist.lake"), shape=1)+
  geom_point(data = filter(e2, O.NET =="ROCK" & key=="River.dist.lake"), shape=19)+
  geom_smooth(data = filter(e2, O.NET =="BUBBS" & key=="Head.river.dist"),method="lm")+
  geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="Head.river.dist"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="EVO"& key=="Head.river.dist"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="KERN"& key=="Head.river.dist"),  method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="ROCK"& key=="Head.river.dist"),  method = "lm")+
  #stat_smooth(data=filter(e2,O.NET=="BUBBS"& key=="River.dist.lake"),method = "lm", formula = y ~ x + I(x^2), size = 1)+
  #stat_smooth(data=filter(e2,O.NET=="CASCADE"& key=="River.dist.lake"),method = "lm", formula = y ~ x + I(x^2), size = 1)+
  #stat_smooth(data=filter(e2,O.NET=="EVO"& key=="River.dist.lake"),method = "lm", formula = y ~ x + I(x^2), size = 1)+
  #stat_smooth(data=filter(e2,O.NET=="KERN"& key=="River.dist.lake"),method = "lm", formula = y ~ x + I(x^2), size = 1)+
  geom_smooth(data=filter(e2,O.NET=="ROCK"& key=="River.dist.lake"),method = "lm")+
  facet_grid(O.NET~key,scales="free", labeller =labeller(key= new_labels)) + #axis.title.x+
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
  geom_point(data = filter(e2, O.NET =="ROCK" & key=="River.dist.lake"), shape=1)+
  #geom_smooth(data = filter(e2, O.NET =="BUBBS" & key=="Head.river.dist"),method="lm")+
  geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="Head.river.dist"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="EVO"& key=="Head.river.dist"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="KERN"& key=="Head.river.dist"),  method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="BUBBS"& key=="River.dist.lake"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="River.dist.lake"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="EVO"& key=="River.dist.lake"),  method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="KERN"& key=="River.dist.lake"),  method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="ROCK"& key=="River.dist.lake"),  method = "lm")+
  facet_grid(  O.NET~key,scales="free_x", labeller =labeller(key= new_labels)) + #axis.title.x+
  ylab("Species Richness")+
  xlab("Log Headwater Distance (m)           Log Distance from Lake (m)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")

#Diversity
g4<-e2%>% 
  ggplot(aes(x = value, y =N1, colour=O.NET )) + #remove , fill=O.NET and see what the grpah looks like, are there t#F8766Dns that both entowrks share together
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="Head.river.dist"), shape=1)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="Head.river.dist"), shape=1)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="Head.river.dist"), shape=1)+
  geom_point(data = filter(e2, O.NET =="KERN" & key=="Head.river.dist"), shape=1)+
  geom_point(data = filter(e2, O.NET =="ROCK" & key=="Head.river.dist"), shape=19)+
  geom_point(data = filter(e2, O.NET =="BUBBS" & key=="River.dist.lake"), shape=1)+
  geom_point(data = filter(e2, O.NET =="CASCADE" & key=="River.dist.lake"), shape=1)+
  geom_point(data = filter(e2, O.NET =="EVO" & key=="River.dist.lake"), shape=19)+
  geom_point(data = filter(e2, O.NET =="KERN" & key=="River.dist.lake"), shape=1)+
  geom_point(data = filter(e2, O.NET =="ROCK" & key=="River.dist.lake"), shape=19)+
  #geom_smooth(data = filter(e2, O.NET =="BUBBS" & key=="Head.river.dist"),method="lm")+
  #geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="Head.river.dist"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="EVO"& key=="Head.river.dist"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="KERN"& key=="Head.river.dist"),  method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="ROCK"& key=="Head.river.dist"),  method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="BUBBS"& key=="River.dist.lake"), method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="CASCADE"& key=="River.dist.lake"), method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="EVO"& key=="River.dist.lake"),  method = "lm")+
  #geom_smooth(data=filter(e2,O.NET=="KERN"& key=="River.dist.lake"),  method = "lm")+
  geom_smooth(data=filter(e2,O.NET=="ROCK"& key=="River.dist.lake"),  method = "lm")+
  facet_grid(  O.NET~key,scales="free", labeller =labeller(key= new_labels)) + #axis.title.x+
  ylab("Species Diversity")+
  xlab("Log Headwater Distance (m)           Log Distance from Lake (m)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(legend.position = "none")


plot_grid(g1,g4, ncol=2)

######
sp_env<-species%>%rownames_to_column(var = "Site")
env_env<-env%>%rownames_to_column(var = "Site")
  
new_env<-sp_env%>%left_join(env_env, by="Site")%>%
  filter(log(River.dist.lake+1)>0)%>%
  filter(log(Pisidium) >0)
  

dog<-lm(log(Pisidium)~log(River.dist.lake+1), data=new_env)
summary(dog)
performance::r2(dog)
pseudoR1 <- ((dog$null.deviance-dog$deviance)/dog$null.deviance)

aaa<-new_env%>%
  filter(log(River.dist.lake+1)>0)%>%
  filter(Simulium >0)%>%
  #filter(Network != "BUBBS" &O.NET != "KERN")%>%
  ggplot(aes(x = log(River.dist.lake+1), y = log(Simulium))) + 
  geom_point()+
  stat_smooth(method = glm,method.args = list(family = gaussian(link = "identity")))+
  #geom_smooth(method = "lm", se=F)+ 
  ggtitle("a)") +
  annotate("text", x = 2, y = .95, label = "R^2 == 0.24", parse = TRUE) +
  xlab(" Log Distance from Upstream Lakes (m)") +labs(y=(("Simulium Log Density")))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

bbb<-new_env%>%
  filter(Pisidium >0)%>%
  #filter(Network != "BUBBS" &O.NET != "KERN")%>%
  ggplot(aes(x = log(River.dist.lake+1), y = log(Pisidium+1))) + 
  ggtitle("b)") +
  geom_point()+
  stat_smooth(method = glm,method.args = list(family = gaussian(link = "identity")))+
  annotate("text", x = 2, y = .95, label = "R^2 == 0.15", parse = TRUE) +
  xlab(" Log Distance from Upstream Lakes (m)") +labs(y=(("Pisidium Log Density")))+theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


plot_grid(aaa,bbb)
