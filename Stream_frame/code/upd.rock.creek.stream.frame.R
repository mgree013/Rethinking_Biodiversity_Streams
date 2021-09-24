#paper title: Rethinking Biodiversity in Stream Ecology Frameworks
#Part 1 of Analysis: Stream Biodiversity Frameworks
#Author: Matthew Douglas Green
#Date: Sept 14,2021
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
dd_specie<-specie%>%filter(River.dist.lake>0)%>%filter(Head.river.dist>2.5)

qqnorm(log(specie$betas.LCBD))
qqline(log(specie$betas.LCBD))
wilcox.test(log(specie$betas.LCBD))

qqnorm((specie$N1))
qqline((specie$N1))
wilcox.test((specie$N1))

mod1<-glm(N1~(River.dist.lake),family=gaussian(link = "identity"),data=dd_specie)
mod2<-glm(N1~(Head.river.dist),family=gaussian(link = "identity"),data=dd_specie)
mod5<-glm(N1~(River.dist.lake)*Head.river.dist,family=gaussian(link = "identity"),data=dd_specie)
null<-glm(N1~1,family=gaussian(link = "identity"),data=dd_specie)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod5, null,weights = TRUE, sort = FALSE)

anova(mod1,mod2,mod5,mod6, null)

pseudoR1 <- ((mod1$null.deviance-mod1$deviance)/mod1$null.deviance)
pseudoR2 <- ((mod2$null.deviance-mod2$deviance)/mod2$null.deviance)
pseudoR0 <- ((mod5$null.deviance-mod5$deviance)/mod5$null.deviance)
pseudoRnull <- ((null$null.deviance-null$deviance)/null$null.deviance)


mod1<-glm(betas.LCBD~(River.dist.lake),family=gaussian(),data=dd_specie)
mod2<-glm(betas.LCBD~(Head.river.dist),family=gaussian(),data=dd_specie)
mod5<-glm(betas.LCBD~River.dist.lake*Head.river.dist,family=gaussian(),data=dd_specie)
null<-glm(betas.LCBD~1,family=gaussian(),data=dd_specie)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod5, null,weights = TRUE, sort = FALSE)

mod1<-betareg(betas.LCBD~(River.dist.lake),data=dd_specie)
mod2<-betareg(betas.LCBD~Head.river.dist,data=dd_specie)
mod5<-betareg(betas.LCBD~River.dist.lake*Head.river.dist,data=dd_specie)
null<-betareg(betas.LCBD~1,data=dd_specie)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod5, null,weights = TRUE, sort = FALSE)
summary(mod1)
anova(mod1,mod2,mod5,mod6, null)

pseudoR1 <- ((mod1$null.deviance-mod1$deviance)/mod1$null.deviance)
pseudoR2 <- ((mod2$null.deviance-mod2$deviance)/mod2$null.deviance)
pseudoR0 <- ((mod5$null.deviance-mod5$deviance)/mod5$null.deviance)
pseudoRnull <- ((null$null.deviance-null$deviance)/null$null.deviance)

mod1<-betareg(betas.LCBD~(River.dist.lake),link = "logit",data=dd_specie)
mod2<-betareg(betas.LCBD~(Head.river.dist),link = "logit",data=dd_specie)
mod5<-betareg(betas.LCBD~River.dist.lake*Head.river.dist,link = "logit",data=dd_specie)
null<-betareg(betas.LCBD~1,link = "logit",data=dd_specie)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod5, null,weights = TRUE, sort = FALSE)


#2)Mixed Models
mod1<-lmer(betas.LCBD~River.dist.lake+ (1|O.NET), data=dd_specie)
mod2<-lmer(betas.LCBD~Head.river.dist+ (1|O.NET),data=dd_specie)
mod3<-lmer(betas.LCBD~Head.river.dist*River.dist.lake+ (1|O.NET),data=dd_specie)
null<-lmer(betas.LCBD~1+ (1|O.NET),data=dd_specie)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod3,null,weights = TRUE, sort = FALSE)
r2(mod3)
rsquared(mod3) 

mod1<-glmmTMB(betas.LCBD~River.dist.lake+ (1|O.NET),family=beta_family(), data=dd_specie)
mod2<-glmmTMB(betas.LCBD~Head.river.dist+ (1|O.NET),family=beta_family(),data=dd_specie)
mod3<-glmmTMB(betas.LCBD~Head.river.dist*River.dist.lake+ (1|O.NET),family=beta_family(),data=dd_specie)
null<-glmmTMB(betas.LCBD~1+ (1|O.NET),family=beta_family(),data=dd_specie)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod3,null,weights = TRUE, sort = FALSE)

performance::r2(mod1)
performance::r2(mod2)
performance::r2(mod3)
performance::r2(null)
summary(mod3)
performance::r2_nakagawa(null,tolerance = 0)
hist(dd_specie$betas.LCBD)


mod1<-glmmTMB(N1~River.dist.lake+ (1|O.NET),family=gaussian(), data=dd_specie)
mod2<-glmmTMB(N1~Head.river.dist+ (1|O.NET),family=gaussian(),data=dd_specie)
mod3<-glmmTMB(N1~River.dist.lake*Head.river.dist+ (1|O.NET),family=gaussian(),data=dd_specie)
null<-glmmTMB(N1~1+ (1|O.NET),family=gaussian(),data=dd_specie)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod3, null,weights = TRUE, sort = F)
summary(mod6)

mod3<-glmer(N1~River.dist.lake*Head.river.dist+ (1|O.NET),family=gaussian(),data=dd_specie)
remotes::install_github("mastoffel/partR2", build_vignettes = TRUE, dependencies = TRUE) 
library(partR2)
partR2(mod3, R2_type = "conditional")
partR2(mod3, R2_type = "marginal")
r.squaredGLMM(mod3)
cor(dd_specie$River.dist.lake,dd_specie$Head.river.dist)
ggplot(dd_specie, aes(x=River.dist.lake,y=Head.river.dist))+geom_point()+geom_smooth(method="lm")
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
