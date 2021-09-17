data<-read.csv("data.csv")

mod1<-glmmTMB(betas.LCBD~River.dist.lake+ (1|O.NET),family=beta_family(), data=data)
mod2<-glmmTMB(betas.LCBD~Head.river.dist+ (1|O.NET),family=beta_family(),data=data)
mod3<-glmmTMB(betas.LCBD~Head.river.dist*River.dist.lake+ (1|O.NET),family=beta_family(),data=data)
null<-glmmTMB(betas.LCBD~1+ (1|O.NET),family=beta_family(),data=data)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod3,null,weights = TRUE, sort = FALSE)

performance::r2(mod1)
performance::r2(mod2)
performance::r2(mod3)
performance::r2(null)
summary(mod3)
performance::r2_nakagawa(null,tolerance = 0)
hist(dd_specie$betas.LCBD)


mod1<-glmmTMB(N1~River.dist.lake+ (1|O.NET),family=gaussian(), data=data)
mod2<-glmmTMB(N1~Head.river.dist+ (1|O.NET),family=gaussian(),data=data)
mod3<-glmmTMB(N1~River.dist.lake*Head.river.dist+ (1|O.NET),family=gaussian(),data=data)
null<-glmmTMB(N1~1+ (1|O.NET),family=gaussian(),data=data)
reported.table2 <- bbmle::AICtab(mod1,mod2,mod3, null,weights = TRUE, sort = F)
summary(mod6)
