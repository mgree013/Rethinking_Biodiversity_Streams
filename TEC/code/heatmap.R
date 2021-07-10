#paper title: Rethinking Biodiversity in Stream Ecology Frameworks
#Plotting Two way  interactions with heatmaps and contour lines
#Author: Matthew Douglas Green
#Date: May 5,2020
#Reference: https://stats.stackexchange.com/questions/307863/level-plot-for-continuous-x-continuous-interaction-with-continuous-response

########################################################################################################################
library(cowplot)
library(viridis)
library(ggplot2)
library(reshape2)
library(metR)


lm.mod <- lm(mpg ~ wt*hp, data = mtcars)
summary(lm.mod)

prepplot <- as.data.frame(matrix(ncol = 3, nrow = 10000))
colnames(prepplot) <- c("hp", "wt", "est.mpg")

prepplot$hp <- rep(seq(52,335, length.out = 100), 100)
prepplot <- prepplot[order(prepplot$hp),]
prepplot$wt <- rep(seq(1.513,5.424, length.out = 100), 100)
prepplot$est.mpg <- 49.80842 - 8.21662*prepplot$wt - 0.12010*prepplot$hp + 
  0.02785*prepplot$wt*prepplot$hp
prepplot$est.mpg<-predict(lm.mod,prepplot, type="response")

ggplot(prepplot, aes(wt, hp, fill = est.mpg)) + 
  geom_tile() +
  ggtitle("A)") +
  xlab("Weight (1000 lbs.)") + ylab("Horsepower") +
  scale_fill_gradientn(colours = c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))
################################################################################################################################################################################################
#TEC FRamework
#Beta Diversity-Space-Env
lm.mod<-glm(betas.LCBD~Spatial*E_PC1,family=gaussian(link = "log"),data=all_big_dat)
summary(lm.mod)

prepplot <- as.data.frame(matrix(ncol = 3, nrow = 10000))
colnames(prepplot) <- c("E_PC1", "Spatial", "est.beta")

prepplot$Spatial <- rep(seq(-2.352044,7.043287, length.out = 100), 100)
prepplot <- prepplot[order(prepplot$Spatial),]
prepplot$E_PC1 <- rep(seq(-4.696132,3.494251, length.out = 100), 100)
prepplot$est.beta <- -4.721324  -0.086734*prepplot$Spatial -0.004319*prepplot$E_PC1 + 
  0.031497*prepplot$Spatial*prepplot$E_PC1
prepplot$est.beta<-predict(lm.mod,prepplot, type="response")

b1<-ggplot(prepplot, aes(E_PC1, Spatial, fill = est.beta)) + 
  geom_tile() +
  ggtitle("A)") +
  xlab("Env") + ylab("Spatial") +
  scale_fill_gradientn(colours = c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Est.BD")

 #Beta-Space*Com size
lm.mod<-glm(betas.LCBD~Spatial*Com.Size.Gradient,family=gaussian(link = "log"),data=all_big_dat)
summary(lm.mod)

prepplot111 <- as.data.frame(matrix(ncol = 3, nrow = 10000))
colnames(prepplot111) <- c("Com.Size.Gradient", "Spatial", "est.beta")

prepplot111$Spatial <- rep(seq(-2.352044,7.043287, length.out = 100), 100)
prepplot111 <- prepplot111[order(prepplot111$Spatial),]
prepplot111$Com.Size.Gradient <- rep(seq(3.170086,9.273556, length.out = 100), 100)
prepplot111$est.beta <- -4.40098  + 0.04790*prepplot111$Spatial -0.04470*prepplot111$Com.Size.Gradient + 
  -0.02189*prepplot111$Spatial*prepplot111$Com.Size.Gradient
prepplot111$est.beta<-predict(lm.mod,prepplot111, type="response")

b2<-ggplot(prepplot111, aes(Com.Size.Gradient, Spatial, fill = est.beta)) + 
  geom_tile() +
  ggtitle("B)") +
  xlab("Com.Size.Gradient") + ylab("Spatial") +
  scale_fill_gradientn(colours = c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Est.BD")

#Beta-Env*Com size
lm.mod<-glm(betas.LCBD~E_PC1*Com.Size.Gradient,family=gaussian(link = "log"),data=all_big_dat)
summary(lm.mod)

prepplot112 <- as.data.frame(matrix(ncol = 3, nrow = 10000))
colnames(prepplot112) <- c("Com.Size.Gradient", "E_PC1", "est.beta")

prepplot112$E_PC1 <- rep(seq(-4.696132,3.494251, length.out = 100), 100)
prepplot112 <- prepplot112[order(prepplot112$E_PC1),]
prepplot112$Com.Size.Gradient <- rep(seq(3.170086,9.273556, length.out = 100), 100)
prepplot112$est.beta <- -4.444071   -0.083698*prepplot112$E_PC1 -0.037954*prepplot112$Com.Size.Gradient + 
  0.005456*prepplot112$E_PC1*prepplot112$Com.Size.Gradient
prepplot112$est.beta<-predict(lm.mod,prepplot112, type="response")

b3<-ggplot(prepplot112, aes(Com.Size.Gradient, E_PC1, fill = est.beta)) + 
  geom_tile() +
  ggtitle("C)") +
  xlab("Com.Size.Gradient") + ylab("Env") +
  scale_fill_gradientn(colours = c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Est.BD")

##################################################
#Local Diversity-Space-Env
lm.mod<-glm(N1~E_PC1*Com.Size.Gradient,family=gaussian(link = "identity"),data=all_big_dat)
summary(lm.mod)

prepplot0 <- as.data.frame(matrix(ncol = 3, nrow = 10000))
colnames(prepplot0) <- c("E_PC1", "Com.Size.Gradient", "est.n1")

str(all_big_dat$E_PC1)
min(all_big_dat$Com.Size.Gradient)

prepplot0$Com.Size.Gradient <- rep(seq(3.170086,9.273556, length.out = 100), 100)
prepplot0 <- prepplot0[order(prepplot0$Com.Size.Gradient),]
prepplot0$E_PC1 <- rep(seq(-4.696132,3.494251, length.out = 100), 100)
prepplot0$est.n1 <-  8.6419   -0.5536*prepplot0$Com.Size.Gradient +3.4244*prepplot0$E_PC1 + 
  -0.6086*prepplot0$Com.Size.Gradient*prepplot0$E_PC1
prepplot0$est.n1<-predict(lm.mod,prepplot0, type="response")


l1<-ggplot(prepplot0, aes(E_PC1, Com.Size.Gradient, fill = est.n1)) + 
  geom_tile() +
  ggtitle("D)") +
  xlab("Env") + ylab("Com.Size.Gradient") +
  scale_fill_gradientn(colours = c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Est.LD")

#Local diversity-Space*Com Size

lm.mod<-glm(N1~Spatial*Com.Size.Gradient,family=gaussian(link = "identity"),data=all_big_dat)
summary(lm.mod)

prepplot01 <- as.data.frame(matrix(ncol = 3, nrow = 10000))
colnames(prepplot01) <- c("Spatial", "Com.Size.Gradient", "SD")

max(all_big_dat$Spatial)
max(all_big_dat$Com.Size.Gradient)

prepplot01$Com.Size.Gradient <- rep(seq(3.170086,9.273556, length.out = 100), 100)
prepplot01 <- prepplot01[order(prepplot01$Com.Size.Gradient),]
prepplot01$Spatial <- rep(seq(-2.352044,7.043287, length.out = 100), 100)
prepplot01$SD <-  8.8852   -0.6209*prepplot01$Com.Size.Gradient +1.7605 *prepplot01$Spatial + 
  -0.2393*prepplot01$Com.Size.Gradient*prepplot01$Spatial
prepplot01$SD<-predict(lm.mod,prepplot01, type="response")


l2<-ggplot(prepplot01, aes(Spatial, Com.Size.Gradient, fill = est.n1)) + 
  geom_tile() +
  ggtitle("E)") +
  xlab("Spatial") + ylab("Com.Size.Gradient") +
  scale_fill_gradientn(colours = c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Est.LD")

#Local diversity-Space*Env

lm.mod<-glm(N1~Spatial*E_PC1,family=gaussian(link = "identity"),data=all_big_dat)
summary(lm.mod)

prepplot016 <- as.data.frame(matrix(ncol = 3, nrow = 10000))
colnames(prepplot016) <- c("Spatial", "E_PC1", "est.n1")



prepplot016$E_PC1 <- rep(seq(-4.696132,3.494251, length.out = 100), 100)
prepplot016 <- prepplot016[order(prepplot016$E_PC1),]
prepplot016$Spatial <- rep(seq(-2.352044,7.043287, length.out = 100), 100)
prepplot016$est.n1 <-  5.11753   -0.23928 *prepplot016$E_PC1 +0.59453  *prepplot016$Spatial + 
  0.03393 *prepplot016$E_PC1*prepplot016$Spatial
prepplot016$est.n1<-predict(lm.mod,prepplot016, type="response")


l3<-ggplot(prepplot016, aes(Spatial, E_PC1, fill = est.n1)) + 
  geom_tile() +
  ggtitle("F)") +
  xlab("Spatial") + ylab("Env") +
  scale_fill_gradientn(colours = c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Est.LD")

plot_grid(b1,b2,b3,l1,l2,l3,nrow=2)
################################################################################################################################################################################################
#Part 1) Stream Eco Frames

#1a)Beta
lm.mod <- lm(betas.LCBD ~ Head.river.dist*(River.dist.lake), data = dd_specie)
summary(lm.mod)

range(specie$Head.river.dist)
lake<-specie%>%drop_na()
range(lake$River.dist.lake)

prepplot1 <- as.data.frame(matrix(ncol = 3, nrow = 10000))
colnames(prepplot1) <- c("River.dist.lake", "Head.river.dist", "est.beta")

prepplot1$Head.river.dist <- rep(seq(5.332341 ,10.063097, length.out = 100), 100)
prepplot1 <- prepplot1[order(prepplot1$Head.river.dist),]
prepplot1$River.dist.lake <- rep(seq(0.6931472,9.6235643, length.out = 100), 100)
prepplot1$est.beta <- -4.99626     +0.21383 *prepplot1$Head.river.dist +0.75407 *prepplot1$River.dist.lake + 
  -0.09176 *prepplot1$Head.river.dist*prepplot1$River.dist.lake
prepplot1$est.beta<-predict(lm.mod,prepplot1, type="response")

aab<-ggplot(prepplot1, aes(River.dist.lake, Head.river.dist, fill = est.beta)) + 
  geom_tile() +
  ggtitle("E)") +
  xlab("Log Distance from Upstream Lakes (m)") + ylab("Log Distance from Headwaters (m)") +
  scale_fill_gradientn(colours = c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Est.BD")

aab+geom_point(data=dd_specie,size= betas.LCBD)
prepplot1_1<-prepplot1%>%add_column(dd_specie$betas.LCBD)%>%complete()
prepplot1_1$beta<-dd_specie$betas.LCBD

#Local
lm.mod <- lm(N1 ~ Head.river.dist*River.dist.lake, data = dd_specie)
summary(lm.mod)

range(dd_specie$Head.river.dist)
lake<-specie%>%drop_na()
range(dd_specie$River.dist.lake)

prepplot2 <- as.data.frame(matrix(ncol = 3, nrow = 10000))
colnames(prepplot2) <- c("River.dist.lake", "Head.river.dist", "est.n1")

prepplot2$Head.river.dist <- rep(seq(5.332341 ,10.063097, length.out = 100), 100)
prepplot2 <- prepplot2[order(prepplot2$Head.river.dist),]
prepplot2$River.dist.lake <- rep(seq(0.6931472,9.6235643, length.out = 100), 100)
prepplot2$est.n1 <- 10.2277    -0.9769 *prepplot2$Head.river.dist   -1.8474  *prepplot2$River.dist.lake + 
  0.2874*prepplot2$Head.river.dist*prepplot2$River.dist.lake
prepplot2$est.n1<-predict(lm.mod,prepplot2, type="response")

aabb<-ggplot(prepplot2, aes(River.dist.lake, Head.river.dist, fill = est.n1)) + 
  geom_tile() +
  ggtitle("F)") +
  xlab("Log Distance from Upstream Lakes (m)") + ylab("Log Distance from Headwaters (m)") +
  scale_fill_gradientn(colours = c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Est.LD")

#Figure 3
plot_grid(d.b1,d.e1,d.b2,d.e2,aab,aabb,nrow=3)

##############################################################################################################
#VIRIDIS
st.n1<-ggplot(prepplot2) +
  geom_tile(aes(River.dist.lake, Head.river.dist, fill = est.n1))+
  ggtitle("f)") +
  scale_fill_viridis_c()+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  ggplot2::xlab("Log Distance from Upstream Lakes (m)") + ylab("Log Distance from Headwaters (m)") +
  ggplot2::labs(fill = "SD")

st.bd<-ggplot(prepplot1) +
  geom_tile(aes(River.dist.lake, Head.river.dist, fill = est.beta))+
  ggtitle("e)") +
  scale_fill_viridis_c()+
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  ggplot2::xlab("Log Distance from Upstream Lakes (m)") + ylab("Log Distance from Headwaters (m)") +
  ggplot2::labs(fill = "BD")

################

tec.n1.1<-ggplot(prepplot0) +
  geom_tile(aes(E_PC1, Com.Size.Gradient, fill = est.n1))+
  scale_fill_viridis_c()+
  ggtitle("b)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  ggplot2::xlab("Environment") + ylab("Log Community Size") +
  ggplot2::labs(fill = "Est.SD")


tec.n1.2<-ggplot(prepplot01) +
  geom_tile(aes(Spatial, Com.Size.Gradient, fill = est.n1))+
  scale_fill_viridis_c()+
  ggtitle("d)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  ggplot2::xlab("Spatial") + ylab("Log Community Size") +
  ggplot2::labs(fill = "Est.SD")

tec.n1.3<-ggplot(prepplot016) +
  geom_tile(aes(Spatial, E_PC1, fill = est.n1))+
  scale_fill_viridis_c()+
  ggtitle("f)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  ggplot2::xlab("Spatial") + ylab("Environment") +
  ggplot2::labs(fill = "Est.SD")

  
tec.bd.3<-ggplot(prepplot) +
  geom_tile(aes( Spatial,E_PC1, fill = est.beta))+
  scale_fill_viridis_c()+
  ggtitle("e)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  ggplot2::ylab("Environment") + xlab("Spatial") +
  ggplot2::labs(fill = "Est.BD")
  
tec.bd.2<-ggplot(prepplot111) +
  geom_tile(aes( Spatial,Com.Size.Gradient, fill = est.beta))+
  scale_fill_viridis_c()+
  ggtitle("c)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  ggplot2::ylab("Log Community Size") + xlab("Spatial") +
  ggplot2::labs(fill = "Est.BD")

tec.bd.1<-ggplot(prepplot112) +
  geom_tile(aes( E_PC1,Com.Size.Gradient, fill = est.beta))+
  scale_fill_viridis_c()+
  ggtitle("a)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  ggplot2::ylab("Log Community Size") + xlab("Environment") +
  ggplot2::labs(fill = "Est.BD")



plot_grid(st.n1,st.bd)
plot_grid(tec.bd.1,tec.n1.1,tec.bd.2,tec.n1.2,tec.bd.3,tec.n1.3,nrow=3)

plot_grid(d.b1,d.e1,d.b2,d.e2,st.bd,st.n1,nrow=3)

plot_grid(d.b1,d.e1,d.b2,d.e2,st_bd_,st_n1_,nrow=3)

plot_grid(tec_bd1_,tec_n1_,tec_bd2_,tec_n2_,tec_bd3_,tec_n3_,nrow=3)


#############################################################################################################################
write.csv(prepplot112,"okay.csv")
write.csv(all_big_dat,"okay1.csv")

oh<-read.csv("okay.csv")

ggplot(data=oh,aes(x = E_PC1, y = Com.Size.Gradient,z=est.beta))+ 
  geom_contour(aes(colour = after_stat(level)))+
  #geom_text_contour(aes(z = est.beta), skip=2, colour = "black", nudge_x =0.2 ) +
  labs(x = "Env", 
       y = "Com.Size.Gradient", 
       z = "LCBD",level="LCBD") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_color_viridis()+
  geom_point(aes(x = ENV, y = Com.Size))

################################################################################################
#Stream Eco Frames Interactive Plots
#write.csv(dd_specie,"stream_data.csv")
library(viridis)
#write.csv(prepplot2,"st.n1.csv")
st_n1<-read.csv("st.n1.csv")

st_n1_<-ggplot(data=st_n1,aes(x = est.River.dist.lake, y = est.Head.river.dist,z=est.n1))+ 
  geom_point(aes(x = River.dist.lake, y = Head.river.dist, colour=SD))+
  geom_contour(aes(colour = after_stat(level)), size=1)+
  ggtitle("f)") +
  #geom_text_contour(aes(z = est.beta), skip=2, colour = "black", nudge_x =0.2 ) +
  labs(x = "Log Distance below Upstream Lakes (m)", 
       y = "Log Distance from Headwaters (m)", 
       z = "LCBD",level="LCBD") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_color_viridis()



#write.csv(prepplot1,"st.bd.csv")
st_bd<-read.csv("st.bd.csv")

st_bd_<-ggplot(data=st_bd,aes(x = est.River.dist.lake, y = est.Head.river.dist,z=est.beta))+ 
  geom_point(aes(x = River.dist.lake, y = Head.river.dist, colour=BD))+
  geom_contour(aes(colour = after_stat(level)), size=1)+
  ggtitle("e)") +
  #geom_text_contour(aes(z = est.beta), skip=2, colour = "black", nudge_x =0.2 ) +
  labs(x = "Log Distance below Upstream Lakes (m)", 
       y = "Log Distance from Headwaters (m)", 
       z = "LCBD",level="LCBD") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_color_viridis()
################################################################################################
#TEC Frames Interactive Plots
#write.csv(all_big_dat,"okay1.csv")

#write.csv(prepplot0,"tec.n1.csv")
tec_n1<-read.csv("tec.n1.csv")

tec_n1_<-ggplot(data=tec_n1,aes(x = est.E_PC1, y = est.Com.Size.Gradient,z=est.n1))+ 
  geom_point(aes(x = E_PC1, y = Com.Size.Gradient, colour=SD))+
  geom_contour(aes(colour = after_stat(level)))+
  ggtitle("b)") +
  #geom_text_contour(aes(z = est.beta), skip=2, colour = "black", nudge_x =0.2 ) +
  labs(x = "Environment", 
       y = "Log Community Size", 
       z = "LCBD",level="LCBD") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_color_viridis()

write.csv(prepplot01,"tec.n2.csv")
tec_n2<-read.csv("tec.n2.csv")

tec_n2_<-ggplot(data=tec_n2,aes(x = Spatial, y = Com.Size.Gradient,z=SD))+ 
  geom_point(aes(x = Spatial, y = Com.Size.Gradient, colour=SD))+
  geom_contour(aes(colour = after_stat(level)))+
  ggtitle("d)") +
  #geom_text_contour(aes(z = est.beta), skip=2, colour = "black", nudge_x =0.2 ) +
  labs(x = "Spatial", 
       y = "Log Community Size", 
       z = "LCBD",level="LCBD") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_color_viridis()


#write.csv(prepplot016,"tec.n3.csv")
tec_n3<-read.csv("tec.n3.csv")

tec_n3_<-ggplot(data=tec_n3,aes(x = est.Spatial, y = est.E_PC1,z=est.n1))+ 
  geom_point(aes(x = Spatial, y = E_PC1, colour=SD))+
  #geom_contour(aes(colour = after_stat(level)))+
  ggtitle("f)") +
  #geom_text_contour(aes(z = est.beta), skip=2, colour = "black", nudge_x =0.2 ) +
  labs(x = "Spatial", 
       y = "Environment", 
       z = "LCBD",level="LCBD") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_color_viridis()





#write.csv(prepplot112,"tec.bd1.csv")
tec_bd1<-read.csv("tec.bd1.csv")

tec_bd1_<-ggplot(data=tec_bd1,aes( x = est.E_PC1,y = est.Com.Size.Gradient,z=est.beta))+ 
  geom_point(aes( x = E_PC1,y = Com.Size.Gradient, colour=BD))+
  geom_contour(aes(colour = after_stat(level)))+
  ggtitle("a)") +
  #geom_text_contour(aes(z = est.beta), skip=2, colour = "black", nudge_x =0.2 ) +
  labs(y = "Log Community Size", 
       x = "Environment", 
       z = "LCBD",level="LCBD") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_color_viridis()

#write.csv(prepplot111,"tec.bd2.csv")
tec_bd2<-read.csv("tec.bd2.csv")

tec_bd2_<-ggplot(data=tec_bd2,aes( x = est.Spatial,y = est.Com.Size.Gradient,z=est.beta))+ 
  geom_point(aes( x = Spatial, y = Com.Size.Gradient,colour=BD))+
  geom_contour(aes(colour = after_stat(level)))+
  ggtitle("c)") +
  #geom_text_contour(aes(z = est.beta), skip=2, colour = "black", nudge_x =0.2 ) +
  labs(y = "Log Community Size", 
       x = "Spatial", 
       z = "LCBD",level="LCBD") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_color_viridis()

#write.csv(prepplot,"tec.bd3.csv")
tec_bd3<-read.csv("tec.bd3.csv")

tec_bd3_<-ggplot(data=tec_bd3,aes(x = est.Spatial,y = est.E_PC1, z=est.beta))+ 
  geom_point(aes( x = Spatial,y = E_PC1, colour=BD))+
  geom_contour(aes(colour = after_stat(level)))+
  ggtitle("e)") +
  #geom_text_contour(aes(z = est.beta), skip=2, colour = "black", nudge_x =0.2 ) +
  labs(y = "Environment", 
       x = "Spatial", 
       z = "LCBD",level="LCBD") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_color_viridis()
##########################################################################################
all_big_dat%>%
  ggplot(aes(x=Spatial,y=log(Head.river.dist+1)))+
  geom_point()+
  geom_smooth(method="lm")

all_big_dat%>%
  ggplot(aes(x=Spatial,y=log(River.dist.lake+1)))+
  geom_point()+
  geom_smooth(method="lm")

all_big_dat%>%
  ggplot(aes(x=Com.Size.Gradient,y=log(Head.river.dist+1)))+
  geom_point()+
  geom_smooth(method="lm")

all_big_dat%>%
  ggplot(aes(x=Com.Size.Gradient,y=log(River.dist.lake+1)))+
  geom_point()+
  geom_smooth(method="lm")

all_big_dat%>%
  ggplot(aes(x=E_PC1,y=log(Head.river.dist+1)))+
  geom_point()+
  geom_smooth(method="lm")

all_big_dat%>%
  filter(log(River.dist.lake+1)>0)%>%
  ggplot(aes(x=E_PC1,y=log(River.dist.lake+1)))+
  geom_point()+
  geom_smooth(method="lm")

all_big_dat%>%
  filter(log(River.dist.lake+1)>0)%>%
  ggplot(aes(x = log(River.dist.lake+1), y = log(Head.river.dist+1),z=E_PC1))+ 
  geom_point(aes(x = log(River.dist.lake+1), y = log(Head.river.dist+1), colour=E_PC1))+
  geom_contour(aes(colour = after_stat(level)))+
  #geom_text_contour(aes(z = est.beta), skip=2, colour = "black", nudge_x =0.2 ) +
  labs(x = "Lake Dist", 
       y = "Head Dist", 
       z = "LCBD",level="LCBD") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_color_viridis()



#Local
all_big_datss<-all_big_dat%>%
  mutate(River.dist.lake=log(1+River.dist.lake),Head.river.dist=log(1+Head.river.dist))

lm.mod <- lm(E_PC1 ~ Head.river.dist*River.dist.lake, data = all_big_datss)
summary(lm.mod)

range(all_big_datss$Head.river.dist)
range(all_big_datss$River.dist.lake)

prepplot211 <- as.data.frame(matrix(ncol = 3, nrow = 10000))
colnames(prepplot211) <- c("River.dist.lake", "Head.river.dist", "est.Env")

prepplot211$Head.river.dist <- rep(seq(2.397895 ,10.063097, length.out = 100), 100)
prepplot211 <- prepplot211[order(prepplot211$Head.river.dist),]
prepplot211$River.dist.lake <- rep(seq(0.000000,9.6235643, length.out = 100), 100)
prepplot211$est.Env <- -6.42886      +0.90477  *prepplot211$Head.river.dist   -0.90876  *prepplot211$River.dist.lake + 
  0.09358*prepplot211$Head.river.dist*prepplot211$River.dist.lake
prepplot211$est.Env<-predict(lm.mod,prepplot211, type="response")

ggplot(prepplot211, aes(River.dist.lake, Head.river.dist, fill = est.Env)) + 
  geom_tile() +
  xlab("Log Distance from Upstream Lakes (m)") + ylab("Log Distance from Headwaters (m)") +
  scale_fill_gradientn(colours = c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Env")

