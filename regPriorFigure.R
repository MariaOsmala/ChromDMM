
#WHITOUT REGULARIZATION
shape<-3
rate<-1
x<-seq(from = 0, to = 5, by = 0.05) 
y<-seq(from = 0, to = 5, by = 0.05)
df<- expand.grid(x=x,y=y)

df$prob <- dgamma(df$x,shape=shape,rate=rate)*dgamma(df$y,shape=shape,rate=rate)
#df$prob <- df$prob/sum(df$prob)

library(ggplot2)
library(colorRamps)      # for matlab.like(...)
library(scales)          # for labels=scientific
library(latex2exp)
setwd("Desktop")

fontsize=12
ggplot(df, aes(x,y))+
  geom_raster(aes(fill=prob), interpolate = TRUE)+
  scale_fill_gradientn(colours=c("white","dark green"),)+scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_continuous(expand =c(0.0, 0.0))+
  labs(title=TeX("$p(\\alpha_1,\\alpha_{2})\\propto\\Gamma(\\alpha_1\\,|\\,\\eta =3,\\nu =1)\\cdot\\Gamma(\\alpha_2\\,|\\,\\eta =3,\\nu =1)$"))+ 
  labs(fill="Unnormalised\ndensity")+
  xlab(TeX("$\\alpha_1$"))+ ylab(TeX("$\\alpha_2$")) +
  theme(plot.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border=element_blank()) +
  theme(axis.line.x=element_line(color="black", size=0.5),
        axis.line.y=element_line(color="black", size=0.5),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        title=element_text(size=14),
        axis.text.x=element_text(size=fontsize),
        axis.text.y=element_text(size=fontsize),
        axis.title.x=element_text(size=fontsize),
        axis.title.y=element_text(size=fontsize))
 ggsave(file="NoRegularization.eps", device="eps")

##############################################################
#WITH REGULARIZATION 
 
 shape<-3
 rate<-1
 shape_h<-1
 rate_h<-4
 
 x<-seq(from = 0, to = 5, by = 0.05) 
 y<-seq(from = 0, to = 5, by = 0.05)
 df<- expand.grid(x=x,y=y)
 h=(df$x-df$y)^2
 df$prob <-dgamma(h, shape=shape_h, rate=rate_h)*dgamma(df$x,shape=shape,rate=rate)*dgamma(df$y,shape=shape,rate=rate)
 #df$prob <- df$prob/sum(df$prob)
 
 
 setwd("Desktop")
 
 fontsize=12
 ggplot(df, aes(x,y))+
   geom_raster(aes(fill=prob), interpolate = TRUE)+
   scale_fill_gradientn(colours=c("white","dark green"))+scale_x_continuous(expand = c(0.0, 0.0)) +
   scale_y_continuous(expand =c(0.0, 0.0))+
   labs(title=TeX("$p(\\alpha_1,\\alpha_{2})\\propto\\Gamma($h$=(\\alpha_1-\\alpha_2)^2 \\,|\\,\\eta_h =1, \\nu_h = 4)\\cdot\\Gamma(\\alpha_1\\,|\\,\\eta =3,\\nu =1)\\cdot\\Gamma(\\alpha_2\\,|\\,\\eta =3,\\nu =1)$"))+ 
   labs(fill="Unnormalised\ndensity")+
   xlab(TeX("$\\alpha_1$"))+ ylab(TeX("$\\alpha_2$")) +
   theme(plot.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border=element_blank()) +
   theme(axis.line.x=element_line(color="black", size=0.5),
         axis.line.y=element_line(color="black", size=0.5),
         legend.title=element_text(size=10),
         legend.text=element_text(size=10),
         title=element_text(size=14),
         axis.text.x=element_text(size=fontsize),
         axis.text.y=element_text(size=fontsize),
         axis.title.x=element_text(size=fontsize),
         axis.title.y=element_text(size=fontsize))
 ggsave(file="WithRegularization.eps")
 