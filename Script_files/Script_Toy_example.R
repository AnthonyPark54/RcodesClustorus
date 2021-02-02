##### Please set the working directory where you saved the provided files
# setwd("./R_Codes_for_submission")

library(MASS)
library(tidyverse)
library(BAMBI)
source('./Source_files/routines.R')

set.seed(20201) 
data <- rbind(rvmsin(n = 100, 50, 3,0,pi,pi),
              rvmsin(n = 100, 3  ,50 ,0,0,pi),
              rvmsin(n = 100, 50, 50,40,-.7,pi-0.7),
              c( 5/4*pi+0.1,pi),
              cbind( runif(10, 0, 2*pi), runif(10, 0, 2*pi) )  )

              

label <- c(rep(1,100), rep(2,100), rep(3,100),rep(4,11))
data <- data-pi
data <- on.torus(data)
data <- cbind(as.data.frame(data) , as.factor(label))
colnames(data) <- c("phi","psi","label")


# plot(dat2, xlim=c(0, 2*pi), ylim=c(0, 2*pi), xlab="phi", ylab="psi")

# Figure ------------------------------------------------------------------
data %>% ggplot(aes(x = phi, y = psi, color = label)) + geom_point() + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle('Data set 1 with true labels') 

icp.torus <- icp.torus.score(as.matrix(data[,1:2]), split.id = NULL,
                method = "mixture",
                mixturefitmethod = "g",
                param = list(J = 3))
ia <- icp.torus.eval(icp.torus, level = 0.1, eval.point = grid.torus())

b <- data.frame(ia$phi,ia$psi, ia$Chat_e == 1)
colnames(b) <- c("phi","psi","C_e")

#----------------------------------------------------------------------------------------
# Figure 4
#----------------------------------------------------------------------------------------
g0 <- ggplot() + geom_contour(aes(phi, psi, z = ifelse(C_e,1,0)), data = b, color = "grey", size = 1,lineend = "round" ) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2]), color = "darkgrey") + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi)))

g0 + annotate("text", label = "x", x = 1/4*pi+0.1, y = -0.15, size = 8, colour = "black") + 
  annotate("text", label = "E1", x = 0, y = 1/2*pi+0.1, size = 8, colour = "black") + 
  annotate("text", label = "E2", x = pi-0.3, y = 0.8, size = 8, colour = "black")+ 
  annotate("text", label = "E3", x = pi-0.3, y = - 1.2, size = 8, colour = "black")

ggsave("./Figures/Figure_4.png", width = 8, height = 8/4*3)

 