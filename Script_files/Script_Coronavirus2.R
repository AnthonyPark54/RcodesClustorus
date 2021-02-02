library(tidyverse)
library(bio3d)
library(cowplot)
library(ClusterR)
library(mclust)  
library(ggdendro)
source('./Source_files/routines.R')

# bio data ------------------------------------------ 
  # pdb.6M15 <- torsion.pdb(read.pdb("6M15"))  # HKU2
  # pdb.6M16 <- torsion.pdb(read.pdb("6M16") ) # SADS-Cov
  # pdb.6VXX <- torsion.pdb(read.pdb("6VXX"))  # SARS-Cov2 (Covid-19)
  # save(pdb.6M15,pdb.6M16,pdb.6VXX,file = "coronavirusdata.RData") 
load(file = "./Corona_Data/coronavirusdata.RData")

dataXX <- torsion_to_dat(pdb.6VXX)
dataXX <- dataXX %>% filter(type %in% c("B"))

dataXX$phi <- dataXX$phi/180*pi
dataXX$psi <- dataXX$psi/180*pi
n <- nrow(dataXX)

data <- cbind(dataXX$phi, dataXX$psi)
split.id <- rep(c(1,2),n/2)

# Clustering by existing methods ------------------------------------------
J <- 3
filename <- "./Figures/"
# source("Ex_Master.R") -------------------

# Perform naive clustering
## By naive K-means

kmeans.out<-KMeans_rcpp(data,clusters = J)
g_kmeans <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(kmeans.out$clusters)) %>% 
  ggplot(aes(phi,psi, shape = membership, color = membership)) + geom_point()  + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle("Naive K-means")


## By Non-angular Gaussian Mixture (using state-of-the-art mclust)
# define angular distance: 

pdist.data2 <- ang.pdist(data) # Use the pairwise L2-angular distance for clustering 

## PAM (Partitioning around medoids - Kaufman and Rousseeuw(1990) )
pam.out <- Cluster_Medoids(as.matrix(pdist.data2),clusters = J)
g_pam <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(pam.out$clusters)) %>% 
  ggplot(aes(phi,psi, color = membership)) + geom_point() + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle("K-means (ignoring the angular constraint)") + ggtitle("Partitioning around medoids")

## K-means in the ambient space with kmeans++ initialization

kmeans.out<-KMeans_rcpp(cbind(cos(data),sin(data)),clusters = J)
g_kmeans2 <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(kmeans.out$clusters)) %>% 
  ggplot(aes(phi,psi, color = membership)) + geom_point() + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle("K-means with chordal distance (in the ambient space) ")

## Hierarchical
hc.complete <- hclust(pdist.data2, method="complete")
# ggdendrogram(hc.complete, rotate=TRUE, size=2) + labs(title="Complete Linkage")
#ggdendrogram(hc.complete,main="Average Linkage", xlab="", sub="", cex=.9)
membership <- cutree(hc.complete, J)
g_hier <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(membership)) %>% 
  ggplot(aes(phi,psi, color = membership)) + geom_point() + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle("Hierachical clustering with average L2-Angular distance")


#----------------------------------------------------------------------------------------
# Figure C.2 - top 2 rows
#----------------------------------------------------------------------------------------

plot_grid(g_kmeans, g_pam, g_kmeans2 + ggtitle("Extrinsic K-means"),
          g_hier + ggtitle("Hierachical clustering"), label_size = 12)
ggsave(paste(filename,"Figure_C2_top.png",sep = ""), width = 9, height = 4*9/6)


# Find alpha and J --------------------------------------------------------
Jvec <- 3:65
l <- list()
set.seed(0)
n <- nrow(data)
n2 <- sum(split.id == 2)
for (j in Jvec){
 l[[j]] <- icp.torus.score(data, split.id = split.id,
                           method = "mixture",
                           mixturefitmethod = "a",
                           param = list(J = j))
}
alphavec <- 1:floor(n2/2) / n2
N <- length(alphavec)

# need a data frame (alpha, J, mu, alpha + mu)
out <- data.frame()
for (j in Jvec){
 Mvec <- alphavec
 a<-icp.torus.eval(l[[j]], level = alphavec, eval.point = grid.torus())
 for (i in 1:N){
   Mvec[i] <- sum(a$Chat_e[,i])/10000
 }
 out <- rbind(out, data.frame(alpha = alphavec, J = j, mu = Mvec, criterion = alphavec + Mvec))
}

out %>% ggplot(aes(x= alpha, y = criterion, color = J)) + geom_point()
# ggsave(paste(filename,"7icp_areaC.png",sep = ""), width = 6, height = 4)


out.index <- which.min(out$criterion)
out[out.index,] 

Jhat <- out[out.index,2]
alphahat <- out[out.index,1]

out %>% filter(alpha == alphahat) %>% ggplot(aes(J,mu)) + geom_point() + geom_line()

out %>% ggplot(aes(x= alpha, y = criterion, color = J)) + geom_point() +
  theme_bw() +  xlab(expression(~ alpha))

# ICP-Clustering by chosen alpha and J ------------------------------------
icp.torus <- l[[Jhat]]

ia <- icp.torus.eval(icp.torus, level = alphahat, eval.point = grid.torus())
b <- data.frame(ia$phi,ia$psi, ia$Chat_mix == 1, ia$Chat_max == 1, ia$Chat_e == 1)
colnames(b) <- c("phi","psi","C_mix","C_max","C_e")
head(b)

b<- b %>%  pivot_longer(3:5, names_to = "Type", values_to = "value")
g0 <- ggplot() + geom_contour(aes(phi, psi, z = ifelse(value,1,0),linetype = Type, color =Type), data = b, size = 1,lineend = "round" ) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(paste("ICP, vM2 mixture with J=",Jhat,", alpha=",alphahat))

c <- cluster.assign.torus(data, icp.torus, level = alphahat) 
g_e <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(c$cluster.id.by.ehat)) %>% 
  ggplot(aes(phi,psi, color = membership)) + geom_point() + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(paste("Clusters, K=", c$ncluster))


membership <- c$cluster.id.outlier 
membership <- ifelse(membership == max(unique(c$cluster.id.outlier)),"out",membership)

g1 <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(membership)) %>%  
  ggplot(aes(phi,psi, color = membership)) + geom_point() + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(paste("Clusters and outliers, K=",length(unique(c$cluster.id.outlier))))


g_mah <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(c$cluster.id.by.Mah.dist)) %>% 
  ggplot(aes(phi,psi,shape = membership, color = membership)) + geom_point() + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(paste("Clustering by Mahalanobis distance, K=",c$ncluster))


g2 <-  data.frame(phi = data[,1], psi = data[,2], membership = as.factor(c$cluster.id.by.ehat)) %>% 
  ggplot(aes(phi,psi, shape = membership, color = membership)) + geom_point() + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(paste("SARS-CoV-2, K=", c$ncluster))
g2_grey <- g2 + scale_color_grey()

level = alphahat
n2 <- icp.torus$n2
ialpha <- floor( (n2 + 1) * level)
t <- icp.torus$mixture$score_ellipse[ialpha]
# Draw.ellipses.bvn.approx(data, icp.torus$mixture$fit$parammat, t, data, c$cluster.id.outlier)

ellipse.param <- icp.torus$mixture$ellipsefit
J <- length(ellipse.param$mu1)

# all_nine_ellipses
theta <- seq(-pi, pi,length.out = 999)
Z <- cbind(cos(theta), sin(theta))

shift <- matrix(0,ncol = 2, nrow = 9)
shift[,1] <- c(0,2*pi,-2*pi)
shift[,2] <- rep(c(0,2*pi,-2*pi), each = 3)


for(j in 1:J){
  mu <- c(ellipse.param$mu1[j], ellipse.param$mu2[j])
  Sinv <- ellipse.param$Sigmainv[[j]]
  c.minus.t <- ellipse.param$c[j] - t
  
  if(c.minus.t < 0){
    cat("skip",j,",")
    next}
  cat("draw",j,",")
  M <- eigen(Sinv/c.minus.t)
  Mmhalf <- M$vectors %*% diag( sqrt(1/M$values) ) %*% t(M$vectors)
  R <- Mmhalf %*% t(Z) 
  for( shift.id in 1:9){
    RR <- R + mu + shift[shift.id,]
    g2 <-   g2 + geom_polygon(aes(x = phi, y = psi, shape = as.factor(1), color =as.factor(1)),
                              color = "black", alpha = 0.1, data = data.frame(phi = RR[1,],psi = RR[2,], value = 1))
    g2_grey <-   g2_grey + geom_polygon(aes(x = phi, y = psi, shape = as.factor(1), color =as.factor(1)),
                                        color = "black", alpha = 0.1, data = data.frame(phi = RR[1,],psi = RR[2,], value = 1))
  }
  
}

#----------------------------------------------------------------------------------------
# Figure 8 - bottom
#----------------------------------------------------------------------------------------
g2 <- g2 + coord_cartesian(xlim = c(-pi, pi), ylim = c(-pi,pi))
g2_grey <- g2_grey + coord_cartesian(xlim = c(-pi, pi), ylim = c(-pi,pi))
g2
ggsave(paste(filename,"Figure_8_bottom.png",sep = ""), width = 6, height = 4)
g2_grey
ggsave(paste(filename,"Figure_8_bottom_grey.png",sep = ""), width = 6, height = 4)



cp.torus.obj<-cp.torus.kde(data, level = 0.1)

b<- cp.torus.obj$cp.torus %>% pivot_longer(3:5, names_to = "Type", values_to = "value")
g_cp <- ggplot() + geom_contour(aes(phi, psi, z = ifelse(value,1,0),linetype = Type, color =Type), data = b, size = 1,lineend = "round" ) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi)))


#----------------------------------------------------------------------------------------
# Figure C.2 - bottom 2 rows
#----------------------------------------------------------------------------------------
plot_grid(g_cp + ggtitle("KDE-based Prediction sets"), g0 + ggtitle("Mixture-based Prediction sets"),
          g_e, g1, label_size = 12)

ggsave(paste(filename,"Figure_C2_bottom.png",sep = ""), width = 9, height = 6)


####################
out[out.index,]
Jhat
alphahat

chooseJ <- Jhat 
icp.torus <- l[[chooseJ]]
p1 <- plot_clustering(icp.torus, choosealpha = 0.02)
p2 <- plot_clustering(icp.torus, choosealpha = 0.03)
p3 <- plot_clustering(icp.torus, choosealpha = 0.05)
p4 <- plot_clustering(icp.torus, choosealpha = 0.10)

plot_grid(p1, p2, p3,p4, label_size = 12)
#ggsave(paste(filename,"ClusterPlot.png",sep=""), width = 12, height = 8)

# kappa and alpha 
# Now with kde-based fit

set.seed(101)
n <- nrow(data)

icp.torus.kde<- icp.torus.score(data, split.id = split.id,
                                method = "kde",
                                mixturefitmethod = "a",
                                param = list(concentration = 25))
icp.kde.region <- icp.torus.eval(icp.torus.kde, level = 0.02, eval.point = grid.torus())
p1<- data.frame(Chat = as.vector(icp.kde.region$Chat_kde), psi = icp.kde.region$psi, phi = icp.kde.region$phi) %>% 
  ggplot() + geom_contour(aes(phi,psi,z = Chat)) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(expression(paste("KDE:", ~kappa, " = 25," ~alpha,  " = 0.02")))

icp.kde.region <- icp.torus.eval(icp.torus.kde, level = 0.1, eval.point = grid.torus())
p2<- data.frame(Chat = as.vector(icp.kde.region$Chat_kde), psi = icp.kde.region$psi, phi = icp.kde.region$phi) %>% 
  ggplot() + geom_contour(aes(phi,psi,z = Chat)) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(expression(paste("KDE:", ~kappa, " = 25," ~alpha,  " = 0.1")))


icp.torus.kde<- icp.torus.score(data, split.id = split.id,
                                method = "kde",
                                mixturefitmethod = "a",
                                param = list(concentration = 100))
icp.kde.region <- icp.torus.eval(icp.torus.kde, level = 0.02, eval.point = grid.torus())
p3<- data.frame(Chat = as.vector(icp.kde.region$Chat_kde), psi = icp.kde.region$psi, phi = icp.kde.region$phi) %>% 
  ggplot() + geom_contour(aes(phi,psi,z = Chat)) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(expression(paste("KDE:", ~kappa, " = 100," ~alpha,  " = 0.02")))

icp.kde.region <- icp.torus.eval(icp.torus.kde, level = 0.1, eval.point = grid.torus())
p4<- data.frame(Chat = as.vector(icp.kde.region$Chat_kde), psi = icp.kde.region$psi, phi = icp.kde.region$phi) %>% 
  ggplot() + geom_contour(aes(phi,psi,z = Chat)) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2,
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(expression(paste("KDE:", ~kappa, " = 100," ~alpha,  " = 0.1")))

plot_grid(p1, p2, p3,p4, label_size = 12)


# choice of alpha and j browse options ------------------------------------
concentration.vec <- 10^( seq(log10(4),log10(150),length.out = 25) )
nkappa <- length(concentration.vec)
mu.vec1 <- concentration.vec
mu.vec2 <- concentration.vec

n <- nrow(data)

for (j in 1:nkappa){
  icp.torus.kde<- icp.torus.score(data,method = "kde",split.id = split.id ,param = list(concentration = concentration.vec[j]))
  icp.kde.region <- icp.torus.eval(icp.torus.kde, level = c(0.05,0.1), eval.point = grid.torus()) 
  mu.vec1[j] <- sum(icp.kde.region$Chat_kde[,1])/10000
  mu.vec2[j] <- sum(icp.kde.region$Chat_kde[,2])/10000
}

pICP_KDE_k <- data.frame(mu1 = mu.vec1, mu2 = mu.vec2, concentration = concentration.vec) %>% 
  pivot_longer(1:2, names_to = "Level") %>% mutate(Level = as.factor ( ifelse(Level=="mu1",0.05,0.1)) ,mu = value) %>% 
  ggplot(aes(concentration, mu, color = Level, group = Level)) + 
  geom_point() + geom_line()  + scale_x_log10() + theme_bw() +
  xlab(expression(~ kappa)) + ylab(expression(~ mu)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5))

# alpha vs j for MIXTURE --------------------------------------------------
out[out.index,]

pICP_MIX_J <- out %>% filter( abs(alpha - 0.1) < 0.001 | abs(alpha - 0.05) < 0.001 ) %>% 
  mutate(Level = as.factor(ifelse( abs(alpha - 0.1) < 0.001, 0.1, 
                                   0.05))) %>%  ggplot(aes(J,mu, color = Level)) +
  geom_point() + geom_line() + theme_bw() +
  ylab(expression(~ mu)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5))

#----------------------------------------------------------------------------------------
# Figure C.3
#----------------------------------------------------------------------------------------
plot_grid(pICP_KDE_k + ggtitle(expression(paste("KDE, over ", ~kappa))), 
          pICP_MIX_J+ ggtitle('Mixture, over J'), label_size = 12)
ggsave(paste(filename,"Figure_C3_top.png",sep=""), width = 6, height = 3)

out %>% group_by(alpha) %>% ggplot(aes(alpha,criterion,color = J)) + theme_bw() +
  geom_point() + scale_colour_distiller(palette="Spectral") + xlab(expression(~ alpha)) +
  geom_abline(slope = 0, intercept = out[out.index,4])
ggsave(paste(filename,"Figure_C3_bottom.png",sep=""), width = 6, height = 3)
 
