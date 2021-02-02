##### Please set the working directory where you saved the provided files
# setwd("./R_Codes_for_submission")

library(tidyverse)
library(bio3d)
library(cowplot)
library(ClusterR)
library(mclust)  
library(ggdendro)
source('./Source_files/routines.R')

# bio data ------------------------------------------ 
   # pdb.6M15 <- torsion.pdb(read.pdb("6M15"))  # HKU2
   # pdb.6M16 <- torsion.pdb(read.pdb("6M16"))  # SADS-CoV
   # pdb.6VXX <- torsion.pdb(read.pdb("6VXX"))  # SARS-CoV-2 (Covid-19)
   # save(pdb.6M15,pdb.6M16,pdb.6VXX, file = "./Corona_Data/coronavirusdata.RData") 
load(file = "./Corona_Data/coronavirusdata.RData")

data16 <- torsion_to_dat(pdb.6M16) %>% filter(type == "B") %>% select(phi,psi, position)

#----------------------------------------------------------------------------------------
# Figure 1
#----------------------------------------------------------------------------------------
data16 %>% ggplot(aes(phi, psi)) + geom_point(size = 1) + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*90, limits = c(-180, 180),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*90, limits = c(-180, 180),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi)))
ggsave("./Figures/Figure_1.png", width = 6, height = 4)


dataXX <- torsion_to_dat(pdb.6VXX) %>% filter(type %in% c("B"))  %>% select(phi,psi, position)
data15 <- torsion_to_dat(pdb.6M15) %>% filter(type == "B") %>% select(phi,psi, position) 

alldata <- rbind(data.frame(data15,Source ="HKU2"),
                 data.frame(data16,Source ="SADS-CoV"),
                 data.frame(dataXX,Source ="SARS-CoV-2"))


#----------------------------------------------------------------------------------------
# Figure C.1
#----------------------------------------------------------------------------------------
alldata %>%   ggplot(aes(phi, psi,color = Source, shape = Source)) + geom_point(size = 1) + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*90, limits = c(-180, 180),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*90, limits = c(-180, 180),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi)))
ggsave("./Figures/Figure_C1.png", width = 6, height = 4)


data <- cbind(data16$phi/180*pi,data16$psi/180*pi)
data <- data[-which(is.na(data[,1])|is.na(data[,2])),]


##### Clustering by existing methods -------------------------------------------------------------------
J <- 3
filename <- "./Figures/"

### Perform naive clustering
## By naive K-means
kmeans.out <- KMeans_rcpp(data,clusters = J)
g_kmeans <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(kmeans.out$clusters)) %>% 
  ggplot(aes(phi,psi, shape = membership, color = membership)) + geom_point()  + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle("Naive K-means")
g_kmeans_grey <- g_kmeans + scale_color_grey()


## By Non-angular Gaussian Mixture (using state-of-the-art mclust)
# define angular distance:
pdist.data2 <- ang.pdist(data) # Use the pairwise L2-angular distance for clustering 

## PAM (Partitioning around medoids - Kaufman and Rousseeuw(1990) )
pam.out <- Cluster_Medoids(as.matrix(pdist.data2),clusters = J)
g_pam <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(pam.out$clusters)) %>% 
  ggplot(aes(phi,psi, shape = membership, color = membership)) + geom_point() + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle("Partitioning around medoids")
g_pam_grey <- g_pam + scale_color_grey()


## K-means in the ambient space with kmeans++ initialization
kmeans.out<-KMeans_rcpp(cbind(cos(data),sin(data)),clusters = J)
g_kmeans2 <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(kmeans.out$clusters)) %>% 
  ggplot(aes(phi,psi, shape = membership, color = membership)) + geom_point() + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle("K-means with chordal distance (in the ambient space) ")
g_kmeans2_grey <- g_kmeans2 + scale_color_grey()


## Hierarchical
hc.complete <- hclust(pdist.data2, method = "complete")
membership <- cutree(hc.complete, J)
g_hier <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(membership)) %>% 
  ggplot(aes(phi,psi, shape = membership, color = membership)) + geom_point() + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle("Hierachical clustering with average L2-Angular distance") 
g_hier_grey <- g_hier + scale_color_grey()


#----------------------------------------------------------------------------------------
# Figure 3
#----------------------------------------------------------------------------------------
plot_grid(g_kmeans, g_pam, g_kmeans2 + ggtitle("Extrinsic K-means"),
          g_hier + ggtitle("Hierachical clustering"), label_size = 12)
ggsave(paste(filename,"Figure_3.png",sep = ""), width = 9, height = 4*9/6)

plot_grid(g_kmeans_grey, g_pam_grey, g_kmeans2_grey + ggtitle("Extrinsic K-means"),
          g_hier_grey + ggtitle("Hierachical clustering"), label_size = 12)
ggsave(paste(filename,"Figure_3_grey.png",sep = ""), width = 9, height = 4*9/6)


##### Conformal Prediction with KDE -------------------------------------------------------

cp.torus.obj<-cp.torus.kde(data, level = 0.1)

b <- cp.torus.obj$cp.torus %>% pivot_longer(4, names_to = "Type", values_to = "value")
g_cp <- ggplot() + geom_contour(aes(phi, psi, z = ifelse(value,1,0), linetype = Type), data = b, size= 1, color = "darkgrey", lineend = "round" ) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2]), size = 0.5)  + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle("KDE-based Prediction sets")
g_cp


b2 <- cp.torus.obj$cp.torus %>% pivot_longer(3:5, names_to = "Type", values_to = "value")
g_cp_color <- ggplot() + geom_contour(aes(phi, psi, z = ifelse(value,1,0),linetype = Type, color =Type), data = b2, size = 0.1,lineend = "round" ) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) + 
  xlab(expression(~ phi)) + 
  ylab(expression(~ psi)) + 
  theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle("KDE-based Prediction sets")
g_cp_color

##### Find alpha and J ---------------------------------------------------------------------

Jvec <- 3:35
l <- list()
set.seed(0) 
n <- nrow(data)  
split.id <- NULL
for (j in Jvec){
  l[[j]] <- icp.torus.score(data, split.id = split.id,
                            method = "mixture",
                            mixturefitmethod = "a",
                            param = list(J = j))
}

n2 <- l[[10]]$n2
alphavec <- 1:floor(n2/2) / n2
N <- length(alphavec)

### need a data frame (alpha, J, mu, alpha + mu)
out <- data.frame()
for (j in Jvec){
  Mvec <- alphavec
  a<-icp.torus.eval(l[[j]], level = alphavec, eval.point = grid.torus())
  for (i in 1:N){
    Mvec[i] <- sum(a$Chat_e[,i])/10000
  }
  out <- rbind(out, data.frame(alpha = alphavec, J = j, mu = Mvec, criterion = alphavec + Mvec))
}

# # CANDIDATES FOR OTHER CRITERION. NOT USED.
# out <- out %>% mutate(criterion2 = alpha*(3/2) + mu, 
#                       criterion3 = sqrt(alpha^2 + mu^2))

out.index <- which.min(out$criterion)
out[out.index,] 

Jhat <- out[out.index,2]
alphahat <- out[out.index,1]

# ICP-Clustering by chosen alpha and J ------------------------------------
icp.torus <- l[[Jhat]]

ia <- icp.torus.eval(icp.torus, level = alphahat, eval.point = grid.torus())
b <- data.frame(ia$phi,ia$psi, ia$Chat_mix == 1, ia$Chat_max == 1, ia$Chat_e == 1)
colnames(b) <- c("phi","psi","C_mix","C_max","C_e")
b <- b %>%  pivot_longer(5, names_to = "Type", values_to = "value")
g0 <- ggplot() + geom_contour(aes(phi, psi, z = ifelse(value,1,0), linetype = Type), data = b, size= 1, color = "darkgrey", lineend = "round" ) + 
  geom_point(mapping = aes(x,y), size = 0.5, data = data.frame(x = data[,1],y =data[,2])) + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle("Mixture-based Prediction sets")
g0

b2 <- data.frame(ia$phi,ia$psi, ia$Chat_mix == 1, ia$Chat_max == 1, ia$Chat_e == 1) 
colnames(b2) <- c("phi","psi","C_mix","C_max","C_e")
b2 <- b2 %>%pivot_longer(3:5, names_to = "Type", values_to = "value") 

g0_color <- ggplot() + geom_contour(aes(phi, psi, z = ifelse(value,1,0),linetype = Type, color =Type), data = b2, size = 0.1,lineend = "round" ) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle("Mixture-based Prediction sets")
g0_color

c <- cluster.assign.torus(data, icp.torus, level = alphahat) 
g_e <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(c$cluster.id.by.ehat)) %>% 
  ggplot(aes(phi,psi, shape = membership, color = membership)) + geom_point()  + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(paste("Clusters, K=", c$ncluster))
g_e_grey <- g_e + scale_color_grey()


membership <- c$cluster.id.outlier
membership <- ifelse(membership == max(unique(c$cluster.id.outlier)),"out",membership)
g1 <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(membership)) %>%  
  ggplot(aes(phi,psi, shape = membership, color = membership)) + geom_point()  + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(paste("Clusters and outliers, K=",length(unique(c$cluster.id.outlier))))
g1_grey <- g1 + scale_color_grey()


g_mah <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(c$cluster.id.by.Mah.dist)) %>% 
  ggplot(aes(phi,psi, shape = membership, color = membership)) + geom_point() + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(paste("Clustering by Mahalanobis distance, K=",c$ncluster))
g_mah_grey <- g_mah + scale_color_grey()


g2 <- g_e + ggtitle(paste("SADS-CoV, K=", c$ncluster))
g2_grey <- g_e_grey + ggtitle(paste("SADS-CoV, K=", c$ncluster))

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
    g2 <- g2 + geom_polygon(aes(x = phi, y = psi, shape = as.factor(1), color =as.factor(1)),
                              color = "black", alpha = 0.1, data = data.frame(phi = RR[1,],psi = RR[2,], value = 1))
    g2_grey <- g2_grey + geom_polygon(aes(x = phi, y = psi, shape = as.factor(1), color =as.factor(1)),
                                        color = "black", alpha = 0.1, data = data.frame(phi = RR[1,],psi = RR[2,], value = 1))
  }
  
}


#----------------------------------------------------------------------------------------
# Figure 8 - top
#----------------------------------------------------------------------------------------
g2 <- g2 + coord_cartesian(xlim = c(-pi, pi), ylim = c(-pi,pi))
g2_grey <- g2_grey + coord_cartesian(xlim = c(-pi, pi), ylim = c(-pi,pi))
g2
ggsave(paste(filename, "Figure_8_top.png",sep = ""), width = 6, height = 4)
g2_grey
ggsave(paste(filename,"Figure_8_top_grey.png",sep = ""), width = 6, height = 4)


#----------------------------------------------------------------------------------------
# Figure 2
#----------------------------------------------------------------------------------------
plot_grid(g_cp, g0, g_e, g1, label_size = 12)
ggsave(paste(filename,"Figure_2.png",sep = ""), width = 9, height = 6)

plot_grid(g_cp, g0, g_e_grey, g1_grey, label_size = 12)
ggsave(paste(filename,"Figure_2_grey.png",sep = ""), width = 9, height = 6)


#----------------------------------------------------------------------------------------
# Supplement figure (All predictions version of Figure 2 top row)
#----------------------------------------------------------------------------------------
plot_grid(g_cp_color,g0_color,label_size = 12)
ggsave(paste(filename,"Figure_Cxx.png",sep = ""), width = 9, height = 4)

#######################
out[out.index,]
Jhat
alphahat

chooseJ <- Jhat 
icp.torus <- l[[chooseJ]]
p1 <- plot_clustering(icp.torus, choosealpha = 0.02)
p2 <- plot_clustering(icp.torus, choosealpha = 0.05)
p3 <- plot_clustering(icp.torus, choosealpha = 0.10)
p4 <- plot_clustering(icp.torus, choosealpha = 0.20)


#----------------------------------------------------------------------------------------
# Figure 6
#----------------------------------------------------------------------------------------
plot_grid(p1, p2, p3, p4, label_size = 12)
ggsave(paste(filename,"Figure_6.png",sep=""), width = 9, height = 6)


### kappa and alpha 
## Now with kde-based fit

set.seed(101)
n <- nrow(data)
split.id <- rep(2,n)
split.id[ sample(n,floor(n/2)) ] <- 1

icp.torus.kde<- icp.torus.score(data, split.id = split.id,
                            method = "kde",
                            mixturefitmethod = "a",
                            param = list(concentration = 25))
icp.kde.region <- icp.torus.eval(icp.torus.kde, level = 0.02, eval.point = grid.torus())
p1<- data.frame(Chat = as.vector(icp.kde.region$Chat_kde), psi = icp.kde.region$psi, phi = icp.kde.region$phi) %>% 
  ggplot() + geom_contour(aes(phi,psi,z = Chat), color = "darkgrey") + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(expression(paste("KDE:", ~kappa, " = 25," ~alpha,  " = 0.02")))

icp.kde.region <- icp.torus.eval(icp.torus.kde, level = 0.1, eval.point = grid.torus())
p2<- data.frame(Chat = as.vector(icp.kde.region$Chat_kde), psi = icp.kde.region$psi, phi = icp.kde.region$phi) %>% 
  ggplot() + geom_contour(aes(phi,psi,z = Chat), color = "darkgrey") + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(expression(paste("KDE:", ~kappa, " = 25," ~alpha,  " = 0.1")))


icp.torus.kde<- icp.torus.score(data, split.id = split.id,
                                method = "kde",
                                mixturefitmethod = "a",
                                param = list(concentration = 100))
icp.kde.region <- icp.torus.eval(icp.torus.kde, level = 0.02, eval.point = grid.torus())
p3<- data.frame(Chat = as.vector(icp.kde.region$Chat_kde), psi = icp.kde.region$psi, phi = icp.kde.region$phi) %>% 
  ggplot() + geom_contour(aes(phi,psi,z = Chat), color = "darkgrey") + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(expression(paste("KDE:", ~kappa, " = 100," ~alpha,  " = 0.02")))

icp.kde.region <- icp.torus.eval(icp.torus.kde, level = 0.1, eval.point = grid.torus())
p4<- data.frame(Chat = as.vector(icp.kde.region$Chat_kde), psi = icp.kde.region$psi, phi = icp.kde.region$phi) %>% 
  ggplot() + geom_contour(aes(phi,psi,z = Chat), color = "darkgrey") + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) + theme_bw() +
  xlab(expression(~ phi)) + ylab(expression(~ psi)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) + 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*pi/2, limits = c(-pi, pi),
                     labels = c(expression(- ~pi),expression(- ~pi/2),"0",expression(~pi/2),expression(~pi))) +
  ggtitle(expression(paste("KDE:", ~kappa, " = 100," ~alpha,  " = 0.1")))


#----------------------------------------------------------------------------------------
# Figure 5
#----------------------------------------------------------------------------------------
plot_grid(p1, p2, p3,p4, label_size = 12)
ggsave(paste(filename,"Figure_5.png",sep=""), width = 9, height = 6)



# choice of alpha and j browse options ------------------------------------
concentration.vec <- 10^( seq(log10(4),log10(150),length.out = 25) )
nkappa <- length(concentration.vec)
mu.vec1 <- concentration.vec
mu.vec2 <- concentration.vec

set.seed(102)
n <- nrow(data)
split.id <- rep(2,n)
split.id[ sample(n,floor(n/2)) ] <- 1
  
  
for (j in 1:nkappa){
  icp.torus.kde<- icp.torus.score(data,method = "kde",split.id = split.id ,param = list(concentration = concentration.vec[j]))
  icp.kde.region <- icp.torus.eval(icp.torus.kde, level = c(0.05,0.1), eval.point = grid.torus()) 
  mu.vec1[j] <- sum(icp.kde.region$Chat_kde[,1])/10000
  mu.vec2[j] <- sum(icp.kde.region$Chat_kde[,2])/10000
}

pICP_KDE_k <- data.frame(mu1 = mu.vec1, mu2 = mu.vec2, concentration = concentration.vec) %>% 
  pivot_longer(1:2, names_to = "Level") %>% mutate(Level = as.factor ( ifelse(Level=="mu1",0.05,0.1)) ,mu = value) %>% 
  ggplot(aes(concentration, mu, color = Level, group = Level)) + theme_bw() + geom_point(aes(shape = Level)) +
  geom_line()  + scale_x_log10() + xlab(expression(~ kappa)) + ylab(expression(~ mu)) +
  theme(axis.title.y = element_text(angle=0, vjust = 0.5))
pICP_KDE_k_grey <- pICP_KDE_k + scale_color_grey()


# alpha vs j for MIXTURE --------------------------------------------------
pICP_MIX_J <- out %>% filter( abs(alpha - 0.1) < 0.001 | abs(alpha - 0.05) < 0.001 ) %>% 
  mutate(Level = as.factor(ifelse( abs(alpha - 0.1) < 0.001, 0.1, 
                          0.05))) %>%  ggplot(aes(J,mu, color = Level)) + theme_bw() +
  geom_point(aes(shape = Level)) + geom_line() + ylab(expression(~ mu)) + theme(axis.title.y = element_text(angle=0, vjust = 0.5))
pICP_MIX_J_grey <- pICP_MIX_J + scale_color_grey()


#----------------------------------------------------------------------------------------
# Figure 7 - top
#----------------------------------------------------------------------------------------
plot_grid(pICP_KDE_k + ggtitle(expression(paste("KDE, over ", ~kappa))), 
          pICP_MIX_J+ ggtitle('Mixture, over J'), label_size = 12)
ggsave(paste(filename,"Figure_7_top.png",sep=""), width = 6, height = 3)

plot_grid(pICP_KDE_k_grey + ggtitle(expression(paste("KDE, over ", ~kappa))), 
          pICP_MIX_J_grey+ ggtitle('Mixture, over J'), label_size = 12)
ggsave(paste(filename,"Figure_7_top_grey.png",sep=""), width = 6, height = 3)

#----------------------------------------------------------------------------------------
# Figure 7 - bottom
#----------------------------------------------------------------------------------------
select_out <- out[which(out$J%%5==2 & out$J <30),]
select_out$J <- as.factor(select_out$J)

select_out %>% group_by(alpha) %>% ggplot(aes(alpha,criterion,color = J, shape = J)) + 
  geom_point() + theme_bw() + xlab(expression(~ alpha)) +
  geom_abline(slope = 0, intercept = out[out.index,4]) + scale_color_grey()
ggsave(paste(filename,"Figure_7_bottom_grey.png",sep=""), width = 6, height = 3)


out %>% ggplot(aes(x= alpha, y = criterion, color = J)) + geom_point() +
  theme_bw() +  xlab(expression(~ alpha)) + scale_colour_distiller(palette="Spectral") +
  geom_abline(slope = 0, intercept = out[out.index,4])
ggsave(paste(filename,"Figure_7_bottom.png",sep = ""), width = 6, height = 3)



best.option <- out %>% as.data.frame()  %>% arrange(criterion) %>% head(1) 

icp.torus <- l[[best.option$J]]

# Empirical coverage of Chat(16) for 6M15 and 6VXX data ----------------------------
Out.dt <- NULL
testing <- data15
testing <- cbind(testing$phi/180*pi,testing$psi/180*pi)
testing <- testing[-which(is.na(testing[,1])|is.na(testing[,2])),]

grid.test <- grid.torus()
testing.n <- nrow(testing)

# given a C, for each testing item, see if it is included in C. 

# aggregate C's  
alphavec <- 1:floor(n2/2) / n2
alphavec <- alphavec[alphavec <= 0.2]
N <- length(alphavec)
CC <- icp.torus.eval(icp.torus, level = alphavec)

C <- CC$Chat_e
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - testing[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_e"))

C <- CC$Chat_mix
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - testing[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_mix"))

C <- CC$Chat_max
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - testing[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_max"))


icp.torus.kde<- icp.torus.score(data, split.id = split.id,
                                method = "kde",
                                mixturefitmethod = "a",
                                param = list(concentration = 25))
icp.kde.region <- icp.torus.eval(icp.torus.kde, level = alphavec, eval.point = grid.torus())

C <- icp.kde.region$Chat_kde
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - testing[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_KDE"))
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, 
                           coverage = (1 - alphavec) - qnorm(1 - (0.02/2))*sqrt(alphavec*(1-alphavec)/testing.n), 
                           Type = "CI(98%)"))
Out.dt2 <- data.frame(Out.dt, Data = "HKU2")

# Empirical coverage (SADS-CoV) -------------------------------------------

Out.dt <- NULL
testing <- data16
testing <- cbind(testing$phi/180*pi,testing$psi/180*pi)
testing <- testing[-which(is.na(testing[,1])|is.na(testing[,2])),]

grid.test <- grid.torus()
testing.n <- nrow(testing)

# given a  C, for each testing item, see if it is included in C. 

# aggregate C's  
alphavec <- 1:floor(n2/2) / n2
alphavec <- alphavec[alphavec <= 0.2]
N <- length(alphavec)
CC <- icp.torus.eval(icp.torus, level = alphavec)

C <- CC$Chat_e
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - testing[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_e"))

C <- CC$Chat_mix
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - testing[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_mix"))

C <- CC$Chat_max
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - testing[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_max"))


icp.torus.kde<- icp.torus.score(data, split.id = split.id,
                                method = "kde",
                                mixturefitmethod = "a",
                                param = list(concentration = 25))
icp.kde.region <- icp.torus.eval(icp.torus.kde, level = alphavec, eval.point = grid.torus())

C <- icp.kde.region$Chat_kde
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - testing[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_KDE"))
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, 
                           coverage = (1 - alphavec) - qnorm(1 - (0.02/2))*sqrt(alphavec*(1-alphavec)/testing.n), 
                           Type = "CI(98%)"))
Out.dt3 <- data.frame(Out.dt, Data = "SADS-CoV")


# Empirical coverage (Covid19) --------------------------------------------
Out.dt <- NULL
testing <- dataXX
testing <- cbind(testing$phi/180*pi,testing$psi/180*pi)

grid.test <- grid.torus()
testing.n <- nrow(testing)

# given a C, for each testing item, see if it is included in C. 

# aggregate C's  
alphavec <- 1:floor(n2/2) / n2
alphavec <- alphavec[alphavec <= 0.2]
N <- length(alphavec)
CC <- icp.torus.eval(icp.torus, level = alphavec)

C <- CC$Chat_e
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - testing[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_e"))

C <- CC$Chat_mix
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - testing[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_mix"))

C <- CC$Chat_max
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - testing[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_max"))


icp.torus.kde<- icp.torus.score(data, split.id = split.id,
                                method = "kde",
                                mixturefitmethod = "a",
                                param = list(concentration = 25))
icp.kde.region <- icp.torus.eval(icp.torus.kde, level = alphavec, eval.point = grid.torus())

C <- icp.kde.region$Chat_kde
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - testing[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_KDE"))
Out.dt <- rbind(Out.dt, 
                data.frame(alpha = alphavec, 
                           coverage = (1 - alphavec) - qnorm(1 - (0.02/2))*sqrt(alphavec*(1-alphavec)/testing.n), 
                           Type = "CI(98%)"))

Out.dt4 <- data.frame(Out.dt, Data = "SARS-CoV-2")


Out_all <- rbind(Out.dt2,Out.dt3,Out.dt4)


q <- seq(0.8,0.99, length.out = 101)
q.se <- q*q
pointwiseCI <- data 



#----------------------------------------------------------------------------------------
# Figure 9
#----------------------------------------------------------------------------------------
Out_all %>% 
  filter(Type != "CI(98%)") %>% 
  ggplot(aes(1-alpha, coverage - 1 +alpha, color=Type, shape = Type)) + 
  guides(color = FALSE) +
  geom_line() + 
  facet_grid(Data ~ . ) + 
  theme_bw() + 
  geom_line(aes(1-alpha, coverage- 1 +alpha), 
            data = Out.dt %>% filter(Type == "CI(98%)"), color = "gray", size = 2) +
  geom_point() + 
  geom_abline(slope =0, intercept = 0) + 
  geom_vline(xintercept =  1- alphahat,linetype = "dotted" ) +
  scale_color_grey() + 
  xlab(expression(1 - ~ alpha)) + 
  ylab(expression(paste("coverage - (1 ",- ~ alpha ,")")))

 
ggsave(paste(filename,"Figure_9.png",sep=""), width = 8, height = 8/3*2)
#ggsave(paste(filename,"Empirical_coverage2a.png",sep=""), width = 8, height = 6)
 