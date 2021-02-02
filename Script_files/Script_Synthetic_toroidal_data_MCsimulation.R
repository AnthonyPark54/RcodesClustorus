##### Please set the working directory where you saved the provided files
# setwd("./R_Codes_for_submission")

library(MASS)
library(tidyverse)
library(cowplot)
library(ClusterR)
library(mclust)  
source('./Source_files/routines.R')

# Clustering by four methods (subfunction definition)------------------------------------------
##### The original function is defined in "Script_Synthetic_toroidal_data.R"
##### This function is a reduced version until the calculation of Rand index

Example_paper_supp_reduced <- function(J, dat1, dat1.test){
  
  data <- dat1[,1:2]
  data.test <- dat1.test
  
  data.test.label <- data.test[,3]
  data.test <- data.test[,1:2]
  
  predicted.label <- matrix(NA, nrow = nrow(data), ncol =4)
  colnames(predicted.label) <- c("kmeans_naive","kmeans_amb","Proposal_ehat","Proposal_out")
  
  # Kmeans_naive 
  kmeans.out <- KMeans_rcpp(data,clusters = J)
  predicted.label[,1] <- predict_KMeans(data.test, kmeans.out$centroids)
  
  # kmeans_ambient space
  kmeans.out <- KMeans_rcpp(cbind(cos(data),sin(data)),clusters = J)
  predicted.label[,2] <- predict_KMeans(cbind(cos(data.test),sin(data.test))
                                        , kmeans.out$centroids)
  
  # Clustering by our method
  
  # 1) Find alpha and J
  Jvec <- 3:35
  l <- list()
  
  # sample splitting; preparing data
  n <- nrow(data) 
  split.id <- rep(2,n)
  split.id[ sample(n,floor(n/2)) ] <- 1 
  
  for (j in Jvec){
    l[[j]] <- icp.torus.score(as.matrix(data), split.id = split.id,
                              method = "mixture",
                              mixturefitmethod = "a",
                              param = list(J = j))
  }
  
  n2 <- l[[10]]$n2 
  alphavec <- 1:floor(n2/2) / n2
  N <- length(alphavec)
  out <- data.frame()
  
  for (j in Jvec){
    Mvec <- alphavec
    a <- icp.torus.eval(l[[j]], level = alphavec, eval.point = grid.torus())
    for (i in 1:N){
      Mvec[i] <- sum(a$Chat_e[,i])/10000
    }
    out <- rbind(out, data.frame(alpha = alphavec, J = j, mu = Mvec, criterion = alphavec + Mvec))
  }
  
  out.index <- which.min(out$criterion)
  out[out.index,]
  
  Jhat <- out[out.index,2]
  alphahat <- out[out.index,1]
  icp.torus <- l[[Jhat]] 
  
  c <- cluster.assign.torus(data.test, icp.torus, level = alphahat)

  predicted.label[,3] <- c$cluster.id.by.ehat
  predicted.label[,4] <- c$cluster.id.outlier
  
  aa <- rep(0,4)
  
  for (j in 1:4){ 
    aa[j] <- adjustedRandIndex(predicted.label[,j],data.test.label)
  }
  
  list(choice = out[out.index,], labels = predicted.label, Rand = aa)
}

##### Simulation: MC experiments with data examples in Figure 10 #####
#####             Produce the result in Table 1                  #####

set.seed(202101) # Global seed for reproduction
REP <- 100
Seed <- sample(10000, REP)
Data1.Result.Rand <- data.frame(Repeat = 1:REP,
                                Naive.Kmeans = rep(0, REP),
                                Extrinsic.Kmeans = rep(0, REP),
                                A_e = rep(0, REP),
                                A_o = rep(0, REP))
Data2.Result.Rand <- data.frame(Repeat = 1:REP,
                                Naive.Kmeans = rep(0, REP),
                                Extrinsic.Kmeans = rep(0, REP),
                                A_e = rep(0, REP),
                                A_o = rep(0, REP))

for (rr in 1:REP){
  set.seed(Seed[rr])
  # Data generation
  # Data set 1 ------------------------------------------------------------
  ## Five true clusters ---------------------------------------------------
  
  Mu1 <- c(3,0)
  Mu2 <- c(2,2)
  Mu3 <- c(1,4)
  Sigma1 <- matrix(c(0.1,0.05,0.05,0.2),2,2)
  Sigma2 <- matrix(c(0.1,0,0,0.01),2,2)
  Sigma3 <- matrix(c(0.01,0,0,0.1),2,2)
  
  unidata <- cbind(2*runif(50, -0.5, 0.5), 0.5*runif(50, -0.5, 0.5))
  data.unif <- cbind(unidata[,1]+ 0, unidata[,2] + 1)
  data.diamond <- t(matrix(c(cos(-pi/4),-sin(-pi/4), 
                             sin(-pi/4),cos(-pi/4)),2,2) %*% t(unidata)) +cbind(rep(5, 50), rep(3, 50))
  
  Example1 <- rbind(mvrnorm(n=70, Mu1, Sigma1), 
                    mvrnorm(n=50, Mu2, Sigma2), 
                    mvrnorm(n=50, Mu3, Sigma3), 
                    data.unif, 
                    data.diamond)
  Example1 <- on.torus(Example1)
  label <- c(rep(1,70),rep(2,50),
             rep(3,50),
             rep(4,50),
             rep(5,50))
  dat1 <- cbind( as.data.frame(Example1) , as.factor(label))
  colnames(dat1) <- c("phi","psi","label")
  
  unidata <- cbind(2*runif(50, -0.5, 0.5), 0.5*runif(50, -0.5, 0.5))
  data.unif <- cbind(unidata[,1]+ 0, unidata[,2] + 1)
  data.diamond <- t(matrix(c(cos(-pi/4),-sin(-pi/4), 
                             sin(-pi/4),cos(-pi/4)),2,2) %*% t(unidata)) +cbind(rep(5, 50), rep(3, 50))
  
  Example1 <- rbind(mvrnorm(n=70, Mu1, Sigma1), 
                    mvrnorm(n=50, Mu2, Sigma2), 
                    mvrnorm(n=50, Mu3, Sigma3), 
                    data.unif, 
                    data.diamond)
  Example1 <- on.torus(Example1)
  label <- c(rep(1,70),rep(2,50),
             rep(3,50),
             rep(4,50),
             rep(5,50))
  dat1.test <- cbind( as.data.frame(Example1) , as.factor(label))
  colnames(dat1.test) <- c("phi","psi","label")
  
  
  # Data set 2 ------------------------------------------------------------
  ## Two true clusters (Ball and L shape) ---------------------------------
  
  Mu <- c(1,5)
  Sigma <- matrix(c(0.05,0,0,0.05),2,2)
  
  unidata.unif <- cbind(runif(50, 0, 2*pi), runif(50, 0, 2*pi))
  data.unif1 <- cbind(1.25*runif(150, -1, 1)+ 3.25, 0.5*runif(150, -1, 1) + 1)
  data.unif2 <- cbind(0.5*runif(200, -1, 1)+ 5, 2*runif(200, -1, 1) + 2.5)
  
  Example1 <- rbind(mvrnorm(n=100, Mu, Sigma), 
                    unidata.unif, 
                    data.unif1, 
                    data.unif2)
  Example1 <- Example1 + 2
  
  Example1 <- on.torus(Example1)
  label <- c(rep(1,100),
             rep(3,50),
             rep(2,350))
  dat2 <- cbind( as.data.frame(Example1) , as.factor(label))
  colnames(dat2) <- c("phi","psi","label")
  unidata.unif <- cbind(runif(50, 0, 2*pi), runif(50, 0, 2*pi))
  data.unif1 <- cbind(1.25*runif(150, -1, 1)+ 3.25, 0.5*runif(150, -1, 1) + 1)
  data.unif2 <- cbind(0.5*runif(200, -1, 1)+ 5, 2*runif(200, -1, 1) + 2.5)
  
  Example1 <- rbind(mvrnorm(n=100, Mu, Sigma), 
                    unidata.unif, 
                    data.unif1, 
                    data.unif2)
  Example1 <- Example1 + 2
  
  Example1 <- on.torus(Example1)
  label <- c(rep(1,100),
             rep(3,50),
             rep(2,350))
  dat2.test <- cbind( as.data.frame(Example1) , as.factor(label))
  colnames(dat2.test) <- c("phi","psi","label")
  
  # Calculation of Rand index (for 4 methods) 
  J <- 5
  dat1.results <- Example_paper_supp_reduced(J, dat1, dat1.test)
  J <- 2
  dat2.results <- Example_paper_supp_reduced(J, dat2, dat2.test)
  
  Data1.Result.Rand[rr, 2:5] <- dat1.results$Rand
  Data2.Result.Rand[rr, 2:5] <- dat2.results$Rand
}

apply(Data1.Result.Rand[,c("Naive.Kmeans", "Extrinsic.Kmeans", "A_e", "A_o")], 2, mean)
apply(Data2.Result.Rand[,c("Naive.Kmeans", "Extrinsic.Kmeans", "A_e", "A_o")], 2, mean)
apply(Data1.Result.Rand[,c("Naive.Kmeans", "Extrinsic.Kmeans", "A_e", "A_o")], 2, sd) / 10
apply(Data2.Result.Rand[,c("Naive.Kmeans", "Extrinsic.Kmeans", "A_e", "A_o")], 2, sd) / 10

write.table(Data1.Result.Rand, "./Simulation_Result/Simulation_example1_Rand.txt")
write.table(Data2.Result.Rand, "./Simulation_Result/Simulation_example2_Rand.txt")



# Median and MAD ----------------------------------------------------------

a1 <- read.table("./Simulation_Result/Simulation_example1_Rand.txt", header = TRUE)
a2 <- read.table("./Simulation_Result/Simulation_example2_Rand.txt", header = TRUE)
apply(a1[,2:5], 2, median)
apply(a1[,2:5], 2, IQR)
apply(a1[,2:5], 2, mean)
apply(a1[,2:5], 2, sd)
apply(a2[,2:5], 2, median)
apply(a2[,2:5], 2, IQR)
apply(a2[,2:5], 2, mean)
apply(a2[,2:5], 2, sd)
