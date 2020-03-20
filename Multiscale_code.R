2# EVA2019 Challenge Data Analysis R Code for Multiscale Team
# A brief description for the method is available in readme.txt. 

##############################################################
##############################################################
######                    Clustering                    ######
##############################################################
##############################################################

# At first, to reflect local feature well, and also for parallel computation, 
# clustering procedure would be done primarily.
# In here, we use PAM based clustering method for extreme data (Bernard et al., 2013).

library(xts) #xts: for time-series data
library(foreach)
library(doParallel) #foreach and doParallel: for parallel data code
### Read data
load('~/Dropbox/Multiscale/DATA_TRAINING.RData')

### Read abnomalie data
# X.min_fill is missing-imputed X.min.
# Missing imputation procedure will be described further. 
load('~/Dropbox/Multiscale/Xmin.Rdata')
load('~/Dropbox/Multiscale/Xmin_fill.Rdata')

# Because of the size of the data, we can execute clusting algorithm with whole dataset subsetting.
X.min_fill <- X.min_fill[c(6571:nrow(X.min_fill)),] #from 1985 to 2003

# We also have to make a subset of locations
set.seed(1)
n.geo <- sample(x=c(1:ncol(X.min_fill)), size=5000) #to reduce the number of points, we use random subset

########################################
##Philippe's method
##Use ClusterMax package
##Clustering of Maxima: Spatial Dependencies among Heavy Rainfall in France, Journal of Climate, 2013
########################################
#for Philippe's method
#install.packages("~/Dropbox/R files/Papers/Philippe/ClusterMax_1.0.tar.gz", type="source", repos=NULL)
library(ClusterMax)
library(cluster) #k-mediod clustering

# K: the number of clusters
# Want to choose optimal number of cluster (modified 20190504)

K.vec <- 30:70 # candidate clustersets, we want to divide the data small enough to compute, large enough to effective tail coefficient estimation
clus.avg.width.list <- list()
avg.width.list <- list()
#result.obj <- list()

# setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

result.obj <- foreach(jj=1:K.vec, .packages=c("ClusterMax", "cluster") ) %dopar%  {
  PAMmado_imsi=PAMfmado.R(X.min_fill[,n.geo],K.vec[jj]) #PAMfmado.R: main function
  result.obj = list(clus.avg.widths=PAMmado_imsi$silinfo$clus.avg.widths, avg.width=PAMmado_imsi$silinfo$avg.width)
}

#stop cluster
stopCluster(cl)

saveRDS(result.obj, "~/Dropbox/EVA2019/optimalcluster(from2003)(2).RDS")
result.obj <- readRDS("~/Dropbox/EVA2019/optimalcluster(from2003).RDS")

# Optimal cluster number is the one 
# which cluster has the mininum value of avg. silhouette coefficient

cluster2003 <- readRDS("~/Dropbox/Multiscale/Methods/Clustering/optimalcluster(from2003)(2).RDS")
which.min(sapply(cluster2003, function(x) return(x$avg.width))) 


########################################
##Geographically constrained version of clustering
########################################
silmado[,"cluster"]
n.geo[as.numeric(rownames(silmado))]

idx.K <- rep(NA, ncol(X.min_fill))
idx.K[n.geo[as.numeric(rownames(silmado.final))]] <- silmado.final[,"cluster"]
coords.idx <- sort(n.geo[as.numeric(rownames(silmado.final))]) 
coords.idx.match <-  match(coords.idx, n.geo[as.numeric(rownames(silmado.final))])
for(jj in 1:length(idx.K)){
  if(is.na(idx.K[jj])){
    #find the nearest point
    result_imsi <- apply(loc[coords.idx,], 1, function(x) sqrt(sum((x-loc[jj,])^2))  )
    idx.imsi <- which(result_imsi==min(result_imsi))
    
    if(length(idx.imsi)!=1){
      #We assign cluser index by searching the majority of nearby five points
      idx.K[jj] <- as.numeric(names(which.max(table(silmado.final[coords.idx.match[order(result_imsi)[c(1:5)]],"cluster"]))) )
      
    }else{
      idx.K[jj] <- silmado.final[coords.idx.match[idx.imsi],"cluster"]
    }
    
  }
}

saveRDS(idx.K, "~/Dropbox/Multiscale/Methods/Clustering/idx(68).RDS")

##############################################################
##############################################################
######                 Missing imputation               ######
##############################################################
##############################################################

# For missing imputation: 
# # Just impute as neighborhood minimum ignoring missing at neighbors.
# # If there isn't any observation at neighborhood, expand the neighbor with radius=100. 
# # If there still isn't any observation at neighborhood, expand the neighbor with radius=150. 
# # If there still isn't any observation at neighborhood, just use the mean of observed X.min. 

X.min_fill2 <- matrix(nrow=nrow(anom.training),ncol=ncol(anom.training))

globalMean = mean(X.min[which(!is.infinite(X.min))], na.rm=TRUE)
for(j in 1:ncol(anom.training)){ # for each location
  print(j) 
  dist.j <- drop(rdist.earth(x1=matrix(loc[j,],nrow=1),x2=loc,miles=FALSE)) # compute the distance between the j-th location and all other locations
  loc.nei <- which(dist.j<=radius) # find its neighbors (within the specified radius of 50km)
  loc.nei2 <- which(dist.j<=radius*2)
  loc.nei3 <- which(dist.j<=radius*3)
  for(i in 1:nrow(anom.training)){ # for each day
    week.i <- i+c(-3:3) # define the week for i-th day (i.e., temporal neighborhood)
    week.i <- week.i[week.i>=1 & week.i<=nrow(anom.training)] # make sure the week contains valid indices
    tmp.min <- min(anom.training[week.i,loc.nei], na.rm=T) #####
    if(is.infinite(tmp.min)){
      week.i2 = i + c(-15:15)
      week.i2 = week.i2[week.i2>=1 & week.i2<=nrow(anom.training)]
      tmp.min = min(anom.training[week.i2,loc.nei2], na.rm=T) #####
      if(is.infinite(tmp.min)){
        week.i3 <- i+c(-30:30)
        week.i3 <- week.i3[week.i3>=1 & week.i3<=nrow(anom.training)]
        tmp.min <- min(anom.training[week.i3,loc.nei3], na.rm=T) #####
        if(is.infinite(tmp.min)){
          tmp.min = globalMean
        }
      }
    }
    X.min_fill2[i,j] <- tmp.min
  }
}





##############################################################
##############################################################
######                    QFM fitting                   ######
##############################################################
##############################################################

# We have used Quantile Factor Model to predict distribution. 

load('~/Dropbox/Laboratory/1Papers/EVA2019/DATA_TRAINING.RData')

# Load missing-imputed dataset X.min_fill2. 
# Missing imputation procedure will be described further.

load('~/Dropbox/Laboratory/1Papers/EVA2019/Multiscale/Xmin_fill2.Rdata')

# Read optimal clustered result
# idx is 16703-dim vector, whose element is a cluster number of each location. 

idx <- readRDS("~/Dropbox/Laboratory/1Papers/EVA2019/Multiscale/Methods/Clustering/idx(68).RDS") 

library(quantreg)
library(foreach) #needed for parallel computing
library(doParallel) #needed for parallel computing
library(abind) #arraybind

# Function for projection matrix 
proj.mat <- function(A){
	return(A %*% solve(t(A) %*% A) %*% t(A))
}

xk <- -1+c(1:400)/100 
n.xk <- length(xk)
n.validation <- length(index.validation)
prediction <- matrix(nrow=n.validation, ncol=n.xk)


# Find time and location index of validation
# Convert index.validation vector to time / location index. 
poo.time <- matrix(rep(1:11315, 16703), ncol=16703)
poo.loc <- matrix(rep(1:16703, 11315), ncol=16703, byrow=T) 

index.time.validation <- poo.time[index.validation]
index.loc.validation <- poo.loc[index.validation]

rm(poo.time); rm(poo.loc) # Cleaning memory up 

# Levels and extreme levels of which quantile value would be estimated.
tau.vec <- seq(from=0.05, to=0.95, by=0.05)
tau.vec.ex <- c(0.96, 0.97, 0.98, 0.99, 0.995, 0.999, 0.9995)
pred.quant.final <- matrix(nrow=n.validation, ncol=length(tau.vec)+length(tau.vec.ex))

rm(anom.training)  # Cleaning memory up 

# Transform quantile process to cdf
# For this, we linearly interpolate quantiles, and find its inverse. 
# cdf.est gets x and return F(x), where F is the estimated cdf. 
cdf.est <- function(x, pred.quant, tau.vec=tau.vec){
	# pred.quant : (T x S) x K, K : number of quantile levels. 
	# in this function, pred.quant is K-dim vector. 
	# Then this function would be applied to the ((ST)XK) object pred.quant. 
	if(x < pred.quant[1]){
		return(0)
	}else if(x > max(pred.quant)){
		return(1)
	}else{
		index <- rank(c(x, pred.quant), ties="min")[1]
		# if(index==50){
			# return(tau.vec[index])
		# }else{
			return(tau.vec[index-1]+(x-pred.quant[index-1])*(tau.vec[index]-tau.vec[index-1])/(pred.quant[index]-pred.quant[index-1]))
		# }
	}
}

# For the estimation of extreme level quantiles, we extrapolate intermediate quantiles. 
# To estimate tail index, we implement Hill estimator (Hill, 1975). 
# gamma.cal is the function calculating Hill estimators. 
# gamma.hat includes estimated tail index at each location.
N <- 16703; T <- 11315
gamma.hat <- vector(length=N)
gamma.cal <- function(data, k=floor(0.5*length(data)^(1/3))){
	data.sorted <- sort(data, dec=TRUE)
	return(mean(log(data.sorted[1:k]))-log(data.sorted[k]))
}
gamma.hat <- apply(X.min_fill2, 2, gamma.cal, k=floor(0.5*T^(1/3)))


# r is the predetermined number of factors. 
# The value is found via proper validation procedure, which is not described in this file. 
# We partition whole region into 5 groups, and in each group find the number r of factors making twCRPS validation score the smallest. 
# Following vectors include member index (idx) in each group. 
cluster.sub1 <- c(54,38,64,20,16,19,2,12,66,41,57,55)
cluster.sub2 <- c(67,7,15,33,11,24,37,61,32,49,44,46,45,65,5,59,17,34)
cluster.sub3 <- c(26,27,30,9,63,13,58,28,48,53,50,4,68,62,14)
cluster.sub4 <- c(3,29,43,36,60,52,42,51,10,8,22,6,47)
cluster.sub5 <- c(40,35,56,23,31,18,25,39,21,1)

# Fitted r values are following: 
# # Group 1 (cluster.sub1) : r = 4 
# # Group 2 (cluster.sub2) : r = 3
# # Group 3 (cluster.sub3) : r = 4 
# # Group 4 (cluster.sub4) : r = 3
# # Group 5 (cluster.sub5) : r = 5

# Group 1
r <- 4
for(clust.num in cluster.sub1){	# cluster number
	#setup parallel backend to use many processors
	cores=detectCores()
	cl <- makeCluster(cores[1]-1) #not to overload your computer
	registerDoParallel(cl)

	# clust.loc : station number which belogs to clust.num-th cluster
	# X.min_fill3 : (missing-imputed) data object at clust.loc stations
	# temp : PCA result of X.min_fill3; initial value of the iteration
	
	clust.loc <- which(idx==clust.num) 
	X.min_fill3 <- X.min_fill2[,clust.loc] #subset (by clustering)
	temp <- princomp(X.min_fill3) 
	#for(ind in 1:length(tau.vec)){
	pred.quant.imsi <- foreach(ind = 1:length(tau.vec)) %dopar% {
	  library(quantreg)
	  print(paste("###########", ind, "th QFM estimation###########", sep=""))
	  tau <- tau.vec[ind]
	  
	  # To minimize the target function of QFM, we use fully iterative approach. 
	  # L.hat: current loading
	  # L.hat.new : updated loading
	  # F.hat : current factor
	  # F.hat.new : updated factor
	  # maximum iteration number: 200 
	  # If spanned space by updated factors is near to that by previous ones enough, stop the iteration.
	  L.hat <- temp$loadings[,1:r]
	  F.hat <- temp$scores[,1:r]
	  
	  F.hat.new <- matrix(0, nrow=nrow(F.hat), ncol=r)
	  L.hat.new <- matrix(0, nrow=nrow(L.hat), ncol=r)
	  
	  iter <- 1
	  
	  while(1){
	    for(t in 1:nrow(F.hat)){
	      F.hat.new[t,] <- rq(X.min_fill3[t,] ~ L.hat-1, tau=tau)$coeff
	    }
	    for(i in 1:nrow(L.hat)){
	      L.hat.new[i,] <- rq(X.min_fill3[,i] ~ F.hat.new-1, tau=tau)$coeff
	    }

	    if(sum((proj.mat(F.hat.new)-proj.mat(F.hat))^2)<=1e-2){
	      break;
	    }
	    iter <- iter+1
	    F.hat <- F.hat.new
	    L.hat <- L.hat.new
	    if(iter==200){
	      print(paste("Max iteration at tau=", tau))
	      break;
	    }
	  }
	  pred.quant.imsi <- list(F=F.hat.new, L=L.hat.new)
	}
	

	
	
	# # pred.quant : QFM result
	# # pred.quant.sort : sorted result (to make quantiles non-crossing) (cf. Chernozhukove et al. (2010))
	
	# N <- dim(X.min_fill3)[2]
	# T <- dim(X.min_fill3)[1]
	# pred.quant <- array(unlist(pred.quant.imsi), dim=c(T,N,length(tau.vec))) 
	# pred.quant.sort <- aperm(apply(pred.quant, c(1,2), sort), c(2,3,1))	
	
	# # pred.quant.ex : extremal level quantile estimates
		
	 # pred.quant.ex <- array(dim=c(T,N,length(tau.vec.ex)))
	  # for(i in 1:length(tau.vec.ex)){
	    # tau.ex <- tau.vec.ex[i]
	    # pred.quant.ex[,,i] <- (((1-tau.vec[length(tau.vec)])/(1-tau.ex))^(matrix(rep(gamma.hat[clust.loc], T), ncol=N, byrow=TRUE))*pred.quant.sort[,,length(tau.vec)])*(pred.quant.sort[,,length(tau.vec)]>=0) + pred.quant.sort[,,length(tau.vec)]*(pred.quant.sort[,,length(tau.vec)]<0)
	  # }
	
	# tau.vec.bind <- c(tau.vec, tau.vec.ex)
	
	# # idx.index.validation : which index belongs to the current cluster among validation indices? 
	# # For the members of idx.index.validation, inversion of estimated quantile process is executed. 
	# # prediction.forsave is the final prediction object at the current cluster.
	
	# idx.index.invalidation <- NULL
	# for(i in 1:length(clust.loc)){ 
		# idx.index <- which(clust.loc[i]==index.loc.validation)
		# idx.index.invalidation <- c(idx.index.invalidation, idx.index)
		# for(k in idx.index){
			# poo <- index.time.validation[k]
			# for(j in 1:length(tau.vec)){
				# pred.quant.final[k, j] <- pred.quant.sort[poo, i, j]
			# }
			# for(j in 1:length(tau.vec.ex)){
				# pred.quant.final[k, length(tau.vec)+j] <- pred.quant.ex[poo, i, j]
			# }
		# prediction[k,] <- sapply(xk, cdf.est, pred.quant=pred.quant.final[k,], tau.vec=tau.vec.bind)
		# }
	# }
	
	# prediction.forsave <- list(idx.index.invalidation=idx.index.invalidation, prediction.idx = prediction[idx.index.invalidation,])
	
	# filename <- paste("prediction(idx", clust.num, "of68).RDS", sep="")
	# saveRDS(prediction.forsave, file=filename)
	
	filename.f <- paste("factor(idx", clust.num, "of68).RDS", sep="")
	saveRDS(pred.quant.imsi, file=filename.f)
	
	stopCluster(cl)
}

# Group 2
r <- 3
for(clust.num in cluster.sub2){	# cluster number
	#setup parallel backend to use many processors
	cores=detectCores()
	cl <- makeCluster(cores[1]-1) #not to overload your computer
	registerDoParallel(cl)

	# clust.loc : station number which belogs to clust.num-th cluster
	# X.min_fill3 : (missing-imputed) data object at clust.loc stations
	# temp : PCA result of X.min_fill3; initial value of the iteration
	
	clust.loc <- which(idx==clust.num) 
	X.min_fill3 <- X.min_fill2[,clust.loc] #subset (by clustering)
	temp <- princomp(X.min_fill3) 
	#for(ind in 1:length(tau.vec)){
	pred.quant.imsi <- foreach(ind = 1:length(tau.vec)) %dopar% {
	  library(quantreg)
	  print(paste("###########", ind, "th QFM estimation###########", sep=""))
	  tau <- tau.vec[ind]
	  
	  # To minimize the target function of QFM, we use fully iterative approach. 
	  # L.hat: current loading
	  # L.hat.new : updated loading
	  # F.hat : current factor
	  # F.hat.new : updated factor
	  # maximum iteration number: 200 
	  # If spanned space by updated factors is near to that by previous ones enough, stop the iteration.
	  L.hat <- temp$loadings[,1:r]
	  F.hat <- temp$scores[,1:r]
	  
	  F.hat.new <- matrix(0, nrow=nrow(F.hat), ncol=r)
	  L.hat.new <- matrix(0, nrow=nrow(L.hat), ncol=r)
	  
	  iter <- 1
	  
	  while(1){
	    for(t in 1:nrow(F.hat)){
	      F.hat.new[t,] <- rq(X.min_fill3[t,] ~ L.hat-1, tau=tau)$coeff
	    }
	    for(i in 1:nrow(L.hat)){
	      L.hat.new[i,] <- rq(X.min_fill3[,i] ~ F.hat.new-1, tau=tau)$coeff
	    }

	    if(sum((proj.mat(F.hat.new)-proj.mat(F.hat))^2)<=1e-2){
	      break;
	    }
	    iter <- iter+1
	    F.hat <- F.hat.new
	    L.hat <- L.hat.new
	    if(iter==200){
	      print(paste("Max iteration at tau=", tau))
	      break;
	    }
	  }

	  pred.quant.imsi <- list(F=F.hat.new, L=L.hat.new)
	}
	

	
	
	# # pred.quant : QFM result
	# # pred.quant.sort : sorted result (to make quantiles non-crossing) (cf. Chernozhukove et al. (2010))
	
	# N <- dim(X.min_fill3)[2]
	# T <- dim(X.min_fill3)[1]
	# pred.quant <- array(unlist(pred.quant.imsi), dim=c(T,N,length(tau.vec))) 
	# pred.quant.sort <- aperm(apply(pred.quant, c(1,2), sort), c(2,3,1))	
	
	# # pred.quant.ex : extremal level quantile estimates
		
	 # pred.quant.ex <- array(dim=c(T,N,length(tau.vec.ex)))
	  # for(i in 1:length(tau.vec.ex)){
	    # tau.ex <- tau.vec.ex[i]
	    # pred.quant.ex[,,i] <- (((1-tau.vec[length(tau.vec)])/(1-tau.ex))^(matrix(rep(gamma.hat[clust.loc], T), ncol=N, byrow=TRUE))*pred.quant.sort[,,length(tau.vec)])*(pred.quant.sort[,,length(tau.vec)]>=0) + pred.quant.sort[,,length(tau.vec)]*(pred.quant.sort[,,length(tau.vec)]<0)
	  # }
	
	# tau.vec.bind <- c(tau.vec, tau.vec.ex)
	
	# # idx.index.validation : which index belongs to the current cluster among validation indices? 
	# # For the members of idx.index.validation, inversion of estimated quantile process is executed. 
	# # prediction.forsave is the final prediction object at the current cluster.
	
	# idx.index.invalidation <- NULL
	# for(i in 1:length(clust.loc)){ 
		# idx.index <- which(clust.loc[i]==index.loc.validation)
		# idx.index.invalidation <- c(idx.index.invalidation, idx.index)
		# for(k in idx.index){
			# poo <- index.time.validation[k]
			# for(j in 1:length(tau.vec)){
				# pred.quant.final[k, j] <- pred.quant.sort[poo, i, j]
			# }
			# for(j in 1:length(tau.vec.ex)){
				# pred.quant.final[k, length(tau.vec)+j] <- pred.quant.ex[poo, i, j]
			# }
		# prediction[k,] <- sapply(xk, cdf.est, pred.quant=pred.quant.final[k,], tau.vec=tau.vec.bind)
		# }
	# }
	
	# prediction.forsave <- list(idx.index.invalidation=idx.index.invalidation, prediction.idx = prediction[idx.index.invalidation,])
	
	# filename <- paste("prediction(idx", clust.num, "of68).RDS", sep="")
	# saveRDS(prediction.forsave, file=filename)
	
	filename.f <- paste("factor(idx", clust.num, "of68).RDS", sep="")
	saveRDS(pred.quant.imsi, file=filename.f)
	
	stopCluster(cl)
}

# Group 3
r <- 4
for(clust.num in cluster.sub3){	# cluster number
	#setup parallel backend to use many processors
	cores=detectCores()
	cl <- makeCluster(cores[1]-1) #not to overload your computer
	registerDoParallel(cl)

	# clust.loc : station number which belogs to clust.num-th cluster
	# X.min_fill3 : (missing-imputed) data object at clust.loc stations
	# temp : PCA result of X.min_fill3; initial value of the iteration
	
	clust.loc <- which(idx==clust.num) 
	X.min_fill3 <- X.min_fill2[,clust.loc] #subset (by clustering)
	temp <- princomp(X.min_fill3) 
	#for(ind in 1:length(tau.vec)){
	pred.quant.imsi <- foreach(ind = 1:length(tau.vec)) %dopar% {
	  library(quantreg)
	  print(paste("###########", ind, "th QFM estimation###########", sep=""))
	  tau <- tau.vec[ind]
	  
	  # To minimize the target function of QFM, we use fully iterative approach. 
	  # L.hat: current loading
	  # L.hat.new : updated loading
	  # F.hat : current factor
	  # F.hat.new : updated factor
	  # maximum iteration number: 200 
	  # If spanned space by updated factors is near to that by previous ones enough, stop the iteration.
	  L.hat <- temp$loadings[,1:r]
	  F.hat <- temp$scores[,1:r]
	  
	  F.hat.new <- matrix(0, nrow=nrow(F.hat), ncol=r)
	  L.hat.new <- matrix(0, nrow=nrow(L.hat), ncol=r)
	  
	  iter <- 1
	  
	  while(1){
	    for(t in 1:nrow(F.hat)){
	      F.hat.new[t,] <- rq(X.min_fill3[t,] ~ L.hat-1, tau=tau)$coeff
	    }
	    for(i in 1:nrow(L.hat)){
	      L.hat.new[i,] <- rq(X.min_fill3[,i] ~ F.hat.new-1, tau=tau)$coeff
	    }

	    if(sum((proj.mat(F.hat.new)-proj.mat(F.hat))^2)<=1e-2){
	      break;
	    }
	    iter <- iter+1
	    F.hat <- F.hat.new
	    L.hat <- L.hat.new
	    if(iter==200){
	      print(paste("Max iteration at tau=", tau))
	      break;
	    }
	  }

	  pred.quant.imsi <- list(F=F.hat.new, L=L.hat.new)
	}
	

	
	
	# # pred.quant : QFM result
	# # pred.quant.sort : sorted result (to make quantiles non-crossing) (cf. Chernozhukove et al. (2010))
	
	# N <- dim(X.min_fill3)[2]
	# T <- dim(X.min_fill3)[1]
	# pred.quant <- array(unlist(pred.quant.imsi), dim=c(T,N,length(tau.vec))) 
	# pred.quant.sort <- aperm(apply(pred.quant, c(1,2), sort), c(2,3,1))	
	
	# # pred.quant.ex : extremal level quantile estimates
		
	 # pred.quant.ex <- array(dim=c(T,N,length(tau.vec.ex)))
	  # for(i in 1:length(tau.vec.ex)){
	    # tau.ex <- tau.vec.ex[i]
	    # pred.quant.ex[,,i] <- (((1-tau.vec[length(tau.vec)])/(1-tau.ex))^(matrix(rep(gamma.hat[clust.loc], T), ncol=N, byrow=TRUE))*pred.quant.sort[,,length(tau.vec)])*(pred.quant.sort[,,length(tau.vec)]>=0) + pred.quant.sort[,,length(tau.vec)]*(pred.quant.sort[,,length(tau.vec)]<0)
	  # }
	
	# tau.vec.bind <- c(tau.vec, tau.vec.ex)
	
	# # idx.index.validation : which index belongs to the current cluster among validation indices? 
	# # For the members of idx.index.validation, inversion of estimated quantile process is executed. 
	# # prediction.forsave is the final prediction object at the current cluster.
	
	# idx.index.invalidation <- NULL
	# for(i in 1:length(clust.loc)){ 
		# idx.index <- which(clust.loc[i]==index.loc.validation)
		# idx.index.invalidation <- c(idx.index.invalidation, idx.index)
		# for(k in idx.index){
			# poo <- index.time.validation[k]
			# for(j in 1:length(tau.vec)){
				# pred.quant.final[k, j] <- pred.quant.sort[poo, i, j]
			# }
			# for(j in 1:length(tau.vec.ex)){
				# pred.quant.final[k, length(tau.vec)+j] <- pred.quant.ex[poo, i, j]
			# }
		# prediction[k,] <- sapply(xk, cdf.est, pred.quant=pred.quant.final[k,], tau.vec=tau.vec.bind)
		# }
	# }
	
	# prediction.forsave <- list(idx.index.invalidation=idx.index.invalidation, prediction.idx = prediction[idx.index.invalidation,])
	
	# filename <- paste("prediction(idx", clust.num, "of68).RDS", sep="")
	# saveRDS(prediction.forsave, file=filename)
	
	filename.f <- paste("factor(idx", clust.num, "of68).RDS", sep="")
	saveRDS(pred.quant.imsi, file=filename.f)
	
	stopCluster(cl)
}

# Group 4
r <- 3
for(clust.num in cluster.sub4){	# cluster number
	#setup parallel backend to use many processors
	cores=detectCores()
	cl <- makeCluster(cores[1]-1) #not to overload your computer
	registerDoParallel(cl)

	# clust.loc : station number which belogs to clust.num-th cluster
	# X.min_fill3 : (missing-imputed) data object at clust.loc stations
	# temp : PCA result of X.min_fill3; initial value of the iteration
	
	clust.loc <- which(idx==clust.num) 
	X.min_fill3 <- X.min_fill2[,clust.loc] #subset (by clustering)
	temp <- princomp(X.min_fill3) 
	#for(ind in 1:length(tau.vec)){
	pred.quant.imsi <- foreach(ind = 1:length(tau.vec)) %dopar% {
	  library(quantreg)
	  print(paste("###########", ind, "th QFM estimation###########", sep=""))
	  tau <- tau.vec[ind]
	  
	  # To minimize the target function of QFM, we use fully iterative approach. 
	  # L.hat: current loading
	  # L.hat.new : updated loading
	  # F.hat : current factor
	  # F.hat.new : updated factor
	  # maximum iteration number: 200 
	  # If spanned space by updated factors is near to that by previous ones enough, stop the iteration.
	  L.hat <- temp$loadings[,1:r]
	  F.hat <- temp$scores[,1:r]
	  
	  F.hat.new <- matrix(0, nrow=nrow(F.hat), ncol=r)
	  L.hat.new <- matrix(0, nrow=nrow(L.hat), ncol=r)
	  
	  iter <- 1
	  
	  while(1){
	    for(t in 1:nrow(F.hat)){
	      F.hat.new[t,] <- rq(X.min_fill3[t,] ~ L.hat-1, tau=tau)$coeff
	    }
	    for(i in 1:nrow(L.hat)){
	      L.hat.new[i,] <- rq(X.min_fill3[,i] ~ F.hat.new-1, tau=tau)$coeff
	    }

	    if(sum((proj.mat(F.hat.new)-proj.mat(F.hat))^2)<=1e-2){
	      break;
	    }
	    iter <- iter+1
	    F.hat <- F.hat.new
	    L.hat <- L.hat.new
	    if(iter==200){
	      print(paste("Max iteration at tau=", tau))
	      break;
	    }
	  }

	  pred.quant.imsi <- list(F=F.hat.new, L=L.hat.new)
	}
	

	
	# # pred.quant : QFM result
	# # pred.quant.sort : sorted result (to make quantiles non-crossing) (cf. Chernozhukove et al. (2010))
	
	# N <- dim(X.min_fill3)[2]
	# T <- dim(X.min_fill3)[1]
	# pred.quant <- array(unlist(pred.quant.imsi), dim=c(T,N,length(tau.vec))) 
	# pred.quant.sort <- aperm(apply(pred.quant, c(1,2), sort), c(2,3,1))	
	
	# # pred.quant.ex : extremal level quantile estimates
		
	 # pred.quant.ex <- array(dim=c(T,N,length(tau.vec.ex)))
	  # for(i in 1:length(tau.vec.ex)){
	    # tau.ex <- tau.vec.ex[i]
	    # pred.quant.ex[,,i] <- (((1-tau.vec[length(tau.vec)])/(1-tau.ex))^(matrix(rep(gamma.hat[clust.loc], T), ncol=N, byrow=TRUE))*pred.quant.sort[,,length(tau.vec)])*(pred.quant.sort[,,length(tau.vec)]>=0) + pred.quant.sort[,,length(tau.vec)]*(pred.quant.sort[,,length(tau.vec)]<0)
	  # }
	
	# tau.vec.bind <- c(tau.vec, tau.vec.ex)
	
	# # idx.index.validation : which index belongs to the current cluster among validation indices? 
	# # For the members of idx.index.validation, inversion of estimated quantile process is executed. 
	# # prediction.forsave is the final prediction object at the current cluster.
	
	# idx.index.invalidation <- NULL
	# for(i in 1:length(clust.loc)){ 
		# idx.index <- which(clust.loc[i]==index.loc.validation)
		# idx.index.invalidation <- c(idx.index.invalidation, idx.index)
		# for(k in idx.index){
			# poo <- index.time.validation[k]
			# for(j in 1:length(tau.vec)){
				# pred.quant.final[k, j] <- pred.quant.sort[poo, i, j]
			# }
			# for(j in 1:length(tau.vec.ex)){
				# pred.quant.final[k, length(tau.vec)+j] <- pred.quant.ex[poo, i, j]
			# }
		# prediction[k,] <- sapply(xk, cdf.est, pred.quant=pred.quant.final[k,], tau.vec=tau.vec.bind)
		# }
	# }
	
	# prediction.forsave <- list(idx.index.invalidation=idx.index.invalidation, prediction.idx = prediction[idx.index.invalidation,])
	
	# filename <- paste("prediction(idx", clust.num, "of68).RDS", sep="")
	# saveRDS(prediction.forsave, file=filename)

	filename.f <- paste("factor(idx", clust.num, "of68).RDS", sep="")
	saveRDS(pred.quant.imsi, file=filename.f)
	
	stopCluster(cl)
}

# Group 5
r <- 5
for(clust.num in cluster.sub5){	# cluster number
	#setup parallel backend to use many processors
	cores=detectCores()
	cl <- makeCluster(cores[1]-1) #not to overload your computer
	registerDoParallel(cl)

	# clust.loc : station number which belogs to clust.num-th cluster
	# X.min_fill3 : (missing-imputed) data object at clust.loc stations
	# temp : PCA result of X.min_fill3; initial value of the iteration
	
	clust.loc <- which(idx==clust.num) 
	X.min_fill3 <- X.min_fill2[,clust.loc] #subset (by clustering)
	temp <- princomp(X.min_fill3) 
	#for(ind in 1:length(tau.vec)){
	pred.quant.imsi <- foreach(ind = 1:length(tau.vec)) %dopar% {
	  library(quantreg)
	  print(paste("###########", ind, "th QFM estimation###########", sep=""))
	  tau <- tau.vec[ind]
	  
	  # To minimize the target function of QFM, we use fully iterative approach. 
	  # L.hat: current loading
	  # L.hat.new : updated loading
	  # F.hat : current factor
	  # F.hat.new : updated factor
	  # maximum iteration number: 200 
	  # If spanned space by updated factors is near to that by previous ones enough, stop the iteration.
	  L.hat <- temp$loadings[,1:r]
	  F.hat <- temp$scores[,1:r]
	  
	  F.hat.new <- matrix(0, nrow=nrow(F.hat), ncol=r)
	  L.hat.new <- matrix(0, nrow=nrow(L.hat), ncol=r)
	  
	  iter <- 1
	  
	  while(1){
	    for(t in 1:nrow(F.hat)){
	      F.hat.new[t,] <- rq(X.min_fill3[t,] ~ L.hat-1, tau=tau)$coeff
	    }
	    for(i in 1:nrow(L.hat)){
	      L.hat.new[i,] <- rq(X.min_fill3[,i] ~ F.hat.new-1, tau=tau)$coeff
	    }

	    if(sum((proj.mat(F.hat.new)-proj.mat(F.hat))^2)<=1e-2){
	      break;
	    }
	    iter <- iter+1
	    F.hat <- F.hat.new
	    L.hat <- L.hat.new
	    if(iter==200){
	      print(paste("Max iteration at tau=", tau))
	      break;
	    }
	  }

	  pred.quant.imsi <- list(F=F.hat.new, L=L.hat.new)
	}
	

	
	
	# # pred.quant : QFM result
	# # pred.quant.sort : sorted result (to make quantiles non-crossing) (cf. Chernozhukove et al. (2010))
	
	# N <- dim(X.min_fill3)[2]
	# T <- dim(X.min_fill3)[1]
	# pred.quant <- array(unlist(pred.quant.imsi), dim=c(T,N,length(tau.vec))) 
	# pred.quant.sort <- aperm(apply(pred.quant, c(1,2), sort), c(2,3,1))	
	
	# # pred.quant.ex : extremal level quantile estimates
		
	 # pred.quant.ex <- array(dim=c(T,N,length(tau.vec.ex)))
	  # for(i in 1:length(tau.vec.ex)){
	    # tau.ex <- tau.vec.ex[i]
	    # pred.quant.ex[,,i] <- (((1-tau.vec[length(tau.vec)])/(1-tau.ex))^(matrix(rep(gamma.hat[clust.loc], T), ncol=N, byrow=TRUE))*pred.quant.sort[,,length(tau.vec)])*(pred.quant.sort[,,length(tau.vec)]>=0) + pred.quant.sort[,,length(tau.vec)]*(pred.quant.sort[,,length(tau.vec)]<0)
	  # }
	
	# tau.vec.bind <- c(tau.vec, tau.vec.ex)
	
	# # idx.index.validation : which index belongs to the current cluster among validation indices? 
	# # For the members of idx.index.validation, inversion of estimated quantile process is executed. 
	# # prediction.forsave is the final prediction object at the current cluster.
	
	# idx.index.invalidation <- NULL
	# for(i in 1:length(clust.loc)){ 
		# idx.index <- which(clust.loc[i]==index.loc.validation)
		# idx.index.invalidation <- c(idx.index.invalidation, idx.index)
		# for(k in idx.index){
			# poo <- index.time.validation[k]
			# for(j in 1:length(tau.vec)){
				# pred.quant.final[k, j] <- pred.quant.sort[poo, i, j]
			# }
			# for(j in 1:length(tau.vec.ex)){
				# pred.quant.final[k, length(tau.vec)+j] <- pred.quant.ex[poo, i, j]
			# }
		# prediction[k,] <- sapply(xk, cdf.est, pred.quant=pred.quant.final[k,], tau.vec=tau.vec.bind)
		# }
	# }
	
	# prediction.forsave <- list(idx.index.invalidation=idx.index.invalidation, prediction.idx = prediction[idx.index.invalidation,])
	
	# filename <- paste("prediction(idx", clust.num, "of68).RDS", sep="")
	# saveRDS(prediction.forsave, file=filename)

	filename.f <- paste("factor(idx", clust.num, "of68).RDS", sep="")
	saveRDS(pred.quant.imsi, file=filename.f)
	
	stopCluster(cl)
}


# Merging the results obtained at each cluster.
# prediction is our final prediction object. 

for(clust.num in 1:68){
	filename <- paste("prediction(idx", clust.num, "of68).RDS", sep="")
	pred.obj <- readRDS(filename)
	index <- pred.obj$idx.index.invalidation
	prediction[index,] <- pred.obj$prediction.idx
}

save(prediction, file="prediction_Multiscale.RData")