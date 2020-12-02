options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
node.idx = as.numeric(args[1])

########## compile ############

setwd("C:/D/BMA/simulation/cov")
source("data_all_cov.rcpp")
source("prior_c.rcpp")

############ set up #############

nSims        = 500                                           ## number of simulation studies to perform;

K0           = 4                                                ## number of baskets (i.e., indications)
covs        = c(0,1,0,0)                                        ## whether there is a covariate, 0 = no, 1 = yes
data.types   = c("Binomial","Binomial","Binomial","Normal")      ## data type for each indication
data.links   = c("Logistic","Logistic","Logistic","Identity")
output       = c("MD","RD","MD","MD")
enr.Parms    = rep(2.0,K0)                                       ## rate parameter for exponential interarrival times for each basket;
out.Parms    = matrix(rep(c(4,0.1),K0),nrow=K0,byrow=T)          ## normal distribution parameters for outcome ascertainment times;

target.n     = c(24,58,34,80)                                   ## targeted sample size

ppEffCrit = rep(0.975,K0)                                        ## posterior probability critical value (efficacy);
ppFutCrit = rep(0.600,K0)                                        ## posterior probability critical value (futility);

expit = function(x){
  exp(x)/(1+exp(x))
}
prior.stats = list()                                             ## sufficient statistics for conjugate priors for each basket in list format;
prior.stats[[1]] <- c(0.10,0.47)
prior.stats[[2]] <- matrix(c(0.18,expit(-1.5163475+1.2083112),0.28,expit(-1.5163475+1.2083112+0.5718859)),2,2,byrow = T)
## mu00, mu01, mu10, mu11
prior.stats[[3]] <- c(0.10,0.40)
prior.stats[[4]] <- list(c(0.00,0.50), c(1,1))

a0 = rep(1,K0)                                                  ## value of power in conjugate power prior;

logit = function(x){
  log(x/(1-x))
}
true.parms = list()                                             ## model parameters for data generation for each basket in list format;
true.parms[[1]] <- c(logit(0.1),0)
true.parms[[2]] <- c(logit(0.18),0,logit(0.28)-logit(0.18))
true.parms[[3]] <- c(logit(0.1),0)
true.parms[[4]] <- c(0.0,0.0,1.00,1.00)


init.theta = true.parms                                          ## initial values of parameters

piZ = c(0.5,0.5,0.5,0.5)                                         ## distribution of treatment ber(piZ)
piX = list()                                                     ## distribution of baseline covariate ber(piX)
piX[[1]] = rep(0.5,4)
piX[[2]] = rep(0.5,4)                                            ## probability for treated group with x = 1

lower_limits = list()
lower_limits[[1]] = rep(-100,2)
lower_limits[[2]] = rep(-100,3)
lower_limits[[3]] = rep(-100,2)
lower_limits[[4]] = c(rep(-100,2),rep(0,2))
upper_limits = rep( 100,4)
slice_widths = rep(1,4)

########### compute the normalized constant #########

set.seed(node.idx)
#set.seed(11)
const = mcmc_blm(data.types, data.links, K0, target.n, prior.stats, a0, init.theta, piX, piZ,
                 lower_limits, upper_limits, slice_widths, covs, nBI = 1000, nC = 5000)
#const[[3]]

##################################

expit<-function(x){exp(x)/(1+exp(x))}

posterior.b<-function(theta, h, j){
  
  stats = prior.stats[[j]]
  n.t = target.n[j]
  pZ = piZ[j]
  pX = c(piX[[1]][j],piX[[2]][j])
  likeli = 1
  
  if(covs[j] == 1){
    for (z in 0:1){
      hh = z*h
      for (x in 0:1){
        eta = theta%*%c(1,z,x)
        if (data.links[j] == "Logistic"){
          mean = expit(eta)
        }else if (data.links[j] == "Identity"){
          mean = eta
        }
        n = n.t*(z*pZ+(1-z)*(1-pZ))*(x*pX[z+1]+(1-x)*(1-pX[z+1]))
        y = n*stats[x+1,hh+1]
        likeli = likeli*(mean^y*(1-mean)^(n-y))^a0[j] #*(mean^(0.001-1)*(1-mean)^(0.001-1))
      }
    }
  }else{
    for (z in 0:1){
      hh = z*h
      eta = theta[1]+theta[2]*z
      if (data.links[j] == "Logistic"){
        mean = expit(eta)
      }else if (data.links[j] == "Identity"){
        mean = eta
      }
      n = n.t*(z*pZ+(1-z)*(1-pZ))
      y = n*stats[hh+1]
      likeli = likeli*(mean^y*(1-mean)^(n-y))^a0[j]*(mean^(0.001-1)*(1-mean)^(0.001-1))
      likeli = likeli*sqrt(exp(2*theta[1]+theta[2]))/(1+exp(theta[1]))/(1+exp(theta[1]+theta[2])) 
    }
  }
  
  likeli
}

posterior.n<-function(theta, h, j){
  
  stats = prior.stats[[j]]
  n.t = target.n[j]
  pZ = piZ[j]
  pX = c(piX[[1]][j],piX[[2]][j])
  likeli = 1
  
  if (covs[j] == 1){
    for (z in 0:1){
      hh = z*h
      for (x in 0:1){
        eta = theta[1:3]%*%c(1,z,x)
        mean = eta
        n = n.t*(z*pZ+(1-z)*(1-pZ))*(x*pX[z+1]+(1-x)*(1-pX[z+1]))
        y = stats[[1]][x+1,hh+1]
        sd = stats[[2]][x+1,hh+1]
        tau = theta[4+z]
        likeli = likeli*((tau/(2*pi))^(n/2)*exp(-tau*((n-1)*sd^2+n*(y-mean)^2)/2))^a0[j]
        likeli = likeli*sqrt( tau^(0.1-1)*exp(-0.1*tau) )
      }
    }
  }else{
    for (z in 0:1){
      hh = z*h
      eta = theta[1:2]%*%c(1,z)
      mean = eta
      n = n.t*(z*pZ+(1-z)*(1-pZ))
      y = stats[[1]][hh+1]
      sd = stats[[2]][hh+1]
      tau = theta[3+z]
      likeli = likeli*((tau/(2*pi))^(n/2)*exp(-tau*((n-1)*sd^2+n*(y-mean)^2)/2))^a0[j]
      likeli = likeli*tau^(0.1-1)*exp(-0.1*tau)
    }
  }
  
  likeli
}

"%^%" <- function(S, power){
  with(eigen(S), vectors %*% (values^power * t(vectors)))
}

norms = function(t){
  sqrt(sum(t^2))
}

r.p = 0.7
K = 100


c.inv = list()            # list of normalized constants
c.inv[[1]] = rep(0,2)
c.inv[[2]] = rep(0,2)
c.inv[[3]] = rep(0,2)
c.inv[[4]] = rep(0,2)

for (j in 1:K0){
  
  #  if(data.links[j] == "Logistic"){
  
  for (h in 0:1){
    
    theta = const[[2*(j-1)+h+1]]
    p = ncol(theta)
    if (data.types[j] == "Normal"){
      phi = cbind(theta[,1:(p-2)], log(theta[,(p-1):p]))
    }else{
      phi = theta
    }
    
    mean = colMeans(phi)
    sigma = cov(phi)
    psi = sweep(phi,2,mean)%*%(sigma%^%(-0.5))
    
    # construct the working parameter space with radius r #
    psi.norm = apply(psi, 1, norms)
    radius = max(psi.norm)*r.p
    w.psi = index1 = {}
    for (t in 1:nrow(psi)){
      if (psi.norm[t]<radius){
        w.psi = rbind(w.psi,c(psi[t,],psi.norm[t]))
        index1 = c(index1,t)
      }
    }
    
    # partition the working space and compute the kernels needed for c#
    order.psi = w.psi[order(w.psi[,p+1]),]
    index2 = index1[order(w.psi[,p+1])]
    Q.star = Q.t = {}
    for (k in 1:K){
      A.k = index3 = {}
      for (w in 1:nrow(order.psi)){
        if (radius*(k-1)/K <= order.psi[w,p+1] && order.psi[w,p+1] < radius*k/K){
          A.k = rbind(A.k,order.psi[w,])
          index3 = c(index3,index2[w])
        }
      }
      
      if ( is.null(A.k) == TRUE ){
        q.k = 0
        q.star = 0
      }else{
        
        if (data.types[j] == "Normal"){
          sweep.temp = sweep(A.k[,1:p]%*%(sigma%^%(0.5)),2,-mean)
          if (nrow(sweep.temp) == 1){
            A.k.temp = t(as.matrix(c(sweep.temp[,1:(p-2)],exp(sweep.temp[,(p-1):p]))))
          }else{
            A.k.temp = cbind(sweep.temp[,1:(p-2)],exp(sweep.temp[,(p-1):p]))
          }
          Jacob = function(x){
            det(sigma%^%(0.5))*exp(x[p-1]+x[p])
          }
          posterior.temp = rep(NA,nrow(A.k.temp))
          for (t in 1:nrow(A.k.temp)){
            posterior.temp[t] = posterior.n(A.k.temp[t,],h,j)*Jacob(phi[index3[t],])
          }
          q.star = min(posterior.temp)
          q.t = function(q){
            q.star/q
          }
        }else{
          J = det(sigma%^%(0.5))
          # psi.star.temp = psi.star[1:p]%*%(sigma%^%(0.5))+mean
          A.k.temp = sweep(A.k[,1:p]%*%(sigma%^%(0.5)),2,-mean)
          posterior.temp = rep(NA,nrow(A.k.temp))
          for (t in 1:nrow(A.k.temp)){
            posterior.temp[t] = posterior.b(A.k.temp[t,],h,j)*J
          }
          q.star = min(posterior.temp)
          q.t = function(q){
            q.star/q
          }
        }
        
        q.k = sum(apply(as.matrix(posterior.temp), 1, q.t))
        
      }
      
      V.k = ((radius*k/K)^p-(radius*(k-1)/K)^p)*pi^(p/2)/gamma(p/2+1)
      Q.star = c(Q.star,q.star*V.k)
      Q.t = c(Q.t,q.k)
    }
    
    d.hat = sum(Q.t)/(nrow(theta)*sum(Q.star))
    c.inv[[j]][h+1] = log(d.hat)
  }
  
}

#########################################

models = matrix(NA,2^K0,K0)             ## generate model matrix (complete model space unless subset);  
t = c(0,1)
for (k in 0:(K0-1))
{
  models[,k+1] = rep( rep(t, each = 2^(K0-k-1), times = 1 ), times = 2^k, each = 1 )
}

M0         = nrow(models)

num1 = apply((models==1),1,sum)
num0 = apply((models==0),1,sum)
num = cbind(num1,num0)


## enthusiastic independence prior model probabilities;
pmp      = (0.675)^num1*(1-0.675)^num0
pmp      = pmp / sum(pmp)

set.seed(node.idx)
ind = BMA_Design(data.types, data.links, target.n, K0, pmp, prior.stats, a0, true.parms, c.inv, piZ, piX, 
               init.theta, lower_limits, upper_limits, slice_widths, enr.Parms, out.Parms, ppEffCrit, 
               ppFutCrit, nSims, covs, output, nBI = 5000, nMC = 50000)

name = paste("C:/D/BMA/simulation/tableA1/IND/log0-",node.idx,".RData",sep = "")
save(ind, file=name)


## positive dependence model probabilities;
pmp   = (1+abs(num1-num0))^5.430     
pmp   = pmp / sum(pmp)

set.seed(node.idx)
dep = BMA_Design(data.types, data.links, target.n, K0, pmp, prior.stats, a0, true.parms, c.inv, piZ, piX, 
               init.theta, lower_limits, upper_limits, slice_widths, enr.Parms, out.Parms, ppEffCrit, 
               ppFutCrit, nSims, covs, output, nBI = 5000, nMC = 50000)

name = paste("C:/D/BMA/simulation/tableA1/DEP/log0-",node.idx,".RData",sep = "")
save(dep, file=name)


## enthusiastic prior model probabilities \w possible pessimistic nugget;
rho4   = 0.45
rho3   = 0.40
rho0   = 0.005
rhoTot = rho0+rho3+rho4
pmp   = (num1==0)*rho0 + (num1==1)*(1-rhoTot)/(nrow(models)-6) + (num1==2)*(1-rhoTot)/(nrow(models)-6) + (num1==3)*rho3/4 + (num1==4)*rho4

set.seed(node.idx)
alt = BMA_Design(data.types, data.links, target.n, K0, pmp, prior.stats, a0, true.parms, c.inv, piZ, piX, 
               init.theta, lower_limits, upper_limits, slice_widths, enr.Parms, out.Parms, ppEffCrit, 
               ppFutCrit, nSims, covs, output, nBI = 5000, nMC = 50000)

name = paste("C:/D/BMA/simulation/tableA1/ALT/log0-",node.idx,".RData",sep = "")
save(alt, file=name)
