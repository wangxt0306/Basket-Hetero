options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
node.idx = as.numeric(args[1]);


if (.Platform$OS.type == "windows") { node.idx = 1 }


if (.Platform$OS.type == "windows") { root.path = "C:/Users/psioda/Documents/Research Papers Temp/bayesDesignBasketHetero/00-github/BHM/SOURCE";          }
if (.Platform$OS.type == "unix")    { root.path = "/proj/psiodalab/projects/bayesDesignBasketHetero/BHM/SOURCE";                                          }

setwd(root.path);
source(".//source-bhm.rcpp");


if (.Platform$OS.type == "windows") { root.path = "C:/Users/psioda/Documents/Research Papers Temp/bayesDesignBasketHetero/00-github/SETTINGS/BINOMIAL";   }
if (.Platform$OS.type == "unix")    { root.path = "/proj/psiodalab/projects/bayesDesignBasketHetero/SETTINGS/BINOMIAL";                                   }

setwd(root.path);
scenarios=read.csv(file="true-parm.csv");
scenarios = rbind(scenarios,scenarios,scenarios,scenarios);


if (.Platform$OS.type == "windows") { root.path = "C:/Users/psioda/Documents/Research Papers Temp/bayesDesignBasketHetero/00-github/BHM/RESULTS_BINOMIAL"; }
if (.Platform$OS.type == "unix")    { root.path = "/proj/psiodalab/projects/bayesDesignBasketHetero/BHM/RESULTS_BINOMIAL";                                 }

setwd(root.path);

set.seed(node.idx);

nSims        =  5000;                                             ## number of simulation studies to perform;
nMC          = 25000;                                             ## number of mcmc samples per analysis;

K0           = 4;                                                 ## number of baskets (i.e., indications)
out.Parms    = matrix(rep(c(4,0.1),K0),nrow=K0,byrow=T);          ## normal distribution parameters for outcome ascertainment times;
data.types   = c("BINARY","BINARY","BINARY","BINARY");            ## data type for each indication
target.n     = c(58,58,58,58);


ppEffCrit = rep(0.975,K0);                                        ## posterior probability critical value (efficacy);
ppFutCrit = rep(0.600,K0);                                        ## posterior probability critical value (futility);



nPerNode = 1;
nrow(scenarios)/nPerNode;

set.start = 1  + nPerNode *(node.idx-1)
set.stop  = min(nPerNode *(node.idx),nrow(scenarios ))


scenarios = scenarios [seq(set.start,set.stop),];
nScenarios = nrow(scenarios );


num.scenarios = nrow(scenarios);

logit = log(0.50/(1-0.50)) - log(0.23/(1-0.23))

offset       = matrix(c(rep(0,K0),rep(logit,K0)),nrow=2,ncol=K0,byrow=T);
	

idx = 0;

for (e in 1:3)          {
for (s in 1:nScenarios) {
for (p in 1:2)          {



	enr.Parms    = c(2.0,2.0,2.0,2.0); 

	r = 1; 
      true.parms = list(K0);
	for (k in 1:K0)
      {
		true.parms[[k]] <- c(unlist(scenarios[s,c(r,r+1)]));
		if (scenarios[s,r]<scenarios[s,r+1] & e==2) { enr.Parms[k] = 4.0; }
		if (scenarios[s,r]>scenarios[s,r+1] & e==3) { enr.Parms[k] = 0.5; }
		r=r+2;
	}

      y = rbind(unlist(true.parms)); colnames(y) <- c("ind1.rr0","ind1.rr1","ind2.rr0","ind2.rr1","ind3.rr0","ind3.rr1","ind4.rr0","ind4.rr1");
	x = unlist(BHM_Design(nSims,nMC,data.types,target.n,ppEffCrit,ppFutCrit,enr.Parms,out.Parms,true.parms,c(offset[p,])));

	z = cbind(e,s,p,y,rbind(x)); colnames(z)[1:3] <- c("enrollment","scenario","pmp");
	idx = idx+1;

	if (idx==1) { all.z = z              }
	else        { all.z = rbind(all.z,z) }
}
}
}

name = paste("./K",K0,"-BHM-",formatC(node.idx, width = 4, format = "d", flag = "0"),".CSV",sep="")
write.csv(all.z,file=name ,row.names=F)

