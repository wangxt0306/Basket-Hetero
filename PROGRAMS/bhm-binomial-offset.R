options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
node.idx = as.numeric(args[1]);


if (.Platform$OS.type == "windows") { node.idx = 1 }


if (.Platform$OS.type == "windows") { root.path = "C:/Users/psioda/Desktop/bayesDesignBasketHetero/BHM/SOURCE";          }
if (.Platform$OS.type == "unix")    { root.path = "/proj/psiodalab/projects/bayesDesignBasketHetero/BHM/SOURCE";                                          }

setwd(root.path);
source(".//source-bhm.rcpp");


if (.Platform$OS.type == "windows") { root.path = "C:/Users/psioda/Desktop/bayesDesignBasketHetero/BHM/RESULTS_BINOMIAL"; }
if (.Platform$OS.type == "unix")    { root.path = "/proj/psiodalab/projects/bayesDesignBasketHetero/BHM/RESULTS_BINOMIAL";                                 }

setwd(root.path);

set.seed(node.idx);

nSims        =   500;                                             ## number of simulation studies to perform;
nMC          = 50000;                                             ## number of mcmc samples per analysis;

K0           = 4;                                                 ## number of baskets (i.e., indications)
out.Parms    = matrix(rep(c(4,0.1),K0),nrow=K0,byrow=T);          ## normal distribution parameters for outcome ascertainment times;
data.types   = c("BINARY","BINARY","BINARY","BINARY");            ## data type for each indication
target.n     = c(58,58,58,58);


ppEffCrit = rep(0.975,K0);                                        ## posterior probability critical value (efficacy);
ppFutCrit = rep(0.600,K0);                                        ## posterior probability critical value (futility);


rr0 = 0.23;
rr1 = 0.50;

true.parms = list(4);                                             ## sufficient statistics for conjugate priors for each basket in list format;
true.parms[[1]] <- c(rr0,rr1);                    
true.parms[[2]] <- c(rr0,rr1); 
true.parms[[3]] <- c(rr0,rr1); 
true.parms[[4]] <- c(rr0,rr1); 


all.true.parms = list();


idx = 0;

for (set in c(1,2,3,4,5))     {
k11 = c(0.23,0.23,0.23,0.23,0.23)  
k12 = c(0.23,0.50,0.50,0.50,0.50) 

k21 = c(0.23,0.23,0.23,0.23,0.23) 
k22 = c(0.23,0.23,0.50,0.50,0.50)  

k31 = c(0.23,0.23,0.23,0.23,0.23)  
k32 = c(0.23,0.23,0.23,0.50,0.50)  

k41 = c(0.23,0.23,0.23,0.23,0.23) 
k42 = c(0.23,0.23,0.23,0.23,0.50)  

   tp = list();
   tp[[1]] = true.parms[[1]];
   tp[[1]][1] = k11[set];
   tp[[1]][2] = k12[set];

   tp[[2]] = true.parms[[2]];
   tp[[2]][1] = k21[set];
   tp[[2]][2] = k22[set];

   tp[[3]] = true.parms[[3]];
   tp[[3]][1] = k31[set];
   tp[[3]][2] = k32[set];

   tp[[4]] = true.parms[[4]];
   tp[[4]][1] = k41[set];
   tp[[4]][2] = k42[set];


   idx = idx + 1;
   all.true.parms[[idx]] <- tp

}

numScenarios =length(all.true.parms);

nPerNode = 5;
numScenarios /nPerNode;

set.start = 1  
set.stop  = numScenarios


select.true.parms = list();
idx=1;
for (s in seq(set.start,set.stop))
{
	select.true.parms[[idx]] <- all.true.parms[[s]];
	idx=idx+1;

}
rm(all.true.parms);
nScenarios =length(select.true.parms);



offset1       = c(0,0,0,0);
offset2       = c(1.208311,1.208311,1.208311,1.208311);
	

idx = 0;

for (e in 1:2)          {
for (s in 1:nScenarios) {

	if (e==1) offset=offset1;
	if (e==2) offset=offset2;

      p = 0;

	enr.Parms    = c(2.0,2.0,2.0,2.0); 


        y = rbind(unlist(select.true.parms[[s]])); colnames(y) <- c("ind1.rr0","ind1.rr1","ind2.rr0","ind2.rr1","ind3.rr0","ind3.rr1","ind4.rr0","ind4.rr1");
	x = unlist(BHM_Design(nSims,nMC,data.types,target.n,ppEffCrit,ppFutCrit,enr.Parms,out.Parms,select.true.parms[[s]],offset));

	z = cbind(e,s,p,y,rbind(x)); colnames(z)[1:3] <- c("enrollment","scenario","pmp");
	idx = idx+1;

	if (idx==1) { all.z = z              }
	else        { all.z = rbind(all.z,z) }
}
}

name = paste("./K",K0,"-OFFSET-",formatC(node.idx, width = 4, format = "d", flag = "0"),".CSV",sep="")
write.csv(all.z,file=name ,row.names=F)
