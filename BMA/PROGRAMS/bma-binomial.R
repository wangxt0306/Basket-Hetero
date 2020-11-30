options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
node.idx = as.numeric(args[1]);


## set note if running on windows platform;
if (.Platform$OS.type == "windows") { node.idx = 1 }

if (.Platform$OS.type == "windows") { root.path = "C:/Users/psioda/Desktop/bayesDesignBasketHetero/BMA/SOURCE";          }
if (.Platform$OS.type == "unix")    { root.path = "/proj/psiodalab/projects/bayesDesignBasketHetero/BMA/SOURCE";                                          }

setwd(root.path);
source(".//source-bma.rcpp");

if (.Platform$OS.type == "windows") { root.path = "C:/Users/psioda/Desktop/bayesDesignBasketHetero/BMA/RESULTS_BINOMIAL"; }
if (.Platform$OS.type == "unix")    { root.path = "/proj/psiodalab/projects/bayesDesignBasketHetero/BMA/RESULTS_BINOMIAL";                                 }

setwd(root.path);


set.seed(node.idx);

nSims        = 100000;                                            ## number of simulation studies to perform;

K0           = 4;                                                 ## number of baskets (i.e., indications)
data.types   = c("BINARY","BINARY","BINARY","BINARY");            ## data type for each indication
enr.Parms    = c(2.0,2.0,2.0,2.0);                                ## rate parameter for exponential interarrival times for each basket;
out.Parms    = matrix(rep(c(4,0.1),K0),nrow=K0,byrow=T);          ## normal distribution parameters for outcome ascertainment times;


target.n     = c(58,58,58,58);


ppEffCrit = rep(0.975,K0);                                        ## posterior probability critical value (efficacy);
ppFutCrit = rep(0.600,K0);                                        ## posterior probability critical value (futility);

rr0 = 0.23;
rr1 = 0.50;

prior.stats = list(4);                                            ## sufficient statistics for conjugate priors for each basket in list format;
prior.stats[[1]] <- c(rr0,rr1);                    
prior.stats[[2]] <- c(rr0,rr1); 
prior.stats[[3]] <- c(rr0,rr1); 
prior.stats[[4]] <- c(rr0,rr1); 
a0 = 1.000;                                                       ## value of power in conjugate power prior;

true.parms = prior.stats;

all.true.parms = list();


idx = 0;

for (k11 in c(0.23,0.50,0.365)) {
for (k12 in c(0.23,0.50,0.365)) {

for (k21 in c(0.23,0.50,0.365)) {
for (k22 in c(0.23,0.50,0.365)) {

for (k31 in c(0.23,0.50,0.365)) {
for (k32 in c(0.23,0.50,0.365)) {

for (k41 in c(0.23,0.50,0.365)) {
for (k42 in c(0.23,0.50,0.365)) {

 if (k12>=k11 & k22>=k21 & k32>=k31 & k42>=k41)
  {

   tp = list();
   tp[[1]] = true.parms[[1]];
   tp[[1]][1] = k11;
   tp[[1]][2] = k12;

   tp[[2]] = true.parms[[2]];
   tp[[2]][1] = k21;
   tp[[2]][2] = k22;

   tp[[3]] = true.parms[[3]];
   tp[[3]][1] = k31;
   tp[[3]][2] = k32;

   tp[[4]] = true.parms[[4]];
   tp[[4]][1] = k41;
   tp[[4]][2] = k42;

   idx = idx + 1;
   all.true.parms[[idx]] <- tp
  }

}}}}}}}}

numScenarios =length(all.true.parms);

nPerNode = 2;
numScenarios /nPerNode;

set.start = 1  + nPerNode *(node.idx-1)
set.stop  = min(nPerNode *(node.idx),numScenarios )


select.true.parms = list();
idx=1;
for (s in seq(set.start,set.stop))
{
	select.true.parms[[idx]] <- all.true.parms[[s]];
	idx=idx+1;

}
rm(all.true.parms);
nScenarios =length(select.true.parms);



models     = model.matix(K0)                                      ## generate model matrix (complete model space unless subset);
M0         = nrow(models);


num1 = apply((models==1),1,sum);
num0 = apply((models==0),1,sum);
num = cbind(num1,num0);

pmp.set = list(6);

## enthusiastic independence prior model probabilities;
pmp      = (0.675)^num1*(1-0.675)^num0;	
pmp      = pmp / sum(pmp); round(pmp,3)
marg.pos = sum( (models[,1]==1)*pmp);marg.pos;  pmp[length(pmp)];
pmp.set[[1]] <-  pmp;

## positive dependence model probabilities;
pmp   = (1+abs(num1-num0))^5.430;	
pmp   = pmp / sum(pmp); round(pmp,3)
marg.pos = sum( (models[,1]==1)*pmp);marg.pos; pmp[length(pmp)];
pmp.set[[2]] <-  pmp;

## positive dependence model probabilities;
rho4   = 0.45;
rho3   = 0.40;
rho0   = 0.005;
rhoTot = rho0+rho3+rho4;
pmp   = c(rho0, rep((1-rhoTot)/(nrow(models)-6),nrow(models)-6),rep(rho3/4,4),rho4);	 sum(pmp)
marg = sum(pmp*c(models[,1])); marg ;
pmp.set[[3]] <-  pmp; 


pmp.set[[4]] <- rep(1/2^K0,2^K0);
pmp.set[[5]] <- rep(1/2^K0,2^K0);
pmp.set[[6]] <- rep(1/2^K0,2^K0);


n.pmp = length(pmp.set)

idx = 0;

for (e in 1:1)          {
for (s in 1:nScenarios) {
for (p in 1:n.pmp)      {

	ppEffCrit = rep(0.975,K0);		## posterior probability critical value (efficacy);
	ppFutCrit = rep(0.600,K0); 		## posterior probability critical value (futility);

	enr.Parms    = c(2.0,2.0,2.0,2.0); 
	m   = models;
	a   = a0;
	pmp = pmp.set[[p]];
	
	
	if (p==5)
	{
		a   = 0.005;
		m   = matrix(models[1,],nrow=1,byrow=T);
		pmp = 1;
	}
	if (p==6)
	{
		a   = 0.005;
		m   = matrix(models[1,],nrow=1,byrow=T);
		pmp = 1;
		ppEffCrit = rep(0.91,K0);                                 ## posterior probability critical value (efficacy);
	}
	
	y = rbind(unlist(select.true.parms[[s]])); colnames(y) <- c("ind1.rr0","ind1.rr1","ind2.rr0","ind2.rr1","ind3.rr0","ind3.rr1","ind4.rr0","ind4.rr1");
	x = unlist(BMA_Design(nSims,data.types,target.n,ppEffCrit,ppFutCrit,enr.Parms,out.Parms,m,pmp,prior.stats,a,select.true.parms[[s]]));

	z = cbind(e,s,p,y,rbind(x)); colnames(z)[1:3] <- c("enrollment","scenario","pmp");
	idx = idx+1;

	if (idx==1) { all.z = z          }
	else        { all.z = rbind(all.z,z) }

}
}
}


name = paste("./K",K0,"-BMA-",formatC(node.idx, width = 4, format = "d", flag = "0"),".CSV",sep="")
write.csv(all.z,file=name ,row.names=F)



