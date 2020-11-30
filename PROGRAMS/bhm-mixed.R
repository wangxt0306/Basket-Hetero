options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
node.idx = as.numeric(args[1]);


if (.Platform$OS.type == "windows") { node.idx = 1 }


if (.Platform$OS.type == "windows") { root.path = "C:/Users/psioda/Documents/Research Papers Temp/bayesDesignBasketHetero/00-github/BHM/SOURCE";          }
if (.Platform$OS.type == "unix")    { root.path = "/proj/psiodalab/projects/bayesDesignBasketHetero/BHM/SOURCE";                                          }

setwd(root.path);
source(".//source-bhm.rcpp");


if (.Platform$OS.type == "windows") { root.path = "C:/Users/psioda/Documents/Research Papers Temp/bayesDesignBasketHetero/00-github/BHM/RESULTS_MIXED"; }
if (.Platform$OS.type == "unix")    { root.path = "/proj/psiodalab/projects/bayesDesignBasketHetero/BHM/RESULTS_MIXED";                                 }

setwd(root.path);

set.seed(node.idx);

nSims        = 10000;                                             ## number of simulation studies to perform;
nMC          = 20000;                                             ## number of mcmc samples per analysis;

K0           = 4;                                                 ## number of baskets (i.e., indications)
out.Parms    = matrix(rep(c(4,0.1),K0),nrow=K0,byrow=T);          ## normal distribution parameters for outcome ascertainment times;
data.types   = c("BINARY","BINARY","BINARY","NORMAL");            ## data type for each indication


target.n     = c(24,58,34,78);


ppEffCrit = rep(0.975,K0);                                        ## posterior probability critical value (efficacy);
ppFutCrit = rep(0.600,K0);                                        ## posterior probability critical value (futility);


true.parms = list(4);                                             ## model parameters for data generation for each basket in list format;
true.parms[[1]] <- c(0.10,0.47);
true.parms[[2]] <- c(0.23,0.50);
true.parms[[3]] <- c(0.10,0.40);
true.parms[[4]] <- list(c(0.00,0.50),c(1.00,1.00));

all.true.parms = list();


idx = 0;

for (k11 in c(0.10,0.47,0.285)) {
for (k12 in c(0.10,0.47,0.285)) {

for (k21 in c(0.23,0.50,0.365)) {
for (k22 in c(0.23,0.50,0.365)) {

for (k31 in c(0.10,0.40,0.25)) {
for (k32 in c(0.10,0.40,0.25)) {

for (k41 in c(0.00,0.25,0.50)) {
for (k42 in c(0.00,0.25,0.50)) {

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
   tp[[4]][[1]][1] = k41;
   tp[[4]][[1]][2] = k42;


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

offset       = c(2.07708,1.208311,1.791759,0.50);
	

idx = 0;
for (e in 1:1)    {
for (s in 1:nScenarios)    {
for (p in 0:1)    {


      if (p==1) { offset2 = offset*0;} else { offset2 = offset; }

	enr.Parms    = c(2.0,2.0,2.0,2.0); 

      y = rbind(unlist(select.true.parms[[s]])); colnames(y) <- c("ind1.rr0","ind1.rr1","ind2.rr0","ind2.rr1","ind3.rr0","ind3.rr1","ind4.mu0","ind4.mu1","ind4.sd0","ind4.sd1");
	x = unlist(BHM_Design(nSims,nMC,data.types,target.n,ppEffCrit,ppFutCrit,enr.Parms,out.Parms,select.true.parms[[s]],offset2));

	z = cbind(e,s,p,y,rbind(x)); colnames(z)[1:3] <- c("enrollment","scenario","pmp");
	idx = idx+1;

	if (idx==1) { all.z = z              }
	else        { all.z = rbind(all.z,z) }
}
}
}

name = paste("./K",K0,"-BHM-",formatC(node.idx, width = 4, format = "d", flag = "0"),".CSV",sep="")
write.csv(all.z,file=name ,row.names=F)

