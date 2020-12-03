# Basket-Hetero
Software for paper "Bayesian Adaptive Design for Concurrent Trials Involving Biologically-Related Diseases"

All programs are setup to be executed on a Linux computing cluster using R (3.6.0).The paths referenced in all programs will need to be updated for the code to work. Once all paths are updated, one can use the SLURM scheduler shell scripts to submit jobs on a SLURM-based computing cluster.

----------------------------------------------------
Folder “BHM” contains programs needed (shell script, R, Rcpp) to generate the results in the main paper using Bayesian hierarchical model. 

Subfolders:

CLUSTER-SCRIPT   --- contains shell script used to submit programs on clusters
                     
		     BATCH-BINOMIAL.sh: Shell script to run R program that generates operating characteristics under an identical endpoints case (for example, binary endpoints)
                     
		     BATCH-MIXED.sh: Shell script to run R program that generates operating characteristics under a different endpoints case

SOURCE           --- contains Rcpp program used to generate operating characteristics like type I error rate, power, etc. It takes the inputs from R and also outputs results to R.
                     
		     souce-bhm.rcpp: Generates estimates of power, type I error rate, bias of posterior mean, etc. for both identical and different endpoints cases. 

PROGRAMS         --- contains R code to call the Rcpp programs and save the desired results
                     
		     bhm-binomial.R: Contains inputs corresponding to an identical endpoints case
                     
		     bhm-mixed.R: Contains inputs corresponding to a different endpoints case

RESULTS_BINOMIAL --- contains an example result file of an identical endpoints case

RESULTS_MIXED    --- contains an example result file of a different endpoints case

Inputs for R programs (bhm_binomial.R and bhm_mixed.R have the same types of inputs): 

    nSims          ----  number of simulation studies to perform

    nMC            ----  number of MCMC samples per analysis

    K0             ----   number of baskets (i.e., indications)

    data.types     ----   data type for each indication

    out.Parms      ----   normal distribution parameters for outcome ascertainment times

    target.n       ----   targeted sample size

    ppEffCrit      ----   posterior probability critical value (efficacy)

    ppFutCrit      ----  posterior probability critical value (futility)

    true.parms     ----   model parameters for data generation for each basket in list format


----------------------------------------------------	
	
Folder “BMA” contains programs needed (shell script, R, Rcpp) to generate the results in the main paper using Bayesian model averaging approach. 

Subfolders:

CLUSTER-SCRIPT   --- contains shell script used to submit programs on clusters
                     
		     BATCH-BINOMIAL.sh: Shell script to run R program that generates operating characteristics under an identical endpoints case (for example, binary endpoints)
                     
		     BATCH-MIXED.sh: Shell script to run R program that generates operating characteristics under a different endpoints case

SOURCE           --- contains Rcpp program used to generate operating characteristics like type I error rate, power, etc. It takes the inputs from R and also outputs results to R.
                     
		     souce-bma.rcpp: Generates estimates of power, type I error rate, bias of posterior mean, etc. for both identical and different endpoints cases. 

PROGRAMS         --- contains R code to call the Rcpp programs and save the desired results
                     
		     bma-binomial.R: Contains inputs corresponding to an identical endpoints case
                     
		     bma-mixed.R: Contains inputs corresponding to a different endpoints case

RESULTS_BINOMIAL --- contains an example result file of an identical endpoints case

RESULTS_MIXED    --- contains an example result file of a different endpoints case 


Inputs for R programs (bma_binomial.R and bma_mixed.R have the same types of inputs): 

    nSims           ----  number of simulation studies to perform

    K0              ----   number of baskets (i.e., indications)

    data.types      ----   data type for each indication

    enr.Parms       ----   rate parameter for exponential interarrival times for each basket

    out.Parms       ----   normal distribution parameters for outcome ascertainment times

    target.n        ----   targeted sample size

    ppEffCrit       ----   posterior probability critical value (efficacy)

    ppFutCrit       ----  posterior probability critical value (futility)

    prior.stats     ----  sufficient statistics for conjugate priors for each basket in list format

    a0              ----   value of power in conjugate power prior;

    true.parms      ----   model parameters for data generation for each basket in list format


--------------------------------------------------------

Folder “AppendixA” contains programs needed (shell script, R, Rcpp) to generate the results in Appendix A of the supplementary materials. 

Subfolders:

CLUSTER-SCRIPT --- contains shell script used to submit programs on clusters
                   
		   AppenA.sh: Shell script to run R program that generates operating characteristics under a different endpoints case (for example, binary endpoints)
                                   
SOURCE         --- contains Rcpp program used to generate operating characteristics like type I error rate, power, etc. It takes the inputs from R and also outputs results to R.
                   
		   prior_c.rcpp: Approximates the marginal likelihoods using partition-weighted kernel estimator
                   
		   souce-bma-A.rcpp: Generates estimates of power, type I error rate, bias of posterior mean, etc. a different endpoints case, using the marginal likelihoods approximated from “prior_c.rcpp”

PROGRAMS       --- contains R code to call the Rcpp programs and save the desired results
                   
		   bma-A.R: Contains inputs corresponding to a different endpoints case and outputs results for TableA1-A2

RESULTS        --- contains an example result file of a different endpoints case

Inputs of bma-A.R:

     nSims          ----  number of simulation studies to perform

     K0             ----   number of baskets (i.e., indications)

     covs           ----  whether there is a binary covariate, 0 = no, 1 = yes

     data.types     ----   data type for each indication

     data.links     ----   GLM links 

     output         ----   output risk/mean difference

     enr.Parms      ----   rate parameter for exponential interarrival times for each basket

     out.Parms      ----   normal distribution parameters for outcome ascertainment times

     target.n       ----   targeted sample size

     ppEffCrit      ----   posterior probability critical value (efficacy)

     ppFutCrit      ----  posterior probability critical value (futility)

     prior.stats    ----  sufficient statistics for conjugate priors for each basket in list format

     a0             ----   value of power in conjugate power prior;

     true.parms     ----   model parameters for data generation for each basket in list format

     init.theta     ----   initial values used for MCMC

     piZ            ----   distribution parameters of treatment: Ber(piZ)

     piX            ----   distribution parameters of baseline covariate: Ber(piX)

     lower_limits   ----   lower limits for slice sampling

     upper_limits   ----  upper limits for slice sampling

     slice_widths   ----  widths for slice sampling
