# Optimum in Convex Hulls (OCH)

OCH is a class of algorithms for extrapolating exclusive RCTs to the broader population by leveraging (partially optimized) physician prescribing patterns seen in inclusive observational data. 

Physicians are like reinforcement learning agents that want to improve patient outcomes over time by giving the best medicines to the appropriate patients. Physicians therefore detect sub-groups of patients who respond well to a treatment, and then give that treatment more often to those patients. In other words, physician prescribing patterns are already partially optimized even before an RCT is conducted, and patients usually do much better over time in observational datasets than in RCTs. OCH formalizes this idea in the potential outcomes framework. 

The OCH algorithms exploit the above concept to extrapolate RCTs. The algorithms can extrapolate the conditinal average treatment effect and the conditional densities of treatment effect to the observational population.

The ``Experiments`` folder contains code to replicate the experimental results for the synthetic data. The real STARD and TEOSS/CATIE datasets used in the paper require approval from the NIMH Data Archive (https://nda.nih.gov/). Please cite the article if you use any of the code in this repository.

# Installation

Press the green button up top and download the zip file. Then:

> library(devtools)

> install_local("Directory of OCH-main.zip")

> library(OCH)

# Run the OCH Algorithms

Generate RCT data with 100 samples, observational data with 1000 samples, 5 predictors, and exclude 90% of patients from the RCT:

> data = get_RCT_OBS_data(nR=100, nO=1000, d=5, prop=0.9)

> tR = data$tR; xR = data$xR; yR = data$yR ## treatment, predictors and response for RCT

> tO = data$tO; xO = data$xO; yO = data$yO; mO = data$mO ## treatment, predictors, response and time steps for observational data

Run OCH_2 and predicting on all patients in observational data:

> pred = OCH12(tR,xR,yR,tO,xO,yO,mO,xT=xO)

Also run OCH_1:

> m1 = which(mO==1) ## only use second time step, i.e., after treatment assignment

> pred = OCH12(tR,xR,yR, tO[m1],xO[m1,],yO[m1],mOt,xT=xO)

Also run OCH_d:

> pred = OCHd(tR,xR,yR,tO,xO,yO,mOt,xT=xO)

Can also use OCHd_LSPC function for discrete response.
