# Optimum in Convex Hulls (OCH)

OCH is a class of algorithms for extrapolating exclusive RCTs to the broader population by leveraging (partially optimized) physician prescribing patterns seen in inclusive observational data. 

Physicians are like reinforcement learning agents that want to improve patient outcomes over time by giving the best medicines to the appropriate patients. Physicians therefore detect sub-groups of patients who respond well to a treatment, and then give that treatment more often to those patients. In other words, physician prescribing patterns are already partially optimized, and patients usually do much better in observational datasets than in RCTs. OCH formalizes this idea in the potential outcomes framework. 

The OCH algorithms exploit the above concept to extrapolate RCTs. The algorithms can extrapolate the conditinal average treatment effect and the conditional densities of treatment effect to the observational population.

The ``Experiments`` folder contains code to replicate the experimental results for the synthetic data. The real STARD and TEOSS/CATIE datasets used in the paper require approval from the NIMH Data Archive (https://nda.nih.gov/). Please cite the article if you use any of the code in this repository.

# Run the OCH Algorithms

