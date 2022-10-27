# Optimum in Convex Hulls (OCH)

OCH is a class of algorithms for extrapolating exclusive RCTs to the broader population by leveraging (partially optimized) physician prescribing patterns seen in observational data. 

Physicians are like reinforcement learning agents that want to improve patient outcomes over time by giving the best medicine to the appropriate patient. Physicians therefore detect sub-groups of patients who respond well to a treatment, and then give that treatment more often to those patients. As a result, patients usually do much better in observational datasets than in RCTs. OCH formalizes this idea in the potential outcomes framework.

The algorithm can extrapolate the conditinal average treatment effect and the conditional densities of treatment effect.

The ``Experiments`` folder contains code to replicate the experimental results for the synthetic data. The real STARD and TEOSS/CATIE datasets used in the paper require approval from the NIMH Data Archive (https://nda.nih.gov/). Please cite the article if you use any of the code in this repository.

# Run the OCH Algorithms

