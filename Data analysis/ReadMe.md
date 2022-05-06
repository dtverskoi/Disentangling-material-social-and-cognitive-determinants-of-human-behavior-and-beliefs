doit_clean.m estimates parameters of the utility function and beliefs dynamics with or without authority messaging.
The functions Init.m (which produces initial conditions for the maximum likelihood procedure) and Collinearity.m 
(which performs multicollinearity detection) are used to produce the above estimates. Note that the function Colli-
nearity.m uses the results of the function vif.m which calculates variance inflation factor (VIF).

The doit_clean.m code also includes:
(1) estimates of the alternative models. To perform this analysis, the function MLEM.m is required;
(2) bootstrap confidence intervals of the parameter estimates;
(3) analysis of gender differences;
(4) effects of the differences in the model parameters between individuals;
(5) analysis of utility losses;
(6) distributions of parameters.
