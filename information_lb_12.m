function d12=information_lb_12(U,th_points,m0,S0)
% Simplified version of the information lower bound
% for two possible values of the parameter theta. See equations (35), (36)
% in the article.
% To find the optimal signal, minimize this cost function with constraints
% U'*U<=u_max or u_min<=U<=u_max
% This function is used in examples 2 and 3
% The initial conditions m0, S0, are assumed to be independent on theta
% U - the input signal
% th_points a matrix n_{theta} x 2 containing two possible parameters
% m0, S0 - initial mean and covariance matrix of the state vector
% On the basis of the article:
% An Approximate Bayesian Approach to Optimal Input Signal 
% Design for System Identification, by Piotr Bania and Anna Wojcik,
% Entropy, 2025.
% doi: 10.20944/preprints202509.1390.v2
% Copyright CC-BY-NC, 2025, by Piotr Bania (pba@agh.edu.pl)
d12=-dij_calc(U,th_points(:,1),th_points(:,2),m0,S0,m0,S0);
