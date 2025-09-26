function ilb=information_lb(U,th_points,p_prior,m0,S0)
% Calculation of information lower bound by Lemma 1, 3, 4
% This is the default version which assume that 
% the initial conditions m0, S0 are independent on theta
% If m0, S0 depend on theta write the function get_init_cond 
% for calculation m0(theta), S0(theta) and proceed according to the coments
% in the code below.
% th_points - n_th x r - matrix of the parameters
% p_prior - prior probability distribution, r x 1
% To find the optimal signal minimize this function with 
% siutable constraints for U, e.g. U'*U<=u_max or u_min<=U<=u_max.
% On the basis of the article:
% An Approximate Bayesian Approach to Optimal Input Signal 
% Design for System Identification, by Piotr Bania and Anna Wojcik,
% Entropy, 2025.
% doi: 10.20944/preprints202509.1390.v2
% Copyright CC-BY-NC, 2025, by Piotr Bania (pba@agh.edu.pl)
r=size(th_points,2);D=eye(r);
for i=1:r
    for j=i:r
        if i~=j
            % If the initial conditions depend on theta
            % write the function get_init_cond for calculation
            % of the initial conditions m0(theta), S0(theta)
            %[m01,S01]=get_init_cond(i);[m02,S02]=get_init_cond(j);
            % dij=dij_calc(U,th_points(:,i),th_points(:,j),m01,S01,m02,S02);
            % The default version assume that m0, S0 are fixed
            dij=dij_calc(U,th_points(:,i),th_points(:,j),m0,S0,m0,S0);
            D(i,j)=exp(-dij);
            D(j,i)=D(i,j);
        end
    end
end
% ilb = -p_prior'*log(D*p_prior); % information lower bound by Lemma 1
% equivalent criterion for minimization
ilb=100*(1+(p_prior'*log(D*p_prior))/log(r));
% the minimum value of the above ilb is 0
