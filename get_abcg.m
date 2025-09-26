function [A,B,G,C,D,Sv,Sv2]=get_abcg(th,u)
% System matrices as a function of parameter theta and control u
% Input: 
% th - parameters n_th x 1 
% u -  control vector n_u x 1
% Output
% Matrices A(th,u), B(th,u), G(th,u), C, D=G(th,u)G(th,u)', Sv, Sv2=Sv*Sv'
% On the basis of the article:
% An Approximate Bayesian Approach to Optimal Input Signal 
% Design for System Identification, by Piotr Bania and Anna Wojcik,
% Entropy, 2025.
% doi: 10.20944/preprints202509.1390.v2
% Copyright CC-BY-NC, 2025, by Piotr Bania (pba@agh.edu.pl)
% The matrices below are calculated according to example 3 in the article
T0=0.005747126436782;bc=8.67e+05;Sv=11.85;Sv2=140.4225;
u1=1+u;u1d=u1*T0;e=exp(-u1d);   
thd=th*T0;c=cos(thd);s=sin(thd);
A=e*[c s;-s c];
DEN=u1^2+th^2;
N1=u1+e*(th*s-u1*c);Z1=N1/DEN;
N2=th-e*(u1*s+th*c);Z2=N2/DEN;
B=[Z2;Z1]*bc*u;
e2=e*e;e2_1=1-e2;
D=e2_1*eye(2);
G=sqrt(e2_1)*eye(2);
C=[0 1];

