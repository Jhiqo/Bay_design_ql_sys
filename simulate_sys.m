function [Y,X]=simulate_sys(th,U,xini)
% System simulator
% Input:
% th - n_th x 1 - parameters vector 
% U - control, input signal, size N x n_u
% xini - initial conditions
% Output:
% X - state trajectory, size N x n_x
% Y - observations, size N x n_y
% On the basis of the article:
% An Approximate Bayesian Approach to Optimal Input Signal 
% Design for System Identification, by Piotr Bania and Anna Wojcik,
% Entropy, 2025.
% doi: 10.20944/preprints202509.1390.v2
% Copyright CC-BY-NC, 2025, by Piotr Bania (pba@agh.edu.pl)
 
 [N,n_u]=size(U); % number of steps, number of controls
 [A,~,G,C,~,Sv,~]=get_abcg(th,zeros(n_u,1)); % matrices A,B,C,G,D,S_v
 nx=size(A,1);ny=size(C,1); nw=size(G,2);
 X=zeros(N,nx);Y=zeros(N,ny);x=xini;
for k=1:N
    [A,B,G,~,~,~]=get_abcg(th,U(k,:)'); % matrices A,B,C,G,D,S_v
    X(k,:)=x';
    Y(k,:)=C*x+Sv*randn(ny,1);
    x=A*x+B+G*randn(nw,1);        
end
% Output: Y, X

