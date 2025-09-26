function [q,log_dets]=get_lfun(th,Y,U,m,S,mp,Sp)
% Calculation of -ln(p(theta|Y)) and ln(det(S(theta,U))) according to formula (59).
% The terms independent on theta are omitted.
% Gaussian prior for theta with the mean mp and covariance Sp.
% To find the MAP estimate of the parameter theta, minimize this function
% by fminbnd for the scalar parameter and by fmincon for vector.
% Input:
% th - n_th x 1 - parameters
% Y - data, size N x n_y 
% U - control signal, size N x n_u
% m, S - initial mean and covariance matrix n x 1, n x n,
% mp, Sp - gaussian prior, mp - n_th x 1, Sp - n_th x n_th
% Output: -ln(p(theta|Y)) and ln(det(S(theta,U)))
% If mp, Sp, are omitted then prior is not taken into account.
% On the basis of the article:
% An Approximate Bayesian Approach to Optimal Input Signal 
% Design for System Identification, by Piotr Bania and Anna Wojcik,
% Entropy, 2025.
% doi: 10.20944/preprints202509.1390.v2
% Copyright CC-BY-NC, 2025, by Piotr Bania (pba@agh.edu.pl)
% This is a version with single output, i.e. y_k\inmathbb{R}, 
% For multiple output change the code according to the comments below.
q=0;log_dets=0;N=size(Y,1); % number of steps
for k=1:N
    [A,B,~,C,D,~,Sv2]=get_abcg(th,U(k,:)'); % matrices A,B,C,G,D,S_v
    SG=Sv2+C*S*C';  % Sigma
    L=S*C'/SG;      % Kalman gain single output
    %L=(SG\(C*S))'; % Kalman gain multiple output
    e=Y(k,:)'-C*m;  % new measurement, prediction error 
    m=m+L*e;        % mean correction
    S=S-(L*SG*L');  % covariance correction
    m=A*m+B;        % mean prediction 
    S=A*S*A'+D;     % covariance prediction   
    lsg=log(SG);    % log det(Sigma), single output
    %lsg=log(det(SG)); % log det(Sigma), multiple output
    en=e*e/SG;         % e'*inv(SG)*e, single output     
    %en=e'*(SG\e);     % e'*inv(SG)*e, multiple output
    al=1-1/k;bt=1-al;
    q=al*q+bt*(en+lsg); % calculation of -ln p(theta|Y)/N
    log_dets=al*log_dets+bt*lsg;    % calculation of log(det(S))/N
end
q=0.5*q;
% Gaussian prior term
if nargin>5
    eth=th-mp;
    q=q+eth*eth/(2*N*Sp);
end
% Output: q=(1/(2N))|Y-F|'*inv(S)*|Y-F|+(1/(2N))log(det(S))-ln(prior)/N
%         log_dets=(1/(2N))log(det(S))

