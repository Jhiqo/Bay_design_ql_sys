function dij=dij_calc(U,th_1,th_2,m01,S01,m02,S02)
% Calculation of d_ij quantity according to Lemma 3 and 4
% U - input signal, design parameter,
% th_1, th_2 - parameters vectors
% m01,S01, - mean and covariance for th_1
% m02,S02, - mean and covariance for th_2
% On the basis of the article:
% An Approximate Bayesian Approach to Optimal Input Signal 
% Design for System Identification, by Piotr Bania and Anna Wojcik,
% Entropy, 2025.
% doi: 10.20944/preprints202509.1390.v2
% Copyright CC-BY-NC, 2025, by Piotr Bania (pba@agh.edu.pl)
nx=size(m01,1);nx1=nx+1;nx2=2*nx;
N=size(U,1); % number of steps
[~,~,~,C1,~,~,Sv2]=get_abcg(th_1,U(1,:)'); % matrices A,B,C,G,D,S_v
C2=C1;C=[C1 -C1]/sqrt(2);A=zeros(nx2,nx2);D=zeros(nx2,nx2);S=zeros(nx2,nx2);
m=[m01;m02];S(1:nx,1:nx)=S01;S(nx1:nx2,nx1:nx2)=S02;S1=S01;S2=S02;
q=0;log_detr=0;log_detsg1=0;log_detsg2=0;

for k=1:N
    [A1,B1,~,~,D1,~,~]=get_abcg(th_1,U(k,:)'); % matrices A,B,C,G,D,S_v
    [A2,B2,~,~,D2,~,~]=get_abcg(th_2,U(k,:)'); % matrices A,B,C,G,D,S_v
    A(1:nx,1:nx)=A1;A(nx1:nx2,nx1:nx2)=A2; % \tilde{A} matrix of the extended system
    D(1:nx,1:nx)=D1;D(nx1:nx2,nx1:nx2)=D2; % diffusion matrix of the extended system
    B=[B1;B2];
    R=Sv2+C*S*C';   % Sigma_k
    lr=log(R);    % log det(Sigma_k), single output
    %lr=log(det(R)); % log det(Sigma_k), multiple output
    L=S*C'/R;      % Kalman gain single output
    %L=(R\(C*S))'; % Kalman gain multiple output
    S=S-(L*R*L');  % covariance correction
    cm=C*m;                        % C*m_k 
    en=cm*cm/R;                    % m_k'*C'*inv(Sigma_k)*C*m_k, single output     
    %en=cm'*(R\cm);                % m_k'*C'*inv(sigma_k)*C*m_k, multiple output
    m=A*(eye(nx2)-L*C)*m+B;        % mean prediction 
    S=A*S*A'+D;                    % covariance prediction
    q=q+en;                   % 
    log_detr=log_detr+lr;     %
    %********* det(S1) ********** 
    SG1=Sv2+C1*S1*C1';  % Sigma_1
    L1=S1*C1'/SG1;      % Kalman gain single output
    %L1=(SG1\(C1*S1))'; % Kalman gain multiple output
    S1=S1-(L1*SG1*L1');  % covariance correction 
    S1=A1*S1*A1'+D1;     % covariance prediction
    log_detsg1=log_detsg1+log(SG1); % log(det(S1)) single output
    %log_detsg1=log_detsg1+log(det(SG1)); % log(det(S1)) multiple output
    
    %********* det(S2) ********** 
    SG2=Sv2+C2*S2*C2';  % Sigma_2
    L2=S2*C2'/SG2;      % Kalman gain single output
    %L2=(SG2\(C1*S2))'; % Kalman gain multiple output
    S2=S2-(L2*SG2*L2');  % covariance correction 
    S2=A2*S2*A2'+D2;     % covariance prediction
    log_detsg2=log_detsg2+log(SG2); % log(det(S2)) single output
    %log_detsg2=log_detsg2+log(det(SG2)); % log(det(S2)) multiple output
end
% Formula (74)
dij=0.25*q+0.5*log_detr-0.25*(log_detsg1+log_detsg2);
