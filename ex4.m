close all
clear all

load ssm_spins.txt -ascii;
load ssm_spins_test.txt -ascii;

%ssm_kalman(xx, y0, Q0, A, Q, C, R, pass)
% SSM_KALMAN - kalman-smoother estimates of SSM state posterior
%
% [Y,V,Vj,L] = SSM_KALMAN(X,Y0,Q0,A,Q,C,R) peforms the Kalman
% smoothing recursions on the (DxT) data matrix X, for the
% LGSSM defined by the following parameters
%  y0 Kx1 - initial latent state
%  Q0 KxK - initial variance
%  A  KxK - latent dynamics matrix
%  Q  KxK - innovariations covariance matrix
%  C  DxK - output loading matrix
%  R  DxD - output noise matrix

A=0.99*[cos(2*pi/180)  -sin(2*pi/180)     0           0;             ...
        sin(2*pi/180)   cos(2*pi/180)     0           0;             ...
            0                  0      cos(2*pi/90)   -sin(2*pi/90);  ...
            0                  0      sin(2*pi/90)   cos(2*pi/90)]; 

Q=eye(4)-A*A';



C=[ 1   0   1   0; ...
    0   1   0   1; ...
    1   0   0   1; ...
    0   0   1   1; ...
   0.5 0.5 0.5 0.5];

R=eye(5);

Q0=eye(4);
Q=abs(randn(4,4))
Q=Q'*Q;
Q=(Q+Q')/4;
Y0=randn(4,1);
X=ssm_spins';

logdet = @(A)(2*sum(log(diag(chol(A)))));
[Y,V,~,L] = ssm_kalman(X,Y0,Q0,A,Q,C,R, 'filt');
plot(Y');
title('Kalman Filtering - Posterior Mean Estimates')
figure;
plot(cellfun(logdet,V));
title('Kalman Filtering - Log-det of posterior variances')

figure;
[Y,V,Vj,L] = ssm_kalman(X,Y0,Q0,A,Q,C,R, 'smooth');
plot(Y');
title('Kalman Smoothing - Posterior Mean Estimates')
figure;
plot(cellfun(logdet,V));
title('Kalman Smoothing - Log-det of posterior variances')