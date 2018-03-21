function [result, autocorr, sigmaInv] = hyperLRxDetectorCorr_M(Mi,K)
%HYPERRX LRX anomaly detector
%   hyperLRxDetector performs the Local RX anomaly detector using Correlation
%   instead of covariance
%
% Usage
%   [result] = hyperRxDetector(M)
% Inputs
%   M  - 2D data matrix (p x N)
%   K  - Size of the kernel window, K x K
% Outputs
%   result - Detector output (1 x N)
%   sigma - Correlation matrix (p x p)
%   sigmaInv - Inverse of correlation matrix (p x p)
[h,w,p]=size(Mi);
M=hyperConvert2d(Mi);
% [p, N] = size(M);
N=h*w;

% Compute correlation matrix of size K
% correlation matrix will be of size p x p
result = zeros(N, 1);


    for j=1:N
        autocorr = hyperCorrK_M(M,h,w,K,j);
        result(j) = M(:,j).' * pinv(autocorr) * M(:,j);
        if(result(j)>=10000)
            weird=1;
        end
    
    end
result = abs(result);

return;