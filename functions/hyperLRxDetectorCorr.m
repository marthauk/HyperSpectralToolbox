function [result, autocorr, sigmaInv] = hyperLRxDetectorCorr(M,K)
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

[p, N] = size(M);


% Compute correlation matrix of size K
% correlation matrix will be of size p x p
autocorr_all_pixels = zeros(p*N,1);
result = zeros(N, 1);
anomalies_detected=zeros(p,N/2);
tresh_LRX = 6.0000e+14;

% for all spectral_bands
%for i= 1:p
    % for all N= m x n pixels
    for j=1:N
        autocorr = hyperCorrK(M,K,p);
        autocorrInv = inv(autocorr);
        %disp(autocorrInv);
        %result(j) = M(:,i).' * autocorrInv;
        result(j) = M(:,j).' * autocorrInv * M(:,j);
        
        %result(j)= result(j) *  M(:,i);
    end
%end

result = abs(result);

return;