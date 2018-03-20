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
result = zeros(N, 1);

% for all spectral_bands
%for i= 1:p
    % for all N= m x n pixels
    h = waitbar(0,'Initializing waitbar ..');

    for j=1:N
        autocorr = hyperCorrK(M,K,j);
        %autocorrInv = inv(autocorr);
        %disp(autocorrInv);
        %result(j) = M(:,i).' * autocorrInv;
        result(j) = M(:,j).' * inv(autocorr) * M(:,j);
        if(result(j)>=10000)
            weird=1;
        end
   %     temp_result=result(j);
        %result(j)= result(j) *  M(:,i);
      waitbar(j/N,h,'Updated LRX progress');

    end
%end

return;

