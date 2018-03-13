function [R_k_k,block_K] = hyperCorrK(M,K, pixel)
% HYPERCORRK Computes the sample autocorrelation matrix of a kernel K x K
% hyperCorr compute the sample autocorrelation matrix of a 2D matrix.
%
% Usage
%   [R] = hyperCorr(M)
%
% Inputs
%   M - 2D matrix 
% Outputs
%   R - Sample autocorrelation matrix


[p, N] = size(M);

lower_limit_matrix = pixel - floor(K/2);
higher_limit_matrix = pixel + floor(K/2);

% Check if index is out of bounds 
if( pixel-floor(K/2) < 1)
    % for edges of the matrix, gonna assume that we just throw out points
    % outside of the edge, and use half the KERNEL
    lower_limit_matrix = 1;   
end
if ( pixel + floor(K/2) > N)
  % M(band, neighbouring_pixels) * (M(band, Neighbouring pixels)
  higher_limit_matrix = N;
end

%take out block of size K x K that I want to work with first
block_K = zeros(p,K*K); 

number_of_rows = 100;
for n=1:K
    block_K(:,1+(n-1)*K:n*K) = M(:,lower_limit_matrix+(n-1)*number_of_rows:higher_limit_matrix+(n-1)*number_of_rows);
end
%R_k_k = (M(:,lower_limit_matrix:higher_limit_matrix)*M(:,lower_limit_matrix:higher_limit_matrix).')/K;
R_k_k = (block_K *block_K.')/K;
%for n =2:K    
%R_k_k = R_k_k + (M(:,lower_limit_matrix*n:higher_limit_matrix*n)*M(:,lower_limit_matrix*n:higher_limit_matrix*n).')/K;
%end