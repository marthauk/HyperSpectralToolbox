function [R_k_k] = hyperCorrK(M,K, pixel)
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

R_k_k = (M(:,lower_limit_matrix:higher_limit_matrix)*M(:,lower_limit_matrix:higher_limit_matrix).')/K;