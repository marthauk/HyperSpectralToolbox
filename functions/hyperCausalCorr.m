function [R_causal_corr] = hyperCausalCorr(M, pixel_n)
% hyperCausalCorr Computes the sample causal autocorrelation matrix taking
% into account the n-1 pixels preceeding pixel n
%
% Usage
%   [R] = hyperCausalCorr(M, pixel_n)
%
% Inputs
%   M - 2D matrix 
%   pixel_n - nth-pixel vector in the matrix
% Outputs
%   R - Sample causal autocorrelation matrix

%sum_causal = 0;
%for i=1:pixel_n
%    sum_causal = sum_causal + M(:,i)* M(:,i).';
%end
R_causal_corr =  M(:,pixel_n)* M(:,pixel_n).';%1/pixel_n * sum_causal;