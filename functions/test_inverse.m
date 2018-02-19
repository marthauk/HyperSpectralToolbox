function [inv_A_gauss, inv_A ] = test_inverse()
%TEST_INVERSE Summary of this function goes here
%   Detailed explanation goes here
A=magic(7);
[p,N] = size(A);
inv_A_gauss = gauss_jordan_inverse(A,p);
inv_A= inv(A);

end

