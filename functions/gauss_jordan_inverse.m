function [A_inv] = gauss_jordan_inverse(A,size_p)
%GAUSS_JORDAN_INVERSE Summary of this function goes here
% This function implements the Gauss-Jordan method for calculating inverse of a square matrix.
% It acts as a high level model for later implementation in hardware.
%   Detailed explanation goes here
% USAGE:
% Inputs:
% A         -   Matrix of size p x p
% size_p    -   column size
% Outputs:
% A_inv     -   inverse matrix

A_inv = eye(size_p);

% Forward elimination to build an upper triangular matrix
for (i=1:1:size_p)
   if (A(i,i) == 0)
      for (j =i+1:1:size_p)
          if (A(j,j)~=0)
              % The operations below will be different in hardware, because
              % of parallell operations
              %temp_i = row(i);
              %row(i) = row(j);
              %row(j) = temp_i;
              temp_i = A(i,:);
              A(i,:) = A(j,:);
              A(j,:) = temp_i;
          end
      end 
   end
   if (A(i,i) ==0)
      error('Matrix is singular'); 
   end
   for (j = i +1:1: size_p)
        % The operations below will be different in hardware, because
        % of parallell operations
        A_j_i_temp =A(j,i);
        A_i_i_temp = A(i,i);
        A(j,:) = A(j,:)- A(i,:)*A(j,i)/A(i,i);
        A_inv(j,:) = A_inv(j,:) - A_inv(i,:)*A_j_i_temp/A_i_i_temp;
   end
end
%triangular_matrix = A_inv;

% Backward elimination to build a diagonal matrix
for(i=size_p:-1:2)
   for( j=i-1:-1: 1)
        % The operations below will be different in hardware, because
        % of parallell operations
        A_j_i_temp =A(j,i);
        A_i_i_temp = A(i,i);
        A(j,:) = A(j,:)-A(i,:)*A(j,i)/A(i,i);
        A_inv(j,:) = A_inv(j,:) - A_inv(i,:)*A_j_i_temp/A_i_i_temp;
   end
end

% Last division to build an identity matrix
for ( i = 1:+1:size_p)
    A_inv(i,:)= A_inv(i,:)*1/A(i,i);
end



end

