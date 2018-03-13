function [A_inv,A_mode_elim,A_mode_elim_inv ] = gauss_jordan_inverse(A,mode)
%function [A_inv,A_mode_elim,A_mode_elim_inv ] = gauss_jordan_inverse(A)
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
[size_p,m]=size(A);
A_inv = eye(size_p);
%A_inv = A;
% Forward elimination to build an upper triangular matrix
if(strcmp(mode,'forward') | strcmp(mode,'all'))
 for (i=1:1:size_p)
   if (A(i,i) == 0)
 %i=1; 
  for (j =i+1:1:size_p)
     disp(j);
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
        %A(j,:) = A(j,:)- A(i,:)*A_j_i_temp/A_i_i_temp;
        %A_inv(j,:) = A_inv(j,:) - A_inv(i,:)*A_j_i_temp/A_i_i_temp;
        for (l= 1:size_p)
         A(j,l) = A(j,l)- A(i,l)*A_j_i_temp/A_i_i_temp; 
        A_inv(j,l) = A_inv(j,l) - A_inv(i,l)*A_j_i_temp/A_i_i_temp;
        end
        disp(j);
   end
   disp(A);
end
end

if (strcmp(mode,'forward'))
    A_mode_elim = A;
    A_mode_elim_inv = A_inv;
end

% Backward elimination to build a diagonal matrix
if(strcmp(mode,'backward') | strcmp(mode,'all'))
    for(i=size_p:-1:2)
%i = 3;
   for( j=i-1:-1: 1)
        % The operations below will be different in hardware, because
        % of parallell operations
        A_j_i_temp =A(j,i);
        A_i_i_temp = A(i,i);
        %A(j,:) = A(j,:)-A(i,:)*cast(cast(A(j,i)/A(i,i),'int32'),'double');
        %A_inv(j,:) = A_inv(j,:) - A_inv(i,:)*cast(cast(A_j_i_temp/A_i_i_temp,'int32'),'double');
        
        %A(j,:) = A(j,:)-A(i,:)*A(j,i)/A(i,i);
        %A_inv(j,:) = A_inv(j,:) - A_inv(i,:)*A_j_i_temp/A_i_i_temp;
        
        for (k=1:size_p)
          A(j,k) = A(j,k)-A(i,k)*A(j,i)/A(i,i);
          A_inv(j,k) = A_inv(j,k) - A_inv(i,k)*A_j_i_temp/A_i_i_temp; 
        end
   end
end

end
if (strcmp(mode,'backward'))
A_mode_elim_inv = A_inv;
A_mode_elim = A;
end
if(strcmp(mode,'identity'))
    A_mode_elim_inv = zeros(3);
    A_mode_elim = zeros(3);
end
if(strcmp(mode,'all'))
    A_mode_elim_inv = zeros(3);
    A_mode_elim = zeros(3);
end
% Last division to build an identity matrix
for ( i = 1:+1:size_p)
       disp(A_inv);
    A_inv(i,:)= A_inv(i,:)*1/A(i,i);
 
end





end

