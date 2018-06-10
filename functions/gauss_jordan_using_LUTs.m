function [A_inv,A_mode_elim,A_mode_elim_inv ] = gauss_jordan_using_LUTs(A,mode)
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

 DIV_PRECISION =17;
 division_lut_values = zeros(1,2^DIV_PRECISION);


for i =1: 2^DIV_PRECISION
    division_lut_values(i)=(2^DIV_PRECISION *1)/i;
end

% divisor_inv set to zero in the beginning
divisor_inv=0;

% Forward elimination to build an upper triangular matrix
if(strcmp(mode,'forward') | strcmp(mode,'all'))
 for (i=1:1:size_p)
   if (A(i,i) == 0)
 %i=1; 
  for (j =i+1:1:size_p)
  %   disp(j);
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
%       error('Matrix is singular'); 
   end
   for (j = i +1:1: size_p)
        % The operations below will be different in hardware, because
        % of parallell operations
        A_j_i_temp =A(j,i);
        A_i_i_temp = A(i,i);
        
        % Multiplying the value up to be able to use it...(Most data in
        % matrix A  = 0.000x
        A_i_i_temp = A_i_i_temp * 1000000;
        % Check if the divisor is negative
        if sign(A_i_i_temp) ==-1 
            A_i_i_temp = A_i_i_temp*-1;
        end
        % look up division
        if A_i_i_temp<= 2^DIV_PRECISION
        A_i_i_temp = uint32(A_i_i_temp);
        divisor_inv = division_lut_values(A_i_i_temp);
        divisor_inv = divisor_inv/1000000;
        else
            divisor_inv=0;
        end
        for (l= 1:size_p)
         A(j,l) = A(j,l)- bitsra(A(i,l)*divisor_inv,DIV_PRECISION); 
         A_inv(j,l) = A_inv(j,l) - bitsra(A_inv(i,l)*A_j_i_temp*divisor_inv,DIV_PRECISION); 
        % A(j,l) = A(j,l)- A(i,l)*A_j_i_temp/A_i_i_temp; 
        %  A_inv(j,l) = A_inv(j,l) - A_inv(i,l)*A_j_i_temp/A_i_i_temp;
        end
    %    disp(j);
   end
 %  disp(A);
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
        A_i_i_temp = A_i_i_temp * 1000000;
        % Check if the divisor is negative
        if sign(A_i_i_temp) ==-1 
            A_i_i_temp = A_i_i_temp*-1;
        end
        % look up division
        if A_i_i_temp<= 2^DIV_PRECISION
        A_i_i_temp = uint32(A_i_i_temp);
        divisor_inv = division_lut_values(A_i_i_temp);
        divisor_inv = divisor_inv/1000000;
        else
            divisor_inv=0;
        end
        
        for (k=1:size_p)
           A(j,k) = A(j,k)- bitsra(A(i,k)*divisor_inv,DIV_PRECISION); 
         A_inv(j,k) = A_inv(j,k) - bitsra(A_inv(i,k)*A_j_i_temp*divisor_inv,DIV_PRECISION); 
        % 
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
     % Multiplying the value up to be able to use it...(Most data in
        % matrix A  = 0.000x
        A_i_i_temp = A_i_i_temp * 1000000;
        % Check if the divisor is negative
        if sign(A_i_i_temp) ==-1 
            A_i_i_temp = A_i_i_temp*-1;
        end
        % look up division
        if A_i_i_temp<= 2^DIV_PRECISION
        A_i_i_temp = uint32(A_i_i_temp);
        divisor_inv = division_lut_values(A_i_i_temp);
        divisor_inv = divisor_inv/1000000;
        else
            divisor_inv=0;
        end
    A_inv(i,:)=bitsra(A(i,i)*divisor_inv,DIV_PRECISION);  
end





end

