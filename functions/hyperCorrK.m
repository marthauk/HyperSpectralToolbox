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
number_of_rows = 30;
number_of_cols=30;
if(mod(pixel,number_of_cols)==0)
    pixel_column =number_of_cols;
else
pixel_column=mod(pixel,number_of_cols);
end
lower_limit_matrix = pixel_column - floor(K/2);
%higher_limit_matrix = pixel + floor(K/2);
higher_limit_matrix = lower_limit_matrix + K-1; %because matlab loops include the last index

% Check if index is out of bounds 
% if ( higher_limit_matrix > number_of_cols)
%   % M(band, neighbouring_pixels) * (M(band, Neighbouring pixels)
%   higher_limit_matrix = number_of_cols;
% end
while higher_limit_matrix>number_of_cols
    higher_limit_matrix = higher_limit_matrix-1;
    lower_limit_matrix = lower_limit_matrix-1;
end
% if( lower_limit_matrix < 1)
%     % for edges of the matrix, gonna assume that we just throw out points
%     % outside of the edge, and use half the KERNEL
%     lower_limit_matrix = 1;
%     higher_limit_matrix=K;
% end
while lower_limit_matrix<1
    lower_limit_matrix = lower_limit_matrix +1;
    higher_limit_matrix = higher_limit_matrix +1;
end

% while higher_limit_matrix-lower_limit_matrix>K-1
%     higher_limit_matrix = higher_limit_matrix-1; 
% end
while higher_limit_matrix-lower_limit_matrix<K-1 & higher_limit_matrix >=number_of_cols
    lower_limit_matrix = lower_limit_matrix-1; 
end


lower_limit_row= ceil(pixel/number_of_cols)-floor(K/2);
higher_limit_row = lower_limit_row+floor(K)-1;
while(lower_limit_row<1)
    lower_limit_row = floor(lower_limit_row) +1;
    higher_limit_row = higher_limit_row+1;
end

while(higher_limit_row>number_of_rows)
    higher_limit_row = higher_limit_row-1;
    lower_limit_row = lower_limit_row-1;
end

%take out block of size K x K that I want to work with first
block_K = zeros(p,(higher_limit_matrix-lower_limit_matrix+1)*(higher_limit_matrix-lower_limit_matrix+1)); 

%for n=1:(higher_limit_matrix-lower_limit_matrix+1)
n=1;
for row_pos =lower_limit_row:higher_limit_row
   % block_K(:,1+(n-1)*(higher_limit_matrix-lower_limit_matrix+1):n*(higher_limit_matrix-lower_limit_matrix+1)) = M(:,lower_limit_matrix+(row_pos-1)*number_of_cols:higher_limit_matrix+(row_pos-1)*number_of_cols);
    block_K(:,1+(n-1)*K:n*K) = M(:,lower_limit_matrix+(row_pos-1)*number_of_cols:higher_limit_matrix+(row_pos-1)*number_of_cols);

    n=n+1;
end
%end
%R_k_k = (M(:,lower_limit_matrix:higher_limit_matrix)*M(:,lower_limit_matrix:higher_limit_matrix).')/K;
%remove local mean
%block_K= block_K - repmat(mean(block_K,2),1,(higher_limit_matrix-lower_limit_matrix+1)*(higher_limit_matrix-lower_limit_matrix+1));
%test_r_k=hyperCorr(block_K);

R_k_k =hyperCorr(block_K);
%R_k_k = (block_K *block_K.')/(K*K);

