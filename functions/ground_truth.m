function [ground_truth_map] = ground_truth(height,width, matrix_image, endmembers)
%GROUND_TRUTH Creates a ground truth map for the input image
% INPUT:
% Height: scalar
% width:  scalar
% matrix_image : 2D-matrix 
%endmembers : 2 d matrix of endmembers
ground_truth_map = zeros(height,width);
goodBands = [10:100 116:150 180:216]; % for AVIRIS with 224 channels
[p,N] = size(matrix_image);
endmembers= endmembers(goodBands,:);
[p_end, N_end] = size(endmembers);
NUMROWS=100;
numcolumns= N/NUMROWS;

for i=1:N
    for k=1:N_end
        %if matrix_image(:,j)== endmembers(:,k)
        matrix_vector= matrix_image(:,i);
        endmbembers_vector = endmembers(:,k);
        if isequal(endmbembers_vector, matrix_vector)
            ground_truth_map(ceil(i/numcolumns),mod(i/numcolumns)) = k;
        end
end 
end

