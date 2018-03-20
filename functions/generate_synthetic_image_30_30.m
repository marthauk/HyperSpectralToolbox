%generating random image based on cuprite scene data
h=30;
w= 30;
load('E:\One Drive\OneDrive for Business\NTNU\Master\ground_truthing_aviris_cuprite\cuprite\groundTruth_Cuprite_end12\groundTruth_Cuprite_nEnd12.mat','-mat');
M_endmembers=M;
goodBands = [10:100 116:150 180:216]; % for AVIRIS with 224 channels
M_endmembers=M(goodBands,:);
[n_bands,k] = size(M_endmembers);
image_30_30 = zeros(30,30,n_bands);
reference_anomaly_map = zeros(30,30);
% Setting background
for i=1:h
    for j=1:w
        dice = randi(6);
        if dice>4
            image_30_30(i,j,:)= M_endmembers(:,1); %setting background to alunite
        elseif dice>2 
            image_30_30(i,j,:)= M_endmembers(:,6); %setting background to Kalonite
        else 
            image_30_30(i,j,:)= M_endmembers(:,10);% setting background to pyrope
        end
    end 
end


%create kernels with anomalies of size 1, 5, 10,15, 20, 25 in columns 5, 20,50,100, 400,
%600, in row 35 and 70 
%column locations
% KERNEL_SIZE_ONE_LOCATION =50;
KERNEL_SIZE_TWO_LOCATION = 15;
% KERNEL_SIZE_FIVE_LOCATION =150;
% KERNEL_SIZE_TEN_LOCATION =250;
% KERNEL_SIZE_FIFTEEN_LOCATION =350;
% KERNEL_SIZE_TWENTY_LOCATION =450;
% KERNEL_SIZE_TWENTYFIVE_LOCATION =550;
image_30_30(KERNEL_SIZE_TWO_LOCATION,KERNEL_SIZE_TWO_LOCATION,:)= M_endmembers(:,3);
reference_anomaly_map(KERNEL_SIZE_TWO_LOCATION+1,KERNEL_SIZE_TWO_LOCATION)=1;
reference_anomaly_map(KERNEL_SIZE_TWO_LOCATION,KERNEL_SIZE_TWO_LOCATION)=1;
reference_anomaly_map(KERNEL_SIZE_TWO_LOCATION+1,KERNEL_SIZE_TWO_LOCATION+1)=1;
reference_anomaly_map(KERNEL_SIZE_TWO_LOCATION,KERNEL_SIZE_TWO_LOCATION+1)=1;

image_30_30(KERNEL_SIZE_TWO_LOCATION+1,KERNEL_SIZE_TWO_LOCATION,:)= M_endmembers(:,3);
image_30_30(KERNEL_SIZE_TWO_LOCATION,KERNEL_SIZE_TWO_LOCATION+1,:)= M_endmembers(:,3);
image_30_30(KERNEL_SIZE_TWO_LOCATION+1,KERNEL_SIZE_TWO_LOCATION+1,:)= M_endmembers(:,3);

imnoise(image_30_30,'gaussian',1);
matrix_test=hyperConvert2d(image_30_30);
%[d_acad, anomaly_map,threshold_check_values] = hyperACAD(matrix,100);
K=25;
r_rlx =hyperLRxDetectorCorr(matrix,K);
%d_acad_2d = hyperConvert3d(d_acad.', 30, 30, 1);
r_rlx_2d = hyperConvert3d(r_rlx.', 30, 30, 1);

%anomaly_map_2d = hyperConvert3d(anomaly_map.', 30, 30, 1);
figure;imagesc(r_rlx_2d);title(['LRX AD detector, K= ' num2str(K) ]); axis image; colorbar;
%figure;imagesc(d_acad_2d); axis image; colorbar;

%figure;imagesc(anomaly_map_2d); axis image; colorbar;

%figure;imagesc(reference_anomaly_map);axis image; colorbar;