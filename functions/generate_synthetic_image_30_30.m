clc; clear;close all;
%generating random image based on cuprite scene data
h=30;
w= 30;
%load('E:\One Drive\OneDrive for Business\NTNU\Master\ground_truthing_aviris_cuprite\cuprite\groundTruth_Cuprite_end12\groundTruth_Cuprite_nEnd12.mat','-mat');
load('groundTruth_Cuprite_nEnd12.mat','-mat');
M_endmembers=M;
goodBands = [10:100 116:150 180:216]; % for AVIRIS with 224 channels
M_endmembers=M(goodBands,:);
[n_bands,k] = size(M_endmembers);
image_30_30 = zeros(30,30,n_bands);
reference_anomaly_map = zeros(30,30);
n_true_anomalies =4;
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
%            rN = rand;
%            image_30_30(i,j,:)= rN * M_endmembers(:,1) + 0.25*M_endmembers(:,3)+0.25* M_endmembers(:,6) +(1-rN)*M_endmembers(:,8);

%          rN= rand;
%          image(i,j,:) = rN*M_endmembers(:,1) +0.2*M_endmembers(:,3)+0.2*M_endmembers(:,4)+0.2*M_endmembers(:,7)+rN*M_endmembers(:,12);

    end 
end


%create kernels with anomalies of size 2x2 with bottom left pixel in 15,15
%column locations

KERNEL_SIZE_TWO_LOCATION = 15;

image_30_30(KERNEL_SIZE_TWO_LOCATION,KERNEL_SIZE_TWO_LOCATION,:)= M_endmembers(:,3);
reference_anomaly_map(KERNEL_SIZE_TWO_LOCATION+1,KERNEL_SIZE_TWO_LOCATION)=1;
reference_anomaly_map(KERNEL_SIZE_TWO_LOCATION,KERNEL_SIZE_TWO_LOCATION)=1;
reference_anomaly_map(KERNEL_SIZE_TWO_LOCATION+1,KERNEL_SIZE_TWO_LOCATION+1)=1;
reference_anomaly_map(KERNEL_SIZE_TWO_LOCATION,KERNEL_SIZE_TWO_LOCATION+1)=1;

image_30_30(KERNEL_SIZE_TWO_LOCATION+1,KERNEL_SIZE_TWO_LOCATION,:)= M_endmembers(:,3);
image_30_30(KERNEL_SIZE_TWO_LOCATION,KERNEL_SIZE_TWO_LOCATION+1,:)= M_endmembers(:,3);
image_30_30(KERNEL_SIZE_TWO_LOCATION+1,KERNEL_SIZE_TWO_LOCATION+1,:)= M_endmembers(:,3);

imnoise(image_30_30,'gaussian',1);
matrix=hyperConvert2d(image_30_30);
%[d_acad, anomaly_map,threshold_check_values] = hyperACAD(matrix,100);
% K is size of kernel
K=5;
treshold = 0.9;
%[r_alrx,anomaly_map,not_used ,not_use] =hyperLRX_anomaly_set_remover(matrix,K,treshold);
[r_alrx,anomaly_map,not_used ] =hyperACAD(matrix,treshold);

%d_acad_2d = hyperConvert3d(d_acad.', 30, 30, 1);
r_alrx_2d = hyperConvert3d(r_alrx.', 30, 30, 1);
anomaly_map_2d = hyperConvert3d(anomaly_map.', 30, 30, 1);

%figure;imagesc(r_alrx_2d);title(['ALRX AD detector, K= ' num2str(K) ]); axis image; colorbar;
figure;imagesc(r_alrx_2d);title(['ACAD result, treshold' num2str(treshold) ]); axis image; colorbar;

figure;imagesc(anomaly_map_2d);title(['ACAD anomaly map, treshold= ' num2str(treshold) ]); axis image; colorbar;

%% Evaluate the performance of the AD by setting objective measures
% find max value outputted from the AD
%max_ad_score = max(r_rlx);
treshold_percentage = 0.75;
predicted_anomalies =0;
% for i=1:w*h 
%     if r_rlx(i)>=treshold_percentage *max_ad_score
%         predicted_anomalies =predicted_anomalies+1;
%     end
% end
n_actual_anomalies=n_true_anomalies;
n_true_anomalies =0;
 for i=1:w
     for j=1:h
    if anomaly_map_2d(i,j)==1
        predicted_anomalies =predicted_anomalies+1;
        if reference_anomaly_map(i,j)==1
         n_true_anomalies=n_true_anomalies+1;
        end
    end
    end
end

false_anomalies = predicted_anomalies-n_true_anomalies;
if predicted_anomalies<n_actual_anomalies
    correctly_predicted_anomalies =n_true_anomalies/n_actual_anomalies;
else
correctly_predicted_anomalies =n_true_anomalies/predicted_anomalies;
end


%figure;imagesc(d_acad_2d); axis image; colorbar;

%figure;imagesc(anomaly_map_2d); axis image; colorbar;

%figure;imagesc(reference_anomaly_map);axis image; colorbar;
