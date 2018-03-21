clear; clc; 
%addpath('../davidkun-HyperSpectralToolbox-cc014a4/functions')

%generating random image based on cuprite scene data
h=30;
w= 30;
%load('E:\One Drive\OneDrive for Business\NTNU\Master\ground_truthing_aviris_cuprite\cuprite\groundTruth_Cuprite_end12\groundTruth_Cuprite_nEnd12.mat','-mat');
load('groundTruth_Cuprite_nEnd12.mat','-mat');
M_endmembers=M;
goodBands = [10:100 116:150 180:216]; % for AVIRIS with 224 channels
M_endmembers=M(goodBands,:);

[n_bands,k] = size(M_endmembers);
image_30_30 = zeros(h,w,n_bands);

% Setting background
plot(goodBands,M_endmembers(:,:))

for i=1:h
    for j=1:w
            rN = rand;
            image_30_30(i,j,:)= rN * M_endmembers(:,1) + 0.25*M_endmembers(:,3)+0.25* M_endmembers(:,6) +(1-rN)*M_endmembers(:,8);

    end 
end
image=reshape(image_30_30,h*w,n_bands);
figure,plot(goodBands,image(1:10:end,:));
%create kernels with anomalies of size 2x2 with bottom left pixel in 15,15
%column locations

A_pos = 15;
A_size=1;
M_anomaly=M_endmembers(:,4);

for i=A_pos:A_pos+A_size
    for j=A_pos:A_pos+A_size
        image_30_30(i,j,:)=M_anomaly;
    end
end
imnoise(image_30_30,'gaussian',1);
image_noise=reshape(image_30_30,h*w,n_bands);
figure,
subplot(211),plot(goodBands,image(1:10:end,:)),title('No Anomaly');
subplot(212),plot(goodBands,image_noise(1:10:end,:)),title('With Anomaly');

K=10;
r_rlx =hyperLRxDetectorCorr_M(image_30_30,K);
r_rlx_2d = hyperConvert3d(r_rlx.', 30, 30, 1);
figure;imagesc(r_rlx_2d);title(['LRX AD detector, K= ' num2str(K) ]); axis image; colorbar;