%%
clear; clc;close all;
load('groundTruth_Cuprite_nEnd12.mat','-mat');
w = 200;
h = 200;

resultsDir= ['M:\Documents\Forked_MATLAB_hyperspectral_toolbox\HyperSpectralToolbox\figures\Hsueh' ,datestr(now, 'dd-mmm-yyyy')];

%resultsDir =regexprep(resultsDir,':d*','')
%resultsDir ='E:\One Drive\OneDrive for Business\NTNU\Master\Anomaly detection results\MATLAB\synthetic_images\lol' ;
[status, msg, msgID] = mkdir(resultsDir);

M_endmembers=M;
goodBands = [10:100 116:150 180:216]; % for AVIRIS with 224 channels
M_endmembers=M(goodBands,:);
[n_bands,k] = size(M_endmembers);
image = zeros(h,w,n_bands);
reference_anomaly_map = zeros(h,w);
BACKGROUND=0.2*M_endmembers(:,1) +0.2*M_endmembers(:,3)+0.2*M_endmembers(:,5)+0.2*M_endmembers(:,7)+0.2*M_endmembers(:,12);
SNR = 20;
% Setting background
for i=1:h
    for j=1:w          
        image(i,j,:) = BACKGROUND;
    end 
end

%column locations
KERNEL_SIZE_TWO_LOCATION =40; % column one
KERNEL_SIZE_TWO_2_LOCATION = 70; %column two, mixed pixels
KERNEL_SIZE_TWO_MIXED_LOCATION =100; % column three, mixed pixel
KERNEL_SIZE_ONE_BKG_MIXED_LOCATION =130; % column four, mixed pixel and background, 50/50
KERNEL_SIZE_ONE_BKG_75_MIX_LOCATION =160; %column five, mixed pixel and background, 25/75
M_anomaly_pure = M_endmembers(:,10);
M_A = M_endmembers(:,1);
M_B = M_endmembers(:,3);
M_K = M_endmembers(:,5);
M_M = M_endmembers(:,7);
M_C = M_endmembers(:,12);

for i=KERNEL_SIZE_TWO_LOCATION:30:KERNEL_SIZE_ONE_BKG_75_MIX_LOCATION
    %first column
    image(i,KERNEL_SIZE_TWO_LOCATION,:) = M_anomaly_pure;
    reference_anomaly_map(i,KERNEL_SIZE_TWO_LOCATION)=1;

    image(i,KERNEL_SIZE_TWO_LOCATION +1,:) = M_anomaly_pure;
    reference_anomaly_map(i,KERNEL_SIZE_TWO_LOCATION +1)=1;

    image(i+1,KERNEL_SIZE_TWO_LOCATION,:) = M_anomaly_pure;
    reference_anomaly_map(i+1,KERNEL_SIZE_TWO_LOCATION)=1;

    image(i+1,KERNEL_SIZE_TWO_LOCATION +1,:) = M_anomaly_pure;
    reference_anomaly_map(i+1,KERNEL_SIZE_TWO_LOCATION +1)=1;
    
    % second column
     image(i,KERNEL_SIZE_TWO_2_LOCATION,:) = M_anomaly_pure;
    reference_anomaly_map(i,KERNEL_SIZE_TWO_2_LOCATION)=1;

    image(i,KERNEL_SIZE_TWO_2_LOCATION +1,:) = M_anomaly_pure;
    reference_anomaly_map(i,KERNEL_SIZE_TWO_2_LOCATION +1)=1;

    image(i+1,KERNEL_SIZE_TWO_2_LOCATION,:) = M_anomaly_pure;
    reference_anomaly_map(i+1,KERNEL_SIZE_TWO_2_LOCATION)=1;

    image(i+1,KERNEL_SIZE_TWO_2_LOCATION +1,:) = M_anomaly_pure;
    reference_anomaly_map(i+1,KERNEL_SIZE_TWO_2_LOCATION +1)=1;
        
    for(j=KERNEL_SIZE_TWO_MIXED_LOCATION:30:KERNEL_SIZE_ONE_BKG_75_MIX_LOCATION)
        switch j
              case KERNEL_SIZE_TWO_MIXED_LOCATION
                  switch i
                      case 40
                          image(i+1,j,:)= 0.5*M_A + 0.5 *M_B;
                          image(i+1,j+1,:)=0.5*M_A + 0.5*M_C;
                          image(i,j,:) = 0.5*M_A + 0.5*M_K;
                          image(i,j+1,:)= 0.5*M_A + 0.5*M_M;
                      case 70
                          image(i+1,j,:)= 0.5*M_A + 0.5 *M_B;
                          image(i+1,j+1,:)=0.5*M_A + 0.5*M_C;
                          image(i,j,:) = 0.5*M_B + 0.5*M_K;
                          image(i,j+1,:)= 0.5*M_B + 0.5*M_M;
                      case 100
                          image(i+1,j,:)= 0.5*M_A + 0.5 *M_C;
                          image(i+1,j+1,:)=0.5*M_B + 0.5*M_C;
                          image(i,j,:) = 0.5*M_C + 0.5*M_K;
                          image(i,j+1,:)= 0.5*M_C + 0.5*M_M;
                      case 130
                          image(i+1,j,:)= 0.5*M_A + 0.5 *M_K;
                          image(i+1,j+1,:)=0.5*M_B +0.5*M_K;
                          image(i,j,:) = 0.5*M_C + 0.5*M_K;
                          image(i,j+1,:)= 0.5*M_K + 0.5*M_M;
                     case 160
                          image(i+1,j,:)= 0.5*M_A + 0.5 *M_M;
                          image(i+1,j+1,:)=0.5*M_B + 0.5*M_M;
                          image(i,j,:) = 0.5*M_C + 0.5*M_M;
                          image(i,j+1,:)= 0.5*M_K + 0.5*M_M;
                  end
                  reference_anomaly_map(i,j)=1;
                  reference_anomaly_map(i+1,j)=1;
                  reference_anomaly_map(i+1,j+1)=1;
                  reference_anomaly_map(i,j+1)=1;
            case  KERNEL_SIZE_ONE_BKG_MIXED_LOCATION
                switch i
                    case 40
                        image(i,j,:) = 0.5*M_A + 0.5*BACKGROUND;
                         reference_anomaly_map(i,j)=1;
                    case 70
                        image(i,j,:) = 0.5*M_B + 0.5*BACKGROUND;
                         reference_anomaly_map(i,j)=1;
                    case 100
                        image(i,j,:) = 0.5*M_C + 0.5*BACKGROUND;
                         reference_anomaly_map(i,j)=1;
                    case 130
                        image(i,j,:) = 0.5*M_K + 0.5*BACKGROUND;
                         reference_anomaly_map(i,j)=1;
                    case 160
                        image(i,j,:) = 0.5*M_M + 0.5*BACKGROUND;
                         reference_anomaly_map(i,j)=1;
                end
            case  KERNEL_SIZE_ONE_BKG_75_MIX_LOCATION
            switch i
                case 40
                    image(i,j,:) = 0.25*M_A + 0.75*BACKGROUND;
                     reference_anomaly_map(i,j)=1;
                case 70
                    image(i,j,:) = 0.25*M_B + 0.75*BACKGROUND;
                     reference_anomaly_map(i,j)=1;
                case 100
                    image(i,j,:) = 0.25*M_C + 0.75*BACKGROUND;
                     reference_anomaly_map(i,j)=1;
                case 130
                    image(i,j,:) = 0.25*M_K + 0.75*BACKGROUND;
                     reference_anomaly_map(i,j)=1;
                case 160
                    image(i,j,:) = 0.25*M_M + 0.75*BACKGROUND;
                     reference_anomaly_map(i,j)=1;
            end
        end
    end
end
for(i=1:n_bands)
    image(:,:,i) = awgn(image(:,:,i),SNR);
end
figure;imagesc(image(:,:,160)); title('Band 160 of Hsueh-mimicked image, with gaussian noise'); axis image;

%% RX testing 
matrix = hyperConvert2d(image);
r_rx=hyperRxDetector(matrix);
r_rx = hyperConvert3d(r_rx.', h, w, 1);
figure; imagesc(reference_anomaly_map); title(['Expected anomaly map'] ); axis image; colorbar;
hyperSaveFigure(gcf, sprintf(['%s\\hsueh_expected_anomaly_map' '.png' ], resultsDir));%
figure; imagesc(r_rx); title(['RX AD '] ); axis image; colorbar;
hyperSaveFigure(gcf, sprintf(['%s\\RX AD' '.png' ], resultsDir));%
max_value_rx= max(r_rx);
max_value_rx = max(max_value_rx);
% set 90% of max_value as an anomaly
treshold_rx = max_value_rx*0.90;
anomaly_map_rx=zeros(h,w);
correctly_predicted_anomalies_rx=0;
for i=1:h
    for j=1:w
        if r_rx(i,j) >=treshold_rx
            anomaly_map_rx(i,j)=1;
            if reference_anomaly_map(i,j)==1
                correctly_predicted_anomalies_rx = correctly_predicted_anomalies_rx+1;
            end
        end
    end
end
figure; imagesc(anomaly_map_rx); title(['RX anomaly map'] ); axis image; colorbar;
hyperSaveFigure(gcf, sprintf(['%s\\RX anomaly_map' '.png' ], resultsDir));%
%check difference RX-anomaly_map and reference anomaly map
difference_from_reference_rx = zeros(h,w);
false_anomalies_hsueh_rx=nnz(anomaly_map_rx) - nnz(reference_anomaly_map);
if(false_anomalies_hsueh_rx<0)
    false_anomalies_hsueh_rx=0;
end

for i=1:h
    for j=1:w
        difference_from_reference_rx(i,j)= (reference_anomaly_map(i,j)- anomaly_map_rx(i,j));
    end
end
figure; imagesc(difference_from_reference_rx); title(['false or undetected anomalies RX'] ); axis image; colorbar;
hyperSaveFigure(gcf, sprintf(['%s\\Undetected anomalies RX' '.png' ], resultsDir));%

nnz_rx = nnz(difference_from_reference_rx);
percent_predicted_anomalies_hsueh_rx = (correctly_predicted_anomalies_rx)/nnz((reference_anomaly_map));


%% LRX without anomaly removal
difference_from_reference_lrx = zeros(h,w);
counter_i=1;
anomaly_map_alrx=zeros(1,h*w);
for K=5:5:30 % 35 bugged
 difference_from_reference_lrx = zeros(h,w);
 r_lrx=hyperLRxDetectorCorr(matrix,K);
 r_lrx = hyperConvert3d(r_lrx.', h, w, 1);
 figure; imagesc(r_lrx); title(['LRX K=' num2str(K)] ); axis image; colorbar;
 hyperSaveFigure(gcf, sprintf(['%s\\LRX K=' num2str(K) '.png' ], resultsDir));%
 max_value_lrx= max(r_lrx);
 max_value_lrx = max(max_value_lrx);
% set 75% of max_value as an anomaly
treshold_lrx = max_value_lrx*0.75;
anomaly_map_lrx=zeros(h,w);
correctly_predicted_anomalies_lrx =0;
for i=1:h
    for j=1:w
        if r_lrx(i,j) >=treshold_lrx
            anomaly_map_lrx(i,j)=1;
            if reference_anomaly_map(i,j)==1
                correctly_predicted_anomalies_lrx = correctly_predicted_anomalies_lrx+1;
            end
        end
    end
end

 for i=1:h
    for j=1:w
        difference_from_reference_lrx(i,j)= (reference_anomaly_map(i,j)- anomaly_map_lrx(i,j));
    end
 end
 figure; imagesc(difference_from_reference_lrx); title(['False or undetected anomalies LRX, K=' num2str(K)] ); axis image; colorbar;
 hyperSaveFigure(gcf, sprintf(['%s\\false anomalies LRX K=' num2str(K) '.png' ], resultsDir));%
false_anomalies_hsueh_lrx(counter_i)=nnz(anomaly_map_lrx) - nnz(reference_anomaly_map);
if(false_anomalies_hsueh_lrx(counter_i)<0)
    false_anomalies_hsueh_lrx(counter_i)=0;
end
percent_predicted_anomalies_hsueh_lrx(counter_i) = correctly_predicted_anomalies_lrx/nnz(reference_anomaly_map);
%percent_predicted_anomalies_hsueh_lrx(counter_i) = (nnz(anomaly_map_lrx) - false_anomalies_hsueh_lrx)/nnz(reference_anomaly_map);
 nnz_lrx(counter_i) = nnz(difference_from_reference_lrx);
 counter_i= counter_i +1;
end
%%
%LRX with anomaly removal
difference_from_reference_lrx_ad_remov = zeros(h,w);
counter_i=1;
correctly_predicted_anomalies_alrx=0;
treshold =250;
for K=5:5:35
 difference_from_reference_lrx_ad_remov = zeros(h,w);   
 [r_lrx_ad_remov,anomaly_map_alrx,location_of_anomalies,lsllsl]=hyperLRX_anomaly_set_remover(matrix,K,treshold);
 r_lrx_ad_remov = hyperConvert3d(r_lrx_ad_remov.', h, w, 1);
 figure; imagesc(r_lrx_ad_remov); title(['ALRX AD K=' num2str(K)] ); axis image; colorbar;
 hyperSaveFigure(gcf, sprintf(['%s\\ALRX AD K=' num2str(K) '.png' ], resultsDir));%
 anomaly_map_alrx_2d = hyperConvert3d(anomaly_map_alrx.',h,w,1);
     for i=1:h
        for j=1:w
            if anomaly_map_alrx_2d(i,j) ==1 && reference_anomaly_map(i,j) ==1
            correctly_predicted_anomalies_alrx = correctly_predicted_anomalies_alrx +1;
            end 
        end
     end
 false_anomalies_hsueh_alrx(counter_i)=nnz(anomaly_map_alrx) - correctly_predicted_anomalies_alrx;
% if(false_anomalies_hsueh_alrx(counter_i)<0)
%     false_anomalies_hsueh_alrx=0;
% end
%percent_predicted_anomalies_hsueh_alrx(counter_i) = (nnz(anomaly_map_alrx) - false_anomalies_hsueh_alrx)/nnz(reference_anomaly_map);    
percent_predicted_anomalies_hsueh_alrx(counter_i) = correctly_predicted_anomalies_alrx/nnz(reference_anomaly_map);    
 
figure; imagesc(difference_from_reference_rx); title(['False or undetected anomalies ALRX AD, K=' num2str(K)] ); axis image; colorbar;
 hyperSaveFigure(gcf, sprintf(['%s\\false anomalies ALRX AD K=' num2str(K) '.png' ], resultsDir));%
 nnz_lrx(counter_i) = nnz(difference_from_reference_lrx_ad_remov);
 counter_i= counter_i +1;
end


%% ACAD
matrix = hyperConvert2d(image);
treshold= 0.9;
[r_acad,anomaly_map,not_used]=hyperACAD(matrix,treshold);
r_acad = hyperConvert3d(r_acad.', h, w, 1);
figure; imagesc(reference_anomaly_map); title(['Expected anomaly map'] ); axis image; colorbar;
figure; imagesc(r_acad); title(['ACAD '] ); axis image; colorbar;
figure; imagesc(anomaly_map); title(['Anomaly map '] ); axis image; colorbar;

%% format data
r_acad_formatted = zeros(w,h);
for i =1: w
    for j=1:h
        if r_acad(i,j) <0
            r_acad_formatted(i,j)=0;
        elseif r_acad(i,j) >0 && ~isinf(r_acad(i,j))
            r_acad_formatted(i,j) = r_acad(i,j);
        else
            r_acad_formatted(i,j)=0;
        end
    end
end
