%for Cuprite scene
clc; close all; clear;
h=100;
w= 614;
load('groundTruth_Cuprite_nEnd12.mat','-mat');
M_endmembers=M;
goodBands = [10:100 116:150 180:216]; % for AVIRIS with 224 channels
M_endmembers=M(goodBands,:);
[n_bands,k] = size(M_endmembers);
image = zeros(h,w,n_bands);
reference_anomaly_map = zeros(h,w);
% Setting background
for i=1:h
    for j=1:w
        dice = randi(6);
        if dice>4
            image(i,j,:)= M_endmembers(:,1); %setting background to alunite
        elseif dice>2 
            image(i,j,:)= M_endmembers(:,6); %setting background to Kalonite
        else 
            image(i,j,:)= M_endmembers(:,10);% setting background to pyrope
        end
%          rN=rand;
%          image(i,j,:) = rN*M_endmembers(:,1) +0.2*M_endmembers(:,3)+0.2*M_endmembers(:,4)+0.2*M_endmembers(:,7)+(1-rN)*M_endmembers(:,12);
     end 
end

%imnoise(image,'gaussian',1);


%setting 50 random pixels to be an anomaly
% for i=1: 50
%     h_index=randi(h);
%     w_index= randi(w);
%     signature_index = randi([2 12]);
%     image(:,h_index,w_index) = M_endmembers(:,signature_index);
%     anomaly_map(h_index,w_index)=1;
% end

%create kernels with anomalies of size 1, 5, 10,15, 20, 25 in columns 5, 20,50,100, 400,
%600, in row 35 and 70 
%column locations
KERNEL_SIZE_ONE_LOCATION =50;
KERNEL_SIZE_TWO_LOCATION = 100;
KERNEL_SIZE_FIVE_LOCATION =150;
KERNEL_SIZE_TEN_LOCATION =250;
KERNEL_SIZE_FIFTEEN_LOCATION =350;
KERNEL_SIZE_TWENTY_LOCATION =450;
KERNEL_SIZE_TWENTYFIVE_LOCATION =550;
for i=1:h
    if(mod(i,35)==0)
        image(i,KERNEL_SIZE_ONE_LOCATION,:) = M_endmembers(:,3);
        reference_anomaly_map(i,KERNEL_SIZE_ONE_LOCATION)=1;
        
        image(i,KERNEL_SIZE_TWO_LOCATION,:) = M_endmembers(:,3);
        reference_anomaly_map(i,KERNEL_SIZE_TWO_LOCATION)=1;
        
        image(i,KERNEL_SIZE_TWO_LOCATION +1,:) = M_endmembers(:,3);
        reference_anomaly_map(i,KERNEL_SIZE_TWO_LOCATION +1)=1;
        
        image(i+1,KERNEL_SIZE_TWO_LOCATION,:) = M_endmembers(:,3);
        reference_anomaly_map(i+1,KERNEL_SIZE_TWO_LOCATION)=1;
        
        image(i+1,KERNEL_SIZE_TWO_LOCATION +1,:) = M_endmembers(:,3);
        reference_anomaly_map(i+1,KERNEL_SIZE_TWO_LOCATION +1)=1;
        for(j=5:5:25)
            for row =i-floor(j/2):i+floor(j/2)
                switch j
                      case 5
                          colcenter =KERNEL_SIZE_FIVE_LOCATION;
                      case 10
                          colcenter =KERNEL_SIZE_TEN_LOCATION;
                      case 15
                          colcenter =KERNEL_SIZE_FIFTEEN_LOCATION;
                      case 20
                          colcenter =KERNEL_SIZE_TWENTY_LOCATION;
                      case 25
                          colcenter =KERNEL_SIZE_TWENTYFIVE_LOCATION;
                end
                for col =colcenter-floor(j/2):colcenter+floor(j/2)
                    image(row,col,:) = M_endmembers(:,3);
                    reference_anomaly_map(row,col)=1;
                end
            end
        end
    end
end
imnoise(image,'gaussian',1);



 matrix=hyperConvert2d(image);
 
 
% r_rx = hyperRxDetector(matrix);
% r_rx_2d = hyperConvert3d(r_rx.', h, w, 1);
% figure; imagesc(reference_anomaly_map); title(['Expected anomaly map'] ); axis image; colorbar;
% figure; imagesc(r_rx_2d); title(['RX AD results'] ); axis image; colorbar;
% 
% max_value_rx= max(r_rx);
% % set 75% of max_value as an anomaly
% treshold_rx = max_value_rx*0.75;
% anomaly_map_rx=zeros(h,w);
% for i=1:h
%     for j=1:w
%         if r_rx_2d(i,j) >=treshold_rx
%             anomaly_map_rx(i,j)=1;
%         end
%     end
% end
% figure; imagesc(anomaly_map_rx); title(['RX anomaly map'] ); axis image; colorbar;
% %check difference RX-anomaly_map and reference anomaly map
% difference_from_reference_rx = zeros(h,w);
% false_anomalies_hsueh_rx=nnz(anomaly_map_rx) - nnz(reference_anomaly_map);
% if(false_anomalies_hsueh_rx<0)
%     false_anomalies_hsueh_rx=0;
% end
% 
% for i=1:h
%     for j=1:w
%         difference_from_reference_rx(i,j)= (reference_anomaly_map(i,j)- anomaly_map_rx(i,j));
%     end
% end
% figure; imagesc(difference_from_reference_rx); title(['false or undetected anomalies RX'] ); axis image; colorbar;
% %hyperSaveFigure(gcf, sprintf(['%s\\Undetected anomalies RX' '.png' ], resultsDir));%
% 
% nnz_rx = nnz(difference_from_reference_rx);
% percent_predicted_anomalies_hsueh_rx = (nnz(anomaly_map_rx) - false_anomalies_hsueh_rx)/nnz(reference_anomaly_map);


% K=25;
% r_rlx =hyperLRxDetectorCorr(matrix,K);
% r_rlx_2d = hyperConvert3d(r_rlx.', h, w, 1);
% 
% %anomaly_map_2d = hyperConvert3d(anomaly_map.', 30, 30, 1);
% figure;imagesc(r_rlx_2d);title(['LRX AD detector, K= ' num2str(K) ]); axis image; colorbar;
% 
% treshold=100;



%%ACAD 
figure;imagesc(reference_anomaly_map);axis image; colorbar;

false_anomalies = zeros(1,10);
true_anomalies = zeros(1,10);
correctly_predicted_anomalies = zeros(1,10);
counter_i=1;
for treshold =0.1:0.1:1
[d_acad, anomaly_map,threshold_check_values] = hyperACAD(matrix,treshold);
d_acad_2d = hyperConvert3d(d_acad.', h, w, 1);
anomaly_map_2d = hyperConvert3d(anomaly_map.', h, w, 1);
figure;imagesc(d_acad_2d); title(['ACAD result, treshold = ' num2str(treshold) ] ); axis image; colorbar;
figure;imagesc(anomaly_map_2d);title(['ACAD anomaly result, treshold = ' num2str(treshold) ] );  axis image; colorbar;
 
for i=1:h
     for j =1: w
         if(anomaly_map_2d(i,j)==1  && reference_anomaly_map(i,j)==1)
             true_anomalies(counter_i) =true_anomalies(counter_i)+1;
         elseif anomaly_map_2d(i,j)==1
             false_anomalies(counter_i) =false_anomalies(counter_i)+1;
         end
     end 
%      end
 correctly_predicted_anomalies(counter_i) = true_anomalies(counter_i)/nnz(reference_anomaly_map);
end
counter_i=counter_i+1;
end
