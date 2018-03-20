generate_synthetic_image;
%Convert to 2-D coordinates;
matrix= hyperConvert2d(image);
%RX testing 
r_rx=hyperRxDetector(matrix);
r_rx = hyperConvert3d(r_rx.', h, w, 1);
% date_now = datetime('now');
% DateString= datestr(date_now);
% str= ['E:\One Drive\OneDrive for Business\NTNU\Master\Anomaly detection results\MATLAB\synthetic_images\' , DateString];
% resultsDir =join(str);
resultsDir= ['E:\One Drive\OneDrive for Business\NTNU\Master\Anomaly detection results\MATLAB\synthetic_images\' ,datestr(now, 'dd-mmm-yyyy')];
%resultsDir =regexprep(resultsDir,':d*','')
%resultsDir ='E:\One Drive\OneDrive for Business\NTNU\Master\Anomaly detection results\MATLAB\synthetic_images\lol' ;
[status, msg, msgID] = mkdir(resultsDir);

figure; imagesc(reference_anomaly_map); title(['Expected anomaly map'] ); axis image; colorbar;
hyperSaveFigure(gcf, sprintf(['%s\\expected_anomaly_map' '.png' ], resultsDir));%
figure; imagesc(r_rx); title(['RX'] ); axis image; colorbar;
hyperSaveFigure(gcf, sprintf(['%s\\RX AD' '.png' ], resultsDir));%
max_value_rx= max(r_rx);
% set 90% of max_value as an anomaly
treshold_rx = max_value_rx*0.9;
anomaly_map_rx=zeros(h,w);
for i=1:h
    for j=1:w
        if r_rx(i,j) >=treshold_rx
            anomaly_map_rx(i,j)=1;
        end
    end
end
figure; imagesc(anomaly_map_rx); title(['RX anomaly map'] ); axis image; colorbar;
hyperSaveFigure(gcf, sprintf(['%s\\RX anomaly_map' '.png' ], resultsDir));%
%check difference RX-anomaly_map and reference anomaly map
difference_from_reference_rx = zeros(h,w);
for i=1:h
    for j=1:w
        difference_from_reference_rx(i,j)= (reference_anomaly_map(i,j)- anomaly_map_rx(i,j));
    end
end
figure; imagesc(difference_from_reference_rx); title(['false or undetected anomalies RX'] ); axis image; colorbar;
hyperSaveFigure(gcf, sprintf(['%s\\Undetected anomalies RX' '.png' ], resultsDir));%
nnz_rx = nnz(difference_from_reference_rx);

%LRX without anomaly removal
difference_from_reference_lrx = zeros(h,w);
counter_i=1;
for K=5:5:35
 r_lrx=hyperLRxDetectorCorr(matrix,K);
 r_lrx = hyperConvert3d(r_lrx.', h, w, 1);
 figure; imagesc(r_lrx); title(['LRX K=' num2str(K)] ); axis image; colorbar;
 hyperSaveFigure(gcf, sprintf(['%s\\LRX K=' num2str(K) '.png' ], resultsDir));%
 max_value_lrx= max(r_lrx);
% set 90% of max_value as an anomaly
treshold_lrx = max_value_lrx*0.9;
anomaly_map_lrx=zeros(h,w);
for i=1:h
    for j=1:w
        if r_lrx(i,j) >=treshold_lrx
            anomaly_map_lrx(i,j)=1;
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
 nnz_lrx(counter_i) = nnz(difference_from_reference_lrx);
 counter_i= counter_i +1;
end

%LRX with anomaly removal
difference_from_reference_lrx_ad_remov = zeros(h,w);
counter_i=1;
for K=5:5:35
 [r_lrx_ad_remov,anomalies_detected,location_of_anomalies,lsllsl]=hyperLRX_anomaly_set_remover(matrix,K,1000);
 r_lrx_ad_remov = hyperConvert3d(r_lrx_ad_remov.', h, w, 1);
 figure; imagesc(r_lrx_ad_remov); title(['LRX AD removal K=' num2str(K)] ); axis image; colorbar;
 hyperSaveFigure(gcf, sprintf(['%s\\LRX AD removalK=' num2str(K) '.png' ], resultsDir));%
     for i=1:h
        for j=1:w
            difference_from_reference_lrx_ad_remov(i,j)= (reference_anomaly_map(i,j)- anomaly_map_rx(i,j));
        end
     end
 figure; imagesc(difference_from_reference_rx); title(['False or undetected anomalies LRX AD removal, K=' num2str(K)] ); axis image; colorbar;
 hyperSaveFigure(gcf, sprintf(['%s\\false anomalies LRX AD removal K=' num2str(K) '.png' ], resultsDir));%
 nnz_lrx(counter_i) = nnz(difference_from_reference_lrx_ad_remov);
 counter_i= counter_i +1;
end



