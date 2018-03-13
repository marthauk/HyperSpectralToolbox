function [result, anomalies_detected,location_of_anomalies] = hyperLRX_anomaly_set_remover(M,K)
% LRX anomaly detector, that also removes the detected anomalous targets
% causaly.
%   hyperLRxDetector performs the Local RX anomaly detector using Correlation
%   instead of covariance
%
% Usage
%   [result] = hyperRxDetector(M)
% Inputs
%   M  - 2D data matrix (p x N)
%   K  - Size of the kernel window, K x K
% Outputs
%   result - Detector output (1 x N)

[p, N] = size(M);
%waitbar for progress monitoring
h = waitbar(0,'Initializing waitbar ..');

% Compute correlation matrix of size K
% correlation matrix will be of size p x p
result = zeros(N, 1);
anomalies_detected=zeros(p,N/2);
anomalies_detected_transpose_sum = zeros(p,p);
%tresh_LRX = 6.0000e+14;
tresh_LRX=4000;
location_of_anomalies= zeros(N/2,1);

local_anomalies_set=0;
ROWS=100;
t_an=1;

    for j=1:N
        autocorr = hyperCorrK(M,K,p);
         adaptive_autocorr_inv = inv(autocorr - anomalies_detected_transpose_sum);
         adaptive_autocorr_inv = inv(autocorr - local_anomalies_set);
%result(j) = M(:,i).' * autocorrInv;
        
        result(j) = M(:,j).' * adaptive_autocorr_inv * M(:,j);
        if result(j) > tresh_LRX  
            % This pixel is an anomaly! Add it to the set of anomalies
            anomalies_detected(:,t_an) = M(:,j);
            location_of_anomalies(t_an)=j;
            anomalies_detected_transpose_sum = M(:,j)* M(:,j).' + anomalies_detected_transpose_sum;
            t_an = t_an + 1;
        end
        %if anomalies_detected_transpose_sum contains elements from outside
        %the KERNEL
        for( i=1:t_an)
           for k=1:K
              if(location_of_anomalies(t_an) >j-floor(K/2)+(k-1)*ROWS & location_of_anomalies(t_an)< j+floor(K/2)+(k-1)*ROWS)
                  local_anomalies_set = local_anomalies_set + anomalies_detected(:,t_an)*anomalies_detected(:,t_an).';
                  break;
              end
           end
            %if location_of_anomalies(t_an)<j-floor(K/2) | location_of_anomalies(t_an)>j+floor(K/2)
            %    anomalies_detected_transpose_sum =anomalies_detected_transpose_sum - M(:,t_an)*M(:,t_an).';
           %end 
        end 
        %result(j)= result(j) *  M(:,i);
        waitbar(j/N,h,'Updated progress');
    end


result = abs(result);

return;