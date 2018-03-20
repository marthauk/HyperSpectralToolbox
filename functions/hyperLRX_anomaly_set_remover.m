function [result, anomalies_detected,location_of_anomalies,last_local_anomalies_set] = hyperLRX_anomaly_set_remover(M,K,treshold)
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
tresh_LRX=treshold;
location_of_anomalies= zeros(N/2,1);

local_anomalies_set=0;
last_local_anomalies_set =0;
ROWS=100;
t_an=1;


flag_local_anomaly_found =0;
    for j=1:N
        autocorr = hyperCorrK(M,K,p);
         %adaptive_autocorr_inv = inv(autocorr - anomalies_detected_transpose_sum);
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
        
        lower_limit_matrix = j - floor(K/2);
        higher_limit_matrix = j + floor(K/2);

        % Check if index is out of bounds 
        if( lower_limit_matrix < 1)
            % for edges of the matrix, gonna assume that we just throw out points
            % outside of the edge, and use half the KERNEL
            lower_limit_matrix = 1;   
        end
        if ( higher_limit_matrix > N)
          % M(band, neighbouring_pixels) * (M(band, Neighbouring pixels)
          higher_limit_matrix = N;
        end
        
        if(any(local_anomalies_set))
            %just to check that it works
            last_local_anomalies_set = local_anomalies_set;
        end
        %resetting local_anomalies_set before using it the next iteration
        local_anomalies_set = 0;
        flag_local_anomaly_found =0;
        for i=1:t_an
           if flag_local_anomaly_found == 0
               for k=1:K
                  if(flag_local_anomaly_found==0)
                      if(location_of_anomalies(i) >lower_limit_matrix+(k-1)*ROWS & location_of_anomalies(i)< higher_limit_matrix+(k-1)*ROWS)
                          local_anomalies_set = local_anomalies_set + anomalies_detected(:,i)*anomalies_detected(:,i).';
                          flag_local_anomaly_found =1;
                          break;
                      end
                  end
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
