function hyperDemo_detectors
% HYPERDEMO_DETECTORS Demonstrates target detector algorithms
clear; clc; dbstop if error; close all;
%--------------------------------------------------------------------------
%% Parameters
%resultsDir ='E:\One Drive\OneDrive for Business\NTNU\Master\Forked_MATLAB_hyperspectral_toolbox\MATLAB_Hyperspectral_toolbox\results';
%dataDir = 'E:\One Drive\OneDrive for Business\NTNU\Master\Forked_MATLAB_hyperspectral_toolbox\MATLAB_DEMO_hyperspectral\f970619t01p02r02c';

resultsDir =['E:\One Drive\OneDrive for Business\NTNU\Master\Anomaly detection results\MATLAB\LRX\real image data Cuprite scene\' ,datestr(now, 'dd-mmm-yyyy')];
dataDir = 'E:\One Drive\OneDrive for Business\NTNU\Master\MATLAB_DEMO_hyperspectral\f970619t01p02r02c';
%--------------------------------------------------------------------------

mkdir(resultsDir);

%% Read part of AVIRIS data file that we will further process
M = hyperReadAvirisRfl(sprintf('%s\\f970619t01p02_r02_sc02.a.rfl', dataDir), [1 100], [1 614], [1 224]);
%M = hyperReadAvirisRfl(sprintf('%s\\f970619t01p02_r02_sc04.a.rfl', dataDir), [1 100], [1 614], [1 224]);

M = hyperNormalize(M);
%% Read AVIRIS .spc file
lambdasNm = hyperReadAvirisSpc(sprintf('%s\\f970619t01p02_r02.a.spc', dataDir));

%% Isomorph
[h, w, p] = size(M);
M = hyperConvert2d(M);
%KSC_2d = hyperConvert2d(KSC);
%M=KSC_2d;
%% Resample AVIRIS image.
desiredLambdasNm = 400:(2400-400)/(224-1):2400;
M = hyperResample(M, lambdasNm, desiredLambdasNm);

%% Remove low SNR bands.
goodBands = [10:100 116:150 180:216]; % for AVIRIS with 224 channels
%goodbands_KSC =[10:100 116:150];

%KSC_2d = KSC_2d(goodbands_KSC,:);
%p = length(goodbands_KSC);
M = M(goodBands, :);
p = length(goodBands);

%% Demonstrate difference spectral similarity measurements
M = hyperConvert3d(M, h, w, p);
target = squeeze(M(11, 77, :));
figure; plot(desiredLambdasNm(goodBands), target); grid on;
    title('Target Signature; Pixel (32, 257)'); 

M = hyperConvert2d(M);
  
%% RX Anomly Detector
%r = hyperRxDetector(M);
%r = hyperRxDetectorCor(M);
K=25;
resultsDir =['E:\One Drive\OneDrive for Business\NTNU\Master\Anomaly detection results\MATLAB\LRX\real image data Cuprite scene\' ,datestr(now, 'dd-mmm-yyyy')];
%r = hyperLRxDetectorCorr(M,K);
%g = ground_truth(h,614, M, M_endmembers);
%figure; imagesc(g);colorbar;
treshold = 500;
%for treshold= 500: 250 :2000
%[r,anomalies_detected, location_of_anomalies,last_local_anomalies_set]=hyperLRX_anomaly_set_remover(M,K,treshold);
%[r,anomalies_detected, location_of_anomalies,last_local_anomalies_set]=hyperLRX_anomaly_set_remover(KSC_2d,K,treshold);

 for treshold= 50: 25: 500
 [r, anomaly_map,location_of_anomalies] = hyperACAD(M,treshold);


r = hyperConvert3d(r.', h, w, 1);
figure; imagesc(r); title(['ACAD Detector Results, tresh =' num2str(treshold) '.'] ); axis image;
    colorbar;
 hyperSaveFigure(gcf, sprintf(['%s\\ACAD Detector Results, tresh' num2str(treshold) '.png' ], resultsDir));%
    
%figure; imagesc(r); title(['LRX removing anomalies,tresh =2000,K=25 .'] ); axis image;
%     colorbar;
%  figure; imagesc(r); title(['LRX Cuprite image data sc02 K=' num2str(K) '.'] ); axis image;
%     colorbar;
% hyperSaveFigure(gcf, sprintf(['%s\\LRX_K=25_treshold_500_KSC' num2str(treshold) '.png' ], resultsDir));%
% 
  anomaly_map = hyperConvert3d(anomaly_map.', h, w, 1);   
figure; imagesc(anomaly_map); title(['Anomaly map ACAD, treshold =' num2str(treshold) ', K=' num2str(K) '.'] ); axis image;
    colorbar;
 hyperSaveFigure(gcf, sprintf(['%s\\Anomaly Map ACAD  ' num2str(treshold) '.png' ], resultsDir));%



end

%% Constrained Energy Minimization (CEM)
r = hyperCem(M, target);
r = hyperConvert3d(r, h, w, 1);
figure; imagesc(abs(r)); title('CEM Detector Results'); axis image;
    colorbar;    
hyperSaveFigure(gcf, sprintf('%s\\cem detector.png', resultsDir));

%% Adaptive Cosine Estimator (ACE)
r = hyperAce(M, target);
r = hyperConvert3d(r, h, w, 1);
figure; imagesc(r); title('ACE Detector Results'); axis image;
    colorbar;
hyperSaveFigure(gcf, sprintf('%s\\ace detector.png', resultsDir));    

%% Signed Adaptive Cosine Estimator (S-ACE)
r = hyperSignedAce(M, target);
r = hyperConvert3d(r, h, w, 1);
figure; imagesc(r); title('Signed ACE Detector Results'); axis image;
    colorbar;
hyperSaveFigure(gcf, sprintf('%s\\signed ace detector.png', resultsDir));  

%% Matched Filter
r = hyperMatchedFilter(M, target);
r = hyperConvert3d(r, h, w, 1);
figure; imagesc(r); title('MF Detector Results'); axis image;
    colorbar;
hyperSaveFigure(gcf, sprintf('%s\\mf detector.png', resultsDir)); 

%% Generalized Likehood Ratio Test (GLRT) detector
r = hyperGlrt(M, target);
r = hyperConvert3d(r, h, w, 1);
figure; imagesc(r); title('GLRT Detector Results'); axis image;
    colorbar;
hyperSaveFigure(gcf, sprintf('%s\\cem detector.png', resultsDir));


%% Estimate background endmembers
U = hyperAtgp(M, 5);

%% Hybrid Unstructured Detector (HUD)
r = hyperHud(M, U, target);
r = hyperConvert3d(r, h, w, 1);
figure; imagesc(abs(r)); title('HUD Detector Results'); axis image;
    colorbar;
hyperSaveFigure(gcf, sprintf('%s\\hud detector.png', resultsDir));
   
%% Adaptive Matched Subspace Detector (AMSD)
r = hyperAmsd(M, U, target);
r = hyperConvert3d(r, h, w, 1);
figure; imagesc(abs(r)); title('AMSD Detector Results'); axis image;
    colorbar;
hyperSaveFigure(gcf, sprintf('%s\\amsd detector.png', resultsDir));    
figure; mesh(r); title('AMSD Detector Results');

%% Orthogonal Subspace Projection (OSP)
r = hyperOsp(M, U, target);
r = hyperConvert3d(r, h, w, 1);
figure; imagesc(abs(r)); title('OSP Detector Results'); axis image;
    colorbar;   
hyperSaveFigure(gcf, sprintf('%s\\osp detector.png', resultsDir));    

