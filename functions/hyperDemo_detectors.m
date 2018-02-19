function hyperDemo_detectors
% HYPERDEMO_DETECTORS Demonstrates target detector algorithms
clear; clc; dbstop if error; close all;
%--------------------------------------------------------------------------
%% Parameters
%resultsDir ='E:\One Drive\OneDrive for Business\NTNU\Master\Forked_MATLAB_hyperspectral_toolbox\MATLAB_Hyperspectral_toolbox\results';
%dataDir = 'E:\One Drive\OneDrive for Business\NTNU\Master\Forked_MATLAB_hyperspectral_toolbox\MATLAB_DEMO_hyperspectral\f970619t01p02r02c';

resultsDir ='E:\One Drive\OneDrive for Business\NTNU\Master\MATLAB_Hyperspectral_toolbox\results\ACAD_range300_1500';
dataDir = 'E:\One Drive\OneDrive for Business\NTNU\Master\MATLAB_DEMO_hyperspectral\f970619t01p02r02c';
%--------------------------------------------------------------------------

mkdir(resultsDir);

%% Read part of AVIRIS data file that we will further process
M = hyperReadAvirisRfl(sprintf('%s\\f970619t01p02_r02_sc02.a.rfl', dataDir), [1 100], [1 614], [1 224]);
M = hyperNormalize(M);
%M_test = imread('E:\One Drive\OneDrive for Business\NTNU\Master\MATLAB_DEMO_hyperspectral\220_band_aviris_june12\aviris_hyperspectral_data\19920612_AVIRIS_IndianPine_EW-line_R.tif');
%M_test = hyperNormalize(M_test);
%% Read AVIRIS .spc file
lambdasNm = hyperReadAvirisSpc(sprintf('%s\\f970619t01p02_r02.a.spc', dataDir));

%% Isomorph
[h, w, p] = size(M);
M = hyperConvert2d(M);

%% Resample AVIRIS image.
desiredLambdasNm = 400:(2400-400)/(224-1):2400;
M = hyperResample(M, lambdasNm, desiredLambdasNm);

%% Remove low SNR bands.
goodBands = [10:100 116:150 180:216];
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
%M_test_to = cast(M,'double');
%r = hyperLRxDetectorCorr(M,K,0);
%for tresh= 1500: 50 :2000
tresh = 200;
for tresh= 325: 25: 500
[r, anomalies_detected,treshold_check_values,location_of_anomalies] = hyperACAD(M,tresh);
N= 61400;
anomaly_map= zeros(1,N);

for i=1:1:N/2
    if (anomalies_detected(1,i)~= 0)
        pixel_pos_anomaly = location_of_anomalies(i);
        anomaly_map(pixel_pos_anomaly) = 100;  
    end
end

r = hyperConvert3d(r.', h, w, 1);
figure; imagesc(r); title(['ACAD Detector Results, tresh =' num2str(tresh) '.'] ); axis image;
    colorbar;
hyperSaveFigure(gcf, sprintf(['%s\\ACAD detector using abs' num2str(tresh) '.png' ], resultsDir));%

anomaly_map = hyperConvert3d(anomaly_map.', h, w, 1);    
figure; imagesc(anomaly_map); title([' Anomaly map, tresh =' num2str(tresh) '.'] ); axis image;
    colorbar;    
hyperSaveFigure(gcf, sprintf(['%s\\ACAD anomaly map' num2str(tresh) '.png' ], resultsDir));%
    
    

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

