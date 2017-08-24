% pRF_plotAveragedCoverage_bilateral.m
%
% This script will plot the bilateral averaged coverage of visual field
% maps V1 through VO1 as seen in Figure 2D
%
% JG 05/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
curdir = pwd; outputdir = fullfile(curdir,'output'); if ~exist(outputdir), mkdir(outputdir); end

matchFlag = true;
% Get subject indices
if matchFlag
    % These subjects are matched for variance explained by the pRF model in
    % V1. These subjects are the those present in the main figures.
    load(fullfile(curdir,'voxel_data','varMatched_indices.mat'));
else
    kidI = zeros(1,53); kidI(1:26)=1;
    adI  = zeros(1,53); adI(28:end)=1;
end


% Load right hemisphere data and flip
savePathData = fullfile(curdir,'voxel_data');
dataFile = fullfile(savePathData,'right_coverage_data_V1-VO2');

load(dataFile);

V1cov = coverage.V1; V2cov = coverage.V2;  V3cov = coverage.V3;
V4cov = coverage.V4; VO1cov= coverage.VO1; VO2cov= coverage.VO2;

rV1cov = flipdim(V1cov,2);
rV2cov = flipdim(V2cov,2);
rV3cov = flipdim(V3cov,2);
rV4cov = flipdim(V4cov,2);
rVO1cov = flipdim(VO1cov,2);
rVO2cov = flipdim(VO2cov,2);

kV1covr = nanmean(rV1cov(:,:,kidI),3);
aV1covr = nanmean(rV1cov(:,:,adI),3);

kV2covr = nanmean(rV2cov(:,:,kidI),3);
aV2covr = nanmean(rV2cov(:,:,adI),3);

kV3covr = nanmean(rV3cov(:,:,kidI),3);
aV3covr = nanmean(rV3cov(:,:,adI),3);

kV4covr = nanmean(rV4cov(:,:,kidI),3);
aV4covr = nanmean(rV4cov(:,:,adI),3);

kVO1covr = nanmean(rVO1cov(:,:,kidI),3);
aVO1covr = nanmean(rVO1cov(:,:,adI),3);

kVO2covr = nanmean(rVO2cov(:,:,kidI),3);
aVO2covr = nanmean(rVO2cov(:,:,adI),3);

% Now load left 
dataFile = fullfile(savePathData,'left_coverage_data_V1-VO2');
load(dataFile);

V1cov = coverage.V1; V2cov = coverage.V2;  V3cov = coverage.V3;
V4cov = coverage.V4; VO1cov= coverage.VO1; VO2cov= coverage.VO2;


kV1covl = nanmean(V1cov(:,:,kidI),3);
aV1covl = nanmean(V1cov(:,:,adI),3);

kV2covl = nanmean(V2cov(:,:,kidI),3);
aV2covl = nanmean(V2cov(:,:,adI),3);

kV3covl = nanmean(V3cov(:,:,kidI),3);
aV3covl = nanmean(V3cov(:,:,adI),3);

kV4covl = nanmean(V4cov(:,:,kidI),3);
aV4covl = nanmean(V4cov(:,:,adI),3);

kVO1covl = nanmean(VO1cov(:,:,kidI),3);
aVO1covl = nanmean(VO1cov(:,:,adI),3);

kVO2covl = nanmean(VO2cov(:,:,kidI),3);
aVO2covl = nanmean(VO2cov(:,:,adI),3);

% Now average across hemispheres
kV1cov = (kV1covl + kV1covr) ./2;
kV2cov = (kV2covl + kV2covr) ./2;
kV3cov = (kV3covl + kV3covr) ./2;
kV4cov = (kV4covl + kV4covr) ./2;
kVO1cov = (kVO1covl + kVO1covr) ./2;
kVO2cov = (kVO2covl + kVO2covr) ./2;

aV1cov = (aV1covl + aV1covr) ./2;
aV2cov = (aV2covl + aV2covr) ./2;
aV3cov = (aV3covl + aV3covr) ./2;
aV4cov = (aV4covl + aV4covr) ./2;
aVO1cov = (aVO1covl + aVO1covr) ./2;
aVO2cov = (aVO2covl + aVO2covr) ./2;

% Now let's plot
f = figure('Position',[100 100 1200 600],'color','w');

% First kids
subplot_tight(2,5,1,[0.005,0.005]  );
s_createCoveragePlot_averaged(kV1cov,'V1_Children'); colorbar off; 

subplot_tight(2,5,2,[0.005,0.005]  );
s_createCoveragePlot_averaged(kV2cov,'V2_Children'); colorbar off; 

subplot_tight(2,5,3,[0.005,0.005]  );
s_createCoveragePlot_averaged(kV3cov,'V3_Children'); colorbar off; 

subplot_tight(2,5,4,[0.005,0.005]  );
s_createCoveragePlot_averaged(kV4cov,'V4_Children'); colorbar off; 

subplot_tight(2,5,5,[0.005,0.005]  );
s_createCoveragePlot_averaged(kVO1cov,'VO1_Children'); colorbar off; 


% Now adults
subplot_tight(2,5,6,[0.005,0.005]  );
s_createCoveragePlot_averaged(aV1cov,'V1_Adults'); colorbar off; 

subplot_tight(2,5,7,[0.005,0.005]  );
s_createCoveragePlot_averaged(aV2cov,'V2_Adults'); colorbar off; 

subplot_tight(2,5,8,[0.005,0.005]  );
s_createCoveragePlot_averaged(aV3cov,'V3_Adults'); colorbar off; 

subplot_tight(2,5,9,[0.005,0.005]  );
s_createCoveragePlot_averaged(aV4cov,'V4_Adults'); colorbar off; 

subplot_tight(2,5,10,[0.005,0.005]  );
s_createCoveragePlot_averaged(aVO1cov,'VO1_Adults'); colorbar off; 


% Now save
savePathFig = fullfile(curdir,'output','pRF_figures'); if ~exist(savePathFig), mkdir(savePathFig); end
saveFigFile = fullfile(savePathFig,'bilateral_averageCoverage_allMaps_kidsVadults.fig');
saveas(gcf,saveFigFile)

