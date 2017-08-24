% pRF_loopSigmaVsEcc.m
%
% This script will plot the ecc vs. pRF size for a list of ROIs and loop
% through a list of subjects, saving out each subject's figures into a
% figure directory. It will also save out the regression information for
% each ROI that a subject (slope, intercept, age, variance, etc). This
% script can't be run for the purposes of figure reproduction, but was
% included as an example of how to loop mrVista code for other labs to use.
%
% JG 05/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
curdir = pwd; outputdir = fullfile(curdir,'output'); if ~exist(outputdir), mkdir(outputdir); end

% Set threshold variables here:
vethresh =  0.05;
eccthresh = [0.5 6.5]; % leave out very close to fovea due to presence of space ship, and we want to minimize edging effects so we'll stop at 6.5
sigthresh = [0.21 15];

% The structure in which we will store all subject data to save out:
lineData = {};

% Define the ROI list. We defined bilaterally and each hemisphere in the
% event a subject only has one hemisphere's ROI
roiList = {'longi_bi_V1' 'longi_bi_V2' 'longi_bi_V3' 'longi_bi_hV4' 'longi_bi_VO1' 'longi_bi_VO2'};
roiListR= {'longi_rh_V1' 'longi_rh_V2' 'longi_rh_V3' 'longi_rh_hV4' 'longi_rh_VO1' 'longi_rh_VO2'};
roiListL= {'longi_lh_V1' 'longi_lh_V2' 'longi_lh_V3' 'longi_lh_hV4' 'longi_lh_VO1' 'longi_lh_VO2'};
saveDir = fullfile(outputdir,'images'); if ~exist(saveDir), mkdir(saveDir); end

% Name the output data file
str = num2str(vethresh); str = str(end-1:end);
saveFile = ['bi_EccVsSigma_lineData_anyHemi_vThresh_' str];

% Where do your subject's live?
dataDir = '/path/to/subjects';

% Where do you want to save regression information?
saveDataDir = fullfile(curdir,'voxel_data');

% We'll spit out problem subjects at the end if someone is missing
% something
errors = {};

% Make your list of subjects to analyze
sessions = {'sub1' 'sub2'};


for s = 1:length(sessions)
    
    cd(fullfile(dataDir,sessions{s}))
    
    % Initialize a hidden gray view and the retinotopy model you've run
    view = initHiddenGray('Averages',1);
    rmPath = fullfile('Gray','Averages','retModel-cssFit-fFit.mat');
    view = rmSelect(view,1,rmPath);
    
    fprintf('\n\n\nProcessing %s\n\n\n',sessions{s})
    
    % Load a subject's bilateral map. If it doesn't exist, load the
    % hemisphere in which they have it
    for r = 1:length(roiList)
        if exist(fullfile(dataDir,sessions{s},'3DAnatomy','ROIs',[roiList{r} '.mat']),'file')
            roiListNew{r} = roiList{r};
        elseif exist(fullfile(dataDir,sessions{s},'3DAnatomy','ROIs',[roiListR{r} '.mat']),'file')
            roiListNew{r} = roiListR{r};
        elseif exist(fullfile(dataDir,sessions{s},'3DAnatomy','ROIs',[roiListL{r} '.mat']),'file')
            roiListNew{r} = roiListL{r};
        end
    end
    
    view = loadROI(view,roiListNew);
    
    list = 1:length(viewGet(view, 'ROIs'));
    data = s_rmPlotMultiEccSigma(view, list, vethresh, eccthresh, sigthresh);
    
    myTitle = title(sessions{s},'fontsize',16); set(myTitle,'interpreter','none');
    print('-r300','-dpng',fullfile(saveDir,sessions{s}))
    
    % Place subject's data into larger structure
    lineData{s} = data;
    
    mrvCleanWorkspace; close all;
    
end

save(fullfile(saveDataDir,saveFile),'lineData')

cd(dataDir)

