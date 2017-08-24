% pRF_loopSigmaVsEcc_VTC.m
%
% This script will plot the ecc vs. pRF size for a list of ROIs and loop
% through a list of subjects, saving out each subject's figures into a
% figure directory. It will also save out the regression information for
% each ROI that a subject (slope, intercept, age, variance, ecc, etc).
% This script calls s_rmPlotMultiEccSigma where all the action happens. This
% script can't be run for the purposes of figure reproduction, but was
% included as an example of how to loop mrVista code for other labs to use.
%
% JG 05/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
curdir = pwd; outputdir = fullfile(curdir,'output'); if ~exist(outputdir), mkdir(outputdir); end

% Set threshold variables here:
vethresh =  0.05;
eccthresh = [0 12];
sigthresh = [0.21 15]; 

% Which ROIs would you like to extract data and plot?
roiList = {'rh_pFus_Faces' 'lh_pFus_Faces' 'rh_mFus_Faces' 'lh_mFus_Faces' 'rh_OTS1_WordsNumbers' 'lh_OTS1_WordsNumbers' 'lh_OTS2_WordsNumbers' 'rh_PPA_Placeshouses' 'lh_PPA_Placeshouses'};

% Where do you want to save regression information?
saveDataDir = fullfile(curdir,'voxel_data');

% What should we call the saved data file?
saveFile = 'bi_EccVsSigma_lineData_anyHemi_Ethresh12_VTC.mat';

% This is where the figures of the slopes for each subject will be saved
saveDir = '/sni-storage/kalanit/biac2/kgs/projects/Longitudinal/FMRI/Retinotopy/results/pRF_figures/EccVsSigma/VTC';
saveDir = fullfile(outputdir,'images'); if ~exist(saveDir), mkdir(saveDir); end


% The structure in which we will store all subject data to save out:
lineData = {};

% Where do your subject's live?
dataDir = '/path/to/subjects';
cd(dataDir);

% We'll spit out problem subjects at the end if someone is missing
% something
errors = {};

% Let's loop through all subjects
sessions = {'sub1' 'sub2'};


for s = 1:length(sessions)

   cd(fullfile(dataDir,sessions{s}))
   load mrSESSION.mat;

    % Initialize a hidden gray view and the retinotopy model you've run
    view = initHiddenGray('Averages',1);
    rmPath = fullfile('Gray','Averages','retModel-cssFit-fFit.mat');
    view = rmSelect(view,1,rmPath);

   fprintf('\n\n\nProcessing %s\n\n\n',sessions{s})

   % The loop below will only include an ROI if it both exists and is
   % greater than 8 voxels (which is ~1 inplane functional voxel size).
   checkVec = zeros(size(roiList));
   for r = 1:length(roiList)
       if  exist(fullfile(dataDir,sessions{s},'3DAnatomy','ROIs',[roiList{r} '.mat']),'file')
           load(fullfile(dataDir,sessions{s},'3DAnatomy','ROIs',[roiList{r} '.mat']));
           if size(ROI.coords,2) >= 9
               checkVec(r) = 1; clear ROI;
           end
       end
   end
   
   roiListNew = roiList(find(checkVec));
   
   if ~isempty(roiListNew)
       
       view = loadROI(view,roiListNew);
      
       list = 1:length(viewGet(view, 'ROIs'));
       
       % Run regression and store data into a structure
       data = s_rmPlotMultiEccSigma(view, list, vethresh, eccthresh, sigthresh);
   
       myTitle = title(sessions{s},'fontsize',16); set(myTitle,'interpreter','none');
       print('-r300','-dpng',fullfile(saveDir,sessions{s}))
       
       
       lineData{s} = data;
   else
       % Put a place holder if the subject doesn't have the ROI
       lineData{s} = '';
   end
   
   
   mrvCleanWorkspace; close all;

end

save(fullfile(saveDataDir,saveFile),'lineData')

cd(dataDir)

