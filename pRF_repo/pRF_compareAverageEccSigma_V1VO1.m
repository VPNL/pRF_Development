% pRF_compareAverageEccSigma
%
% This code will load data files saved out form pRF_loopSigmaVsEcc and the
% _VTC version of it as well in order to compare average sigma and average
% eccentricity for each ROI across groups. Will produce Fig. 2AB.
%
% JG 06/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
curdir = pwd; outputdir = fullfile(curdir,'output'); if ~exist(outputdir), mkdir(outputdir); end

% Match for mean variance explained in V1?
matchFlag = true;

% How should we threshold?
eThresh = 7; % Let's not look at centers beyond our stimulus coverage
vThresh = 0.05; % Only look above this variance explained

% Which file would you like to visualize (made from pRF_loopSigmaVsEcc)? 
fileName = 'bi_EccVsSigma_lineData_anyHemi_vThresh_05.mat';

% Where is the data file you want to visualize stored?
dataDir = fullfile(curdir,'voxel_data');

% Where should we save output?
saveDir = fullfile(outputdir,'pRF_figures'); if ~exist(saveDir), mkdir(saveDir); end

% Let's start with retinotopic maps before doing VTC rois
load(fullfile(dataDir,fileName));

roiNames = {'V1','V2','V3','hV4','VO1'};
% Set up data structure, rows per subject, columns per ROI 
% Column 1 = v1, 2 = v2, 3= v3, 4 = v4
% 5 = vo1, 6= vo2, 7 =rpfus, 8 lpfus, 
% 9 = rots1, 10 = lots1, 
ecc    = NaN(length(lineData),5);
sig    = NaN(length(lineData),5);
varexp = NaN(1,length(lineData)); % we use this to variance match subjects


%% First let's extract V1 through VO1 information
for i = 1:numel(lineData)
    
    for m = 1:numel(lineData{1,i})
        
        if strcmp('V1',lineData{1,i}(1,m).roi(end-1:end))
            variance = lineData{1,i}(1,m).variance;
            eccent   = lineData{1,i}(1,m).ecc;
            sigma    = lineData{1,i}(1,m).sigma;
            eccent(variance<=vThresh)=NaN; eccent(eccent>eThresh)=NaN;
            sigma(variance<=vThresh) =NaN; sigma(isnan(eccent)) = NaN; 
            sigma(sigma<0.21) = NaN; sigma(sigma>15)=NaN;
            ecc(i,1) = nanmean(eccent);
            sig(i,1) = nanmean(sigma);
        
        elseif strcmp('V2',lineData{1,i}(1,m).roi(end-1:end))
            variance = lineData{1,i}(1,m).variance;
            eccent   = lineData{1,i}(1,m).ecc;
            sigma    = lineData{1,i}(1,m).sigma;
            eccent(variance<=vThresh)=NaN; eccent(eccent>eThresh)=NaN;
            sigma(variance<=vThresh) =NaN; sigma(isnan(eccent)) = NaN; 
            sigma(sigma<0.21) = NaN; sigma(sigma>15)=NaN;
            ecc(i,2) = nanmean(eccent);
            sig(i,2) = nanmean(sigma);
  
            variance = variance(variance>=0.1);
            varexp(1,i) = nanmean(variance); % We will store variance for V2 for variance matching
            
        elseif strcmp('V3',lineData{1,i}(1,m).roi(end-1:end))
            variance = lineData{1,i}(1,m).variance;
            eccent   = lineData{1,i}(1,m).ecc;
            sigma    = lineData{1,i}(1,m).sigma;
            eccent(variance<=vThresh)=NaN; eccent(eccent>eThresh)=NaN;
            sigma(variance<=vThresh) =NaN; sigma(isnan(eccent)) = NaN; 
            sigma(sigma<0.21) = NaN; sigma(sigma>15)=NaN;
            ecc(i,3) = nanmean(eccent);
            sig(i,3) = nanmean(sigma);
            
        elseif strcmp('V4',lineData{1,i}(1,m).roi(end-1:end))
            variance = lineData{1,i}(1,m).variance;
            eccent   = lineData{1,i}(1,m).ecc;
            sigma    = lineData{1,i}(1,m).sigma;
            eccent(variance<=vThresh)=NaN; eccent(eccent>eThresh)=NaN;
            sigma(variance<=vThresh) =NaN; sigma(isnan(eccent)) = NaN; 
            sigma(sigma<0.21) = NaN; sigma(sigma>15)=NaN;
            ecc(i,4) = nanmean(eccent);
            sig(i,4) = nanmean(sigma);
            
        elseif strcmp('VO1',lineData{1,i}(1,m).roi(end-2:end))
            variance = lineData{1,i}(1,m).variance;
            eccent   = lineData{1,i}(1,m).ecc;
            sigma    = lineData{1,i}(1,m).sigma;
            eccent(variance<=vThresh)=NaN; eccent(eccent>eThresh)=NaN;
            sigma(variance<=vThresh) =NaN; sigma(isnan(eccent)) = NaN; 
            sigma(sigma<0.21) = NaN; sigma(sigma>15)=NaN;
            ecc(i,5) = nanmean(eccent);
            sig(i,5) = nanmean(sigma);
            
        end
    end
end


%% Now let's plot

% Let us first calculate the subject indices
if matchFlag
    % These subjects are matched for variance explained by the pRF model in
    % V1. These subjects are the those present in the main figures.
    load(fullfile(curdir,'voxel_data','varMatched_indices.mat'));
else
    kidI = zeros(1,53); kidI(1:26)=1;
    adI  = zeros(1,53); adI(28:end)=1;
end

% Let's make a 2 panel plot with pRF size on the top, eccentricity on the
% bottom. 
kEcc = (ecc(kidI,:)); 
aEcc = (ecc(adI,:));  
kEcc(end:end+(length(aEcc)-length(kEcc)),:) = NaN;
allEcc = [kEcc aEcc];
meanEcc = nanmean(allEcc,1);
steEcc  = nanste(allEcc,1);

kSig = (sig(kidI,:)); 
aSig = (sig(adI,:));  
kSig(end:end+(length(aSig)-length(kSig)),:) = NaN;
allSig = [kSig aSig];
meanSig = nanmean(allSig,1);
steSig  = nanste(allSig,1);

grouping = [1 2 3 4 5 1 2 3 4 5];
dotPlots = [0.85 1.85 2.85 3.85 4.85 1.15 2.15 3.15 4.15 5.15];
colors = [1    0.6  0.6;
          0.8  0    0;
          1    0.8  0.6
          0.8  0.4  0
          1    1    0.4
          0.6  0.6  0
          0.7  1    0.4
          0.4  0.8  0
          0.7  1    1
          0    0.8  0.8];

f = figure('Position',[100 100 700 900]); 
   
colorby = 'bar'; width = 0.3;
s1 = subplot_tight(2,1,1,[0.06 0.08]);
[ph, eh] = errorbargraph(meanSig,steSig,grouping,colors,width,'bar');
copyobj(ph,s1); copyobj(eh,s1); close(gcf); hold on
plot(dotPlots,allSig,'o','MarkerEdgeColor',[0.6 0.6 0.6])  
ylabel('Mean pRF size (^o)','FontSize',20)
set(gca,'xtick',[1:1:5],'xticklabel',{''},'FontSize',20);
box off

s2 = subplot_tight(2,1,2,[0.06 0.08]);
[ph, eh] = errorbargraph(meanEcc,steEcc,grouping,colors,width,'bar');
copyobj(ph,s2); copyobj(eh,s2); close(gcf); hold on
plot(dotPlots,allEcc,'o','MarkerEdgeColor',[0.6 0.6 0.6])  

lastBar = length(unique(grouping));
set(gca,'xtick',[1:1:5],'xticklabel',{'V1', 'V2', 'V3', 'hV4', 'VO1'},'FontSize',20)
ylabel('Mean pRF ecc. (^o)','FontSize',20)
box off


saveFile = ['pRF_mean_size_and_eccentricity_upTo' num2str(eThresh) 'deg_V1-VO1_kidsVadults.fig'];
saveas(gcf,fullfile(saveDir,saveFile))

