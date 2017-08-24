% pRF_compareAverageEccSigma
%
% This code will load data files saved out form pRF_loopSigmaVsEcc and the
% _VTC version of it as well in order to compare average sigma and average
% eccentricity for each ROI across groups. Produces Fig 3B.
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

% Which file would you like to visualize? 
fileNameVTC = 'bi_EccVsSigma_lineData_anyHemi_Ethresh12_VTC.mat';

dataDir = fullfile(curdir,'voxel_data');
dataDirVTC = dataDir;

% Where should we save output?
saveDir = fullfile(outputdir,'pRF_figures'); if ~exist(saveDir), mkdir(saveDir); end

roiNames = {'r-pFus','l-pFus','r-OTS1','l-OTS1'};

load(fullfile(dataDirVTC,fileNameVTC));

ecc    = NaN(length(lineData),4);
sig    = NaN(length(lineData),4);

%% Now extract VTC roi information

% Some subjects didn't have functional localizer data, so there are empty
% cells in lineData that we can skip (since the data structure into which
% we are placing data is NaN, we don't have to do anything). 

for i = 1:numel(lineData)

    for m = 1:numel(lineData{1,i})
        
        % Check if subject has any ROIs
        if isempty(lineData{1,i})
            continue
        else
            roi_name = lineData{1,i}(1,m).roi;
        end
        
        if strcmp('rh_pFus_Faces',roi_name)
            variance = lineData{1,i}(1,m).variance;
            eccent   = lineData{1,i}(1,m).ecc;
            sigma    = lineData{1,i}(1,m).sigma;
            eccent(variance<=vThresh)=NaN; eccent(eccent>eThresh)=NaN;
            sigma(variance<=vThresh) =NaN; sigma(isnan(eccent)) = NaN; 
            sigma(sigma<0.21) = NaN; sigma(sigma>15)=NaN;
            ecc(i,1) = nanmean(eccent);
            sig(i,1) = nanmean(sigma);
        
        elseif strcmp('lh_pFus_Faces',roi_name)
            variance = lineData{1,i}(1,m).variance;
            eccent   = lineData{1,i}(1,m).ecc;
            sigma    = lineData{1,i}(1,m).sigma;
            eccent(variance<=vThresh)=NaN; eccent(eccent>eThresh)=NaN;
            sigma(variance<=vThresh) =NaN; sigma(isnan(eccent)) = NaN; 
            sigma(sigma<0.21) = NaN; sigma(sigma>15)=NaN;
            ecc(i,2) = nanmean(eccent);
            sig(i,2) = nanmean(sigma);

            
        elseif strcmp('rh_OTS1_WordsNumbers',roi_name)
            variance = lineData{1,i}(1,m).variance;
            eccent   = lineData{1,i}(1,m).ecc;
            sigma    = lineData{1,i}(1,m).sigma;
            eccent(variance<=vThresh)=NaN; eccent(eccent>eThresh)=NaN;
            sigma(variance<=vThresh) =NaN; sigma(isnan(eccent)) = NaN; 
            sigma(sigma<0.21) = NaN; sigma(sigma>15)=NaN;
            ecc(i,3) = nanmean(eccent);
            sig(i,3) = nanmean(sigma);
            
        elseif strcmp('lh_OTS1_WordsNumbers',roi_name)
            variance = lineData{1,i}(1,m).variance;
            eccent   = lineData{1,i}(1,m).ecc;
            sigma    = lineData{1,i}(1,m).sigma;
            eccent(variance<=vThresh)=NaN; eccent(eccent>eThresh)=NaN;
            sigma(variance<=vThresh) =NaN; sigma(isnan(eccent)) = NaN; 
            sigma(sigma<0.21) = NaN; sigma(sigma>15)=NaN;
            ecc(i,4) = nanmean(eccent);
            sig(i,4) = nanmean(sigma);
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

roiNames = {'r-pFus','l-pFus','r-OTS1','l-OTS1'};
grouping = [4 1 5 2 4 1 5 2];
dotPlots = [3.85 0.85 4.85 1.85 4.15 1.15 5.15 2.15];
      
colors = [1    0.6  0.6;
          0.8  0    0;
          0.7  0.7    1;
          0    0.2  0.8;
          1    0.6  0.6;
          0.8  0    0;
          0.7  0.7    1;
          0    0.2  0.8];

f =  figure('Position',[100 100 1000 300]); 
colorby = 'bar'; width = 0.3;
s1 = subplot_tight(1,2,2,[0.12 0.06]);
[ph, eh] = errorbargraph(meanSig,steSig,grouping,colors,width,'bar');
copyobj(ph,s1); copyobj(eh,s1); close(gcf); hold on
plot(dotPlots,allSig,'o','MarkerEdgeColor',[0.6 0.6 0.6])  
ylabel('pRF size (dva)','FontSize',16)
set(gca,'xtick',[1:1:5],'xticklabel',{'pFus', 'pOTS', '', 'pFus', 'pOTS'},'FontSize',16)
box off

s2 = subplot_tight(1,2,1,[0.12 0.06]);
[ph, eh] = errorbargraph(meanEcc,steEcc,grouping,colors,width,'bar');
copyobj(ph,s2); copyobj(eh,s2); close(gcf); hold on
plot(dotPlots,allEcc,'o','MarkerEdgeColor',[0.6 0.6 0.6])  
ylabel('pRF eccentricity (dva)','FontSize',16)

lastBar = length(unique(grouping));
set(gca,'xtick',[1:1:5],'xticklabel',{'pFus', 'pOTS', '', 'pFus', 'pOTS'},'FontSize',16)
box off

saveFile = ['pRF_mean_size_and_eccentricity_upTo' num2str(eThresh) 'deg_FaceWord_kidsVadults.fig'];
saveas(gcf,fullfile(saveDir,saveFile))


