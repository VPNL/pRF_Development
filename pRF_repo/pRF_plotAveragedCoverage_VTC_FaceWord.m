% pRF_plotAveragedCoverage_VTC_FaceWord.m
%
% This script will reproduce Figure 4 panels A-D
%
% JG 05/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
curdir = pwd; outputdir = fullfile(curdir,'output'); if ~exist(outputdir), mkdir(outputdir); end

% Recalculate map coverage (true) or load previously save coverage (false):
calculate = false;

% Match subject variance when plotting?
matchFlag = true;
varThresh = 0.05; % Only take voxels with this variance explained or higher
sigThresh = 0.21; % Only use voxels with sigma size larger than this
eThresh = [0 12];

% where do we want to save figures and data?
savePathData = fullfile(curdir,'voxel_data');
savePathFig = fullfile(curdir,'output','pRF_figures'); if ~exist(savePathFig), mkdir(savePathFig); end


%% Calculate map coverage IF needed, otherwise next section loads previous data
if calculate
    
    maps = {'rh_pFus_Faces' 'lh_pFus_Faces' 'rh_mFus_Faces' 'lh_mFus_Faces' 'rh_OTS1_WordsNumbers' 'lh_OTS1_WordsNumbers' 'lh_OTS2_WordsNumbers'};
    
    
    for i = 1:length(sessions)
        
        cd(fullfile(sessionDir,sessions{i}))
        load mrSESSION.mat;
        
        % We need to load Averages or Averages_nomo for each sub
        nomoFlag = false;
        for d = 1:numel(dataTYPES)
            if strmatch('Averages_nomo', dataTYPES(1,d).name)
                nomoFlag = true;
            end
        end
        
        clear dataTYPES mrSESSION vANATOMYPATH;
        
        
        for m = 1:length(maps) % if sub has map, put pRF_cov into data, otherwise fill with NaNs
            
            % Now load the appropriate datatype
            if nomoFlag
                view = initHiddenGray('Averages_nomo',1);
                rmPath = fullfile('Gray','Averages_nomo','retModel-cssFitNoMo-fFit.mat');
                view = rmSelect(view,1,rmPath);
            else
                view = initHiddenGray('Averages',1);
                rmPath = fullfile('Gray','Averages','retModel-cssFit-fFit.mat');
                view = rmSelect(view,1,rmPath);
            end
            
            if exist(fullfile('3DAnatomy','ROIs',[maps{m} '.mat']))
                view = loadROI(view,maps{m});
                if dilate && i<=26
                    radius = 2; name = [maps{m} '_' num2str(radius) 'mmDilated']; scriptFlag = true;
                    [view roiDil] = roiDilate(view, view.ROIs(1), radius, name, 'k', scriptFlag)
                elseif disk
                    if m==1, radius=6;, elseif m==2, radius=5;, end;
                    name = [maps{m} '_' num2str(radius) 'mmDisk']; select=1; color='r';, grayCoordStart = 'roi'; addFlag=true;
                    [view,ROIdisk, layers] = makeROIdiskGray(view, radius, name, select, color, grayCoordStart,addFlag);
                end
                
                pRF_COV = s_rmPlotCoverage_VTC(view,'prf_size',1,'fieldRange',15,'eccthresh',eThresh,'method','maximum profile','nboot',50,'normalizeRange',1,'smoothSigma',1,'cothresh',varThresh,'weight','variance explained','sigmathresh',sigThresh,'addcenters',1,'newfig',-1);
                
                if strmatch('rh_pFus_Faces',maps{m})
                    coverage.rpfus(:,:,i) = pRF_COV;
                elseif strmatch('lh_pFus_Faces',maps{m})
                    coverage.lpfus(:,:,i) = pRF_COV;
                elseif strmatch('rh_mFus_Faces',maps{m})
                    coverage.rmfus(:,:,i) = pRF_COV;
                elseif strmatch('lh_mFus_Faces',maps{m})
                    coverage.lmfus(:,:,i) = pRF_COV;
                elseif strmatch('rh_OTS1_WordsNumbers',maps{m})
                    coverage.rots1(:,:,i) = pRF_COV;
                elseif strmatch('lh_OTS1_WordsNumbers',maps{m})
                    coverage.lots1(:,:,i) = pRF_COV;
                elseif strmatch('lh_OTS2_WordsNumbers',maps{m})
                    coverage.lots2(:,:,i) = pRF_COV;
                end
            else
                if strmatch('rh_pFus_Faces',maps{m})
                    coverage.rpfus(:,:,i) = NaN(128,128);
                elseif strmatch('lh_pFus_Faces',maps{m})
                    coverage.lpfus(:,:,i) = NaN(128,128);
                elseif strmatch('rh_mFus_Faces',maps{m})
                    coverage.rmfus(:,:,i) = NaN(128,128);
                elseif strmatch('lh_mFus_Faces',maps{m})
                    coverage.lmfus(:,:,i) = NaN(128,128);
                elseif strmatch('rh_OTS1_WordsNumbers',maps{m})
                    coverage.rots1(:,:,i) = NaN(128,128);
                elseif strmatch('lh_OTS1_WordsNumbers',maps{m})
                    coverage.lots1(:,:,i) = NaN(128,128);
                elseif strmatch('lh_OTS2_WordsNumbers',maps{m})
                    coverage.lots2(:,:,i) = NaN(128,128);
                end
            end
        end
        
        mrvCleanWorkspace; clear pRF_COV;
        
    end
    
    % Let's save the coverage info so we don't have to recalculate
    
    str1 = num2str(varThresh);
    str2 = num2str(eThresh(2));
    saveFile = fullfile(savePathData,['coverage_data_VTC_Vthresh' str1(end-1:end) '_Ethresh' str2]);
    
    save(saveFile,'coverage');
    
    %% Let's load previously saved coverage if already calculated
else
    
    str1 = num2str(varThresh);
    str2 = num2str(eThresh(2));
    
    dataFile = fullfile(savePathData,['coverage_data_VTC_Vthresh' str1(end-1:end) '_Ethresh' str2]);
    load(dataFile);
    
end

%% Now plot

% Get subject indices
if matchFlag
    % These subjects are matched for variance explained by the pRF model in
    % V1. These subjects are the those present in the main figures.
    load(fullfile(curdir,'voxel_data','varMatched_indices.mat'));
else
    kidI = zeros(1,53); kidI(1:26)=1;
    adI  = zeros(1,53); adI(28:end)=1;
end

% Now we need to extract the coverage maps for kids and adults and
% count the number of subjects going into each coverage average
rpfusCov = coverage.rpfus;
kNrpfus =  sum(~isnan(squeeze(rpfusCov(1,1,kidI))));
aNrpfus =  sum(~isnan(squeeze(rpfusCov(1,1,adI))));

lpfusCov = coverage.lpfus;
kNlpfus =  sum(~isnan(squeeze(lpfusCov(1,1,kidI))));
aNlpfus =  sum(~isnan(squeeze(lpfusCov(1,1,adI))));

rmfusCov = coverage.rmfus;
kNrmfus =  sum(~isnan(squeeze(rmfusCov(1,1,kidI))));
aNrmfus =  sum(~isnan(squeeze(rmfusCov(1,1,adI))));

lmfusCov = coverage.lmfus;
kNlmfus =  sum(~isnan(squeeze(lmfusCov(1,1,kidI))));
aNlmfus =  sum(~isnan(squeeze(lmfusCov(1,1,adI))));

rots1Cov = coverage.rots1;
kNrots1 =  sum(~isnan(squeeze(rots1Cov(1,1,kidI))));
aNrots1 =  sum(~isnan(squeeze(rots1Cov(1,1,adI))));

lots1Cov = coverage.lots1;
kNlots1 =  sum(~isnan(squeeze(lots1Cov(1,1,kidI))));
aNlots1 =  sum(~isnan(squeeze(lots1Cov(1,1,adI))));

lots2Cov = coverage.lots2;
kNlots2 =  sum(~isnan(squeeze(lots2Cov(1,1,kidI))));
aNlots2 =  sum(~isnan(squeeze(lots2Cov(1,1,adI))));


% Now let's produce the average maximum coverage plots
krpfuscov = nanmean(rpfusCov(:,:,kidI),3);
arpfuscov = nanmean(rpfusCov(:,:,adI),3);

klpfuscov = nanmean(lpfusCov(:,:,kidI),3);
alpfuscov = nanmean(lpfusCov(:,:,adI),3);

krmfuscov = nanmean(rmfusCov(:,:,kidI),3);
armfuscov = nanmean(rmfusCov(:,:,adI),3);

klmfuscov = nanmean(lmfusCov(:,:,kidI),3);
almfuscov = nanmean(lmfusCov(:,:,adI),3);

krots1cov = nanmean(rots1Cov(:,:,kidI),3);
arots1cov = nanmean(rots1Cov(:,:,adI),3);

klots1cov = nanmean(lots1Cov(:,:,kidI),3);
alots1cov = nanmean(lots1Cov(:,:,adI),3);

klots2cov = nanmean(lots2Cov(:,:,kidI),3);
alots2cov = nanmean(lots2Cov(:,:,adI),3);


%%%%%% Now we will plot the coverages  
f = figure('Position',[100 100 1400 700],'color','w');

subplot_tight(2,4,1,[0.03,0.03]);
s_createCoveragePlot_averaged(klots1cov,['lOTS1-Words_Children' num2str(kNlots1)]); colorbar off;

subplot_tight(2,4,2,[0.03,0.03]);
s_createCoveragePlot_averaged(alots1cov,['lOTS1-Words_Adults' num2str(aNlots1)]); colorbar off;

subplot_tight(2,4,5,[0.03,0.03]);
s_createCoveragePlot_averaged(klpfuscov,['lpFus-Faces_Children ' num2str(kNlpfus)]); colorbar off;

subplot_tight(2,4,6,[0.03,0.03]);
s_createCoveragePlot_averaged(alpfuscov,['lpFus-Faces_Adults ' num2str(aNlpfus)]); colorbar off;

subplot_tight(2,4,3,[0.03,0.03]);
s_createCoveragePlot_averaged(krots1cov,['rOTS1-Words_Children' num2str(kNrots1)]); colorbar off;

subplot_tight(2,4,4,[0.03,0.03]);
s_createCoveragePlot_averaged(arots1cov,['rOTS1-Words_Adults' num2str(aNrots1)]); colorbar off;

subplot_tight(2,4,7,[0.03,0.03]);
s_createCoveragePlot_averaged(krpfuscov,['rpFus-Faces_Children ' num2str(kNrpfus)]); colorbar off;

subplot_tight(2,4,8,[0.03,0.03]);
s_createCoveragePlot_averaged(arpfuscov,['rpFus-Faces_Adults ' num2str(aNrpfus)]); colorbar off;


% Now save the figures
str = num2str(varThresh);

saveFigFile = fullfile(savePathFig,['averageCoverage_VTC_FaceWord_kidsVadults_thresh' str(end-1:end) '_notMaxed.fig']);
saveFigFile2 = fullfile(savePathFig,['averageCoverage_VTC_FaceWord_kidsVadults_thresh' str(end-1:end) '_notMaxed.png']);

saveas(gcf,saveFigFile)
saveas(gcf,saveFigFile2)
