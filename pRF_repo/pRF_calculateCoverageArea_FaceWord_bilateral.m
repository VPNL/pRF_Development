% pRF_calculateCoverageArea_FaceWord_bilateral.m
%
% This script was used to calculate the coverage area of the bilateral
% face and word regions
%
% JG 01/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
curdir = pwd; outputdir = fullfile(curdir,'output'); if ~exist(outputdir), mkdir(outputdir); end

% Recalculate map coverage (true) or load previously save coverage (false):
calculate = false;

% What is the AREA of your visual display in degrees of visual angle?
area = 49*pi;

% What threshold contour do you want to use?
thresh = 0.01;

% Match subject variance when plotting?
matchFlag = true;
varThresh = 0.05; % Only take voxels with this variance explained or higher
sigThresh = 0.21;
eThresh = [0 12];

% where do we want to save figures and data?
savePathData = fullfile(curdir,'voxel_data');
savePathFig = fullfile(curdir,'output','pRF_figures'); if ~exist(savePathFig), mkdir(savePathFig); end

%% Calculate map coverage IF needed, otherwise next section loads previous data
if calculate
    
facemaps = {'rh_pFus_Faces' 'lh_pFus_Faces'};
wordmaps = {'rh_OTS1_WordsNumbers' 'lh_OTS1_WordsNumbers'};

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
     
%%%%FACE COVERAGE CALCULATION%%%%
    checkVec = zeros(1,length(facemaps));
    for r = 1:length(facemaps)
        checkVec(r) = exist(fullfile('3DAnatomy','ROIs',[facemaps{r} '.mat']),'file');
    end
    
    % If the subject has none of the ROIs, set vals to NaN and move on to
    % word rois
    if sum(checkVec) == 0
        coverage.pfus(:,:,i) = NaN(128,128);
        coverage.pfus_size(i)= NaN;
    else
        % Now we will only consider those that exist for that subject
        roiListNew = facemaps(find(checkVec));
        
        % Now load these ROIs into the hiddenGray
        for r = 1:length(roiListNew)
            view = loadROI(view,roiListNew{r});
        end
        
        rois = {}; % rois need to be stored in a cell
        for r=1:length(roiListNew)
            rois{r} = view.ROIs(r);
        end
        [view, roi, ~] = combineROIs(view, rois, 'union', 'pFus_faces_bilateral');
        
        % now calculate the coverage 
        pRF_COV = s_rmPlotCoverage_VTC(view,'prf_size',1,'fieldRange',15,'eccthresh',eThresh,'method','density','nboot',50,'normalizeRange',0,'smoothSigma',1,'cothresh',varThresh,'weight','variance explained','sigmathresh',sigThresh,'addcenters',1,'newfig',-1);
        coverage.pfus(:,:,i) = pRF_COV;
        coverage.pfus_size(i) = length(view.ROIs(end).coords);
    end
    
    clear pRF_COV;
        
    
    
    
    
%%%%WORD COVERAGE CALCULATION%%%%

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

    checkVec = zeros(1,length(wordmaps));
    for r = 1:length(wordmaps)
        checkVec(r) = exist(fullfile('3DAnatomy','ROIs',[wordmaps{r} '.mat']),'file');
    end

    % If the subject has none of the ROIs, set vals to NaN and move to next sub
    if sum(checkVec) == 0
        coverage.ots1(:,:,i) = NaN(128,128);
        coverage.ots1(i) = NaN;
    else
        % Now we will only consider those that exist for that subject
        roiListNew = wordmaps(find(checkVec));
        
        % Now load these ROIs into the hiddenGray
        for r = 1:length(roiListNew)
            view = loadROI(view,roiListNew{r});
        end
        
        rois = {}; % rois need to be stored in a cell
        for r=1:length(roiListNew)
            rois{r} = view.ROIs(r);
        end
        [view, roi, ~] = combineROIs(view, rois, 'union', 'OTS1_chars_bilateral');
        
        % now calculate the coverage
        pRF_COV = s_rmPlotCoverage_VTC(view,'prf_size',1,'fieldRange',15,'eccthresh',eThresh,'method','maximum profile','nboot',50,'normalizeRange',0,'smoothSigma',1,'cothresh',varThresh,'weight','variance explained','sigmathresh',sigThresh,'addcenters',1,'newfig',-1);
        coverage.ots1(:,:,i) = pRF_COV;
        coverage.ots1_size(i) = length(view.ROIs(end).coords);
    end
    
    mrvCleanWorkspace; clear pRF_COV view roiListNew rois roi;
    
end

% Let's save the coverage info so we don't have to recalculate

str1 = num2str(varThresh);
str2 = num2str(eThresh(2));
saveFile = fullfile(savePathData,['coverageForDensity_data_FaceWordBilateral_Vthresh' str1(end-1:end) '_Ethresh' str2]);

save(saveFile,'coverage');

%% Let's load previously saved coverage if already calculated
else
    
str1 = num2str(varThresh);
str2 = num2str(eThresh(2));
dataFile = fullfile(savePathData,['coverageForDensity_data_FaceWordBilateral_Vthresh' str1(end-1:end) '_Ethresh' str2]);

load(dataFile);

end

%% Calculate coverage area

% get indices of matched subjects
fileName = 'bi_EccVsSigma_lineData_anyHemi.mat';
dataDir = '/sni-storage/kalanit/biac2/kgs/projects/Longitudinal/FMRI/Retinotopy/results/pRF_data/EccVsSigma';

load(fullfile(dataDir,fileName));
for i = 1:numel(lineData)
    age(i,1) = lineData{1,i}(1,1).age;
    variance = lineData{1,i}(1,2).variance; % We're thresholding off V2 at 0.45 variance explained
    variance = variance(variance>=0.1);
    varexp(i,1) = nanmean(variance);
end

% Now let's define subject indices
if matchFlag
    % Just so happens the worst kids and best adults easily separate in V2,
    % so we will use our cutoff thresholds there 
    kidI = age<16 & varexp>0.45; % This gets rid of 8 worst kids
    adI  = age>18 & varexp<0.63; % This gets rid of 3 best adults
else
    kidI = age<16;
    adI  = age>18;
end

% Now we need to extract the coverage maps for kids and adults and
% count the number of subjects going into each coverage average
pfusCov = coverage.pfus; 
kNpfus =  sum(~isnan(squeeze(pfusCov(1,1,kidI))));
aNpfus =  sum(~isnan(squeeze(pfusCov(1,1,adI))));

ots1Cov = coverage.ots1; 
kNots1 =  sum(~isnan(squeeze(ots1Cov(1,1,kidI))));
aNots1 =  sum(~isnan(squeeze(ots1Cov(1,1,adI))));

% Now let's produce the average coverage plots
kpfuscov = nanmean(pfusCov(:,:,kidI),3);
apfuscov = nanmean(pfusCov(:,:,adI),3);

kots1cov = nanmean(ots1Cov(:,:,kidI),3);
aots1cov = nanmean(ots1Cov(:,:,adI),3);


% First we will exclude subjects with fewer than 10 pRFs and normalize
for i = 1:size(pfusCov,3)
    if coverage.pfus_size(i)<10
        pfusCov(:,:,i) = NaN;
    end
    pfusCov(:,:,isnan(coverage.pfus_size)) = NaN;
    
    
    if coverage.ots1_size < 10
        ots1Cov(:,:,i) = NaN;
    end
    ots1Cov(:,:,isnan(coverage.ots1_size)) = NaN;
end

% Now threshhold density in each subject
for i = 1:size(pfusCov,3)
    faceTmp = pfusCov(:,:,i); 
    faceTmp(faceTmp>=thresh)=1; 
    faceTmp(faceTmp<=thresh)=0;
    faceCoverage(i) = sum(sum(faceTmp));
    
    wordTmp = ots1Cov(:,:,i); 
    wordTmp(wordTmp>=thresh)=1; 
    wordTmp(wordTmp<=thresh)=0;
    wordCoverage(i) = sum(sum(wordTmp));
end

% We can convert the normalized coverage to square degrees-visual-angle
for i = 1:length(faceCoverage)
faceFoV(i) = area*((faceCoverage(i))/(128^2));
wordFoV(i) = area*((wordCoverage(i))/(128^2));
end

% Calculate density
densityFace = coverage.pfus_size ./ faceFoV;
% get rid of infinity values
densityFace(densityFace>500)=NaN;
densityWord = coverage.ots1_size ./ wordFoV;

savefile = fullfile(savePathData,'face_word_coverage_calculation_corrected.mat');
save(savefile,'faceFoV','wordFoV')

% do an ANOVA to see if there is a main effect of age group
data = [faceFoV wordFoV];

roiLabel = [zeros(1,length(faceFoV)) ones(1,length(wordFoV))];
kidLabel = ones(1,length(age(age<18)))+1; adultLabel = ones(1,length(age(age>18)))+2;
ageLabel = [kidLabel adultLabel kidLabel adultLabel];
groups = [roiLabel' ageLabel'];

p = anovan(data,groups,'interaction');
%% Plot the coverage results
kCovF = faceFoV(kidI); 
aCovF = faceFoV(adI);  
kCovF(end:end+(length(aCovF)-length(kCovF))) = NaN;
kCovW = (wordFoV(kidI)); 
aCovW = (wordFoV(adI)); 
kCovW(end:end+(length(aCovW)-length(kCovW))) = NaN;

allCov = [kCovW' aCovW' kCovF' aCovF'];
CovMeans = [nanmean(kCovF) nanmean(aCovF) nanmean(kCovW) nanmean(aCovW)];

% try box plots
g = figure('Position',[100 100 900 600],'Color',[ 1 1 1]); set(gca,'ylim',[0 7]);

subplot_tight(1,1,1,[0.06 0.08]);
boxplot(allCov,'Symbol','o','Widths',0.5);
set(gca,'xticklabel',{'Child','Adult','Child','Adult'},'FontSize',18)
ylabel('Coverage area (dva^2)')

savefile2 = fullfile(savePathFig,'face_word_coverage_boxplot_bilateral.fig')
savefig(savefile2)
