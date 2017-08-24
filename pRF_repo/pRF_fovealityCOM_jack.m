% pRF_fovealityCOM_jack.m
%
% This script will load data saved from pRF_plotAveragedCoverage_VTC and
% calculate the distance of the coverage's center of mass from the center 
% of the visual field. Produces Figure 4E.
%
% JG 06/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
curdir = pwd; outputdir = fullfile(curdir,'output'); if ~exist(outputdir), mkdir(outputdir); end

% match subjects for variance?
matchFlag = true;

% Let's load the saved coverage maps file:
savePathData = fullfile(curdir,'voxel_data');
fileName = 'coverage_data_VTC_Vthresh05_Ethresh12.mat';
load(fullfile(savePathData,fileName));

% Get subject indices
if matchFlag
    % These subjects are matched for variance explained by the pRF model in
    % V1. These subjects are the those present in the main figures.
    load(fullfile(curdir,'voxel_data','varMatched_indices.mat'));
else
    kidI = zeros(1,53); kidI(1:26)=1;
    adI  = zeros(1,53); adI(28:end)=1;
end

% Now extract coverage for the individual ROIs:
rpfusCov = coverage.rpfus; 
lpfusCov = coverage.lpfus; 
rots1Cov = coverage.rots1; 
lots1Cov = coverage.lots1;   

%% right pFus-faces
kN =  sum(~isnan(squeeze(rpfusCov(1,1,kidI))));
aN =  sum(~isnan(squeeze(rpfusCov(1,1,adI))));

% Let's get subjects that have coverage
covI = ~isnan(squeeze(rpfusCov(1,1,:)));
for q = 1:length(kidI), if kidI(q)==1 && covI(q) ==1, kInd(q)=1; else kInd(q)=0; end, end
for q = 1:length(adI), if adI(q)==1 && covI(q) ==1, aInd(q)=1; else aInd(q)=0; end, end    

kInd = logical(kInd); aInd = logical(aInd);

kCov = rpfusCov(:,:,kInd);
aCov = rpfusCov(:,:,aInd);

% now below we're going to jackknife each time drawing kN or aN from each
% distribution of kCov or aCov:

    coordsA = NaN(128*128,2,aN);
    coordsK = NaN(128*128,2,kN);

for boot = 1:kN
    kboot = [1:1:kN]; kboot(boot)=[];
    k = nanmean(kCov(:,:,kboot),3);
   
    totalWeightK(boot) = sum(sum(k));

    ind = 0;
    for i = 1:128
        for j = 1:128
            ind = ind+1;
            coordsK(ind,1,boot) = k(i,j)*i;
            coordsK(ind,2,boot) = k(i,j)*j;
        end
    end
end

    kSum = sum(coordsK,1);

    for z = 1:kN
        kCOM(1,:,z) = kSum(1,:,z) ./ totalWeightK(z);
    end
    
    % Now let's find radial distance in DVA from center/fovea. Every 9.15
    % pixels is a degree of visual angle (128pixels/14 degree diameter)
    fovea = [64 64];
    for z = 1:kN
        kDist(1,:,z) = (kCOM(1,:,z) - fovea) ./ 9.15;
    end
    
    % Now let's convert to radial distance:
    for z = 1:kN
        coordK = kDist(1,:,z);
        kDistRadial(z) = sqrt(coordK(1)^2 + coordK(2)^2);
    end

% Rename for plotting:;
kDist_rpfus_avg = nanmean(kDistRadial); kDist_rpfus_ste = 1.96*(nanstd(kDistRadial)/sqrt(kN));
kDist_rpfus = kDistRadial;
clear kDistRadial;

% Now do adults
for boot = 1:aN
    aboot = [1:1:aN]; aboot(boot)=[];
    a = nanmean(aCov(:,:,aboot),3);
   
    totalWeightA(boot) = sum(sum(a));

    ind = 0;
    for i = 1:128
        for j = 1:128
            ind = ind+1;
            coordsA(ind,1,boot) = a(i,j)*i;
            coordsA(ind,2,boot) = a(i,j)*j;
        end
    end
end

    aSum = sum(coordsA,1);

    for z = 1:aN
        aCOM(1,:,z) = aSum(1,:,z) ./ totalWeightA(z);
    end
    
    % Now let's find radial distance in DVA from center/fovea. Every 9.15
    % pixels is a degree of visual angle (128pixels/14 degree diameter)
    fovea = [64 64];
    for z = 1:aN
        aDist(1,:,z) = (aCOM(1,:,z) - fovea) ./ 9.15;
    end
    
    % Now let's convert to radial distance:
    for z = 1:aN
        coordA = aDist(1,:,z);
        aDistRadial(z) = sqrt(coordA(1)^2 + coordA(2)^2);
    end

% Rename for plotting:;
aDist_rpfus_avg = nanmean(aDistRadial); aDist_rpfus_ste = 1.96*(nanstd(aDistRadial)/sqrt(aN));
aDist_rpfus = aDistRadial;
clear aDistRadial;

%% left pFus-faces
kN =  sum(~isnan(squeeze(lpfusCov(1,1,kidI))));
aN =  sum(~isnan(squeeze(lpfusCov(1,1,adI))));

% Let's get subjects that have coverage
covI = ~isnan(squeeze(lpfusCov(1,1,:)));
for q = 1:length(kidI), if kidI(q)==1 && covI(q) ==1, kInd(q)=1; else kInd(q)=0; end, end
for q = 1:length(adI), if adI(q)==1 && covI(q) ==1, aInd(q)=1; else aInd(q)=0; end, end    

kInd = logical(kInd); aInd = logical(aInd);

kCov = lpfusCov(:,:,kInd);
aCov = lpfusCov(:,:,aInd);

% now below we're going to jackknife each time drawing kN or aN from each
% distribution of kCov or aCov:

    coordsA = NaN(128*128,2,aN);
    coordsK = NaN(128*128,2,kN);

for boot = 1:kN
    kboot = [1:1:kN]; kboot(boot)=[];
    k = nanmean(kCov(:,:,kboot),3);
   
    totalWeightK(boot) = sum(sum(k));

    ind = 0;
    for i = 1:128
        for j = 1:128
            ind = ind+1;
            coordsK(ind,1,boot) = k(i,j)*i;
            coordsK(ind,2,boot) = k(i,j)*j;
        end
    end
end

    kSum = sum(coordsK,1);

    for z = 1:kN
        kCOM(1,:,z) = kSum(1,:,z) ./ totalWeightK(z);
    end
    
    % Now let's find radial distance in DVA from center/fovea. Every 9.15
    % pixels is a degree of visual angle (128pixels/14 degree diameter)
    fovea = [64 64];
    for z = 1:kN
        kDist(1,:,z) = (kCOM(1,:,z) - fovea) ./ 9.15;
    end
    
    % Now let's convert to radial distance:
    for z = 1:kN
        coordK = kDist(1,:,z);
        kDistRadial(z) = sqrt(coordK(1)^2 + coordK(2)^2);
    end

% Rename for plotting:;
kDist_lpfus_avg = nanmean(kDistRadial); kDist_lpfus_ste = 1.96*(nanstd(kDistRadial)/sqrt(kN));
kDist_lpfus = kDistRadial;
clear kDistRadial;

% Now do adults
for boot = 1:aN
    aboot = [1:1:aN]; aboot(boot)=[];
    a = nanmean(aCov(:,:,aboot),3);
    
    totalWeightA(boot) = sum(sum(a));

    ind = 0;
    for i = 1:128
        for j = 1:128
            ind = ind+1;
            coordsA(ind,1,boot) = a(i,j)*i;
            coordsA(ind,2,boot) = a(i,j)*j;
        end
    end
end

    aSum = sum(coordsA,1);

    for z = 1:aN
        aCOM(1,:,z) = aSum(1,:,z) ./ totalWeightA(z);
    end
    
    % Now let's find radial distance in DVA from center/fovea. Every 9.15
    % pixels is a degree of visual angle (128pixels/14 degree diameter)
    fovea = [64 64];
    for z = 1:aN
        aDist(1,:,z) = (aCOM(1,:,z) - fovea) ./ 9.15;
    end
    
    % Now let's convert to radial distance:
    for z = 1:aN
        coordA = aDist(1,:,z);
        aDistRadial(z) = sqrt(coordA(1)^2 + coordA(2)^2);
    end

% Rename for plotting:;
aDist_lpfus_avg = nanmean(aDistRadial); aDist_lpfus_ste = 1.96*(nanstd(aDistRadial)/sqrt(aN));
aDist_lpfus = aDistRadial;
clear aDistRadial;

%% right OTS1-words
kN =  sum(~isnan(squeeze(rots1Cov(1,1,kidI))));
aN =  sum(~isnan(squeeze(rots1Cov(1,1,adI))));

% Let's get subjects that have coverage
covI = ~isnan(squeeze(rots1Cov(1,1,:)));
for q = 1:length(kidI), if kidI(q)==1 && covI(q) ==1, kInd(q)=1; else kInd(q)=0; end, end
for q = 1:length(adI), if adI(q)==1 && covI(q) ==1, aInd(q)=1; else aInd(q)=0; end, end    

kInd = logical(kInd); aInd = logical(aInd);

kCov = rots1Cov(:,:,kInd);
aCov = rots1Cov(:,:,aInd);

% now below we're going to jackknife each time drawing kN or aN from each
% distribution of kCov or aCov:

    coordsA = NaN(128*128,2,aN);
    coordsK = NaN(128*128,2,kN);

    
for boot = 1:kN
    kboot = [1:1:kN]; kboot(boot)=[];
    k = nanmean(kCov(:,:,kboot),3);
   
    totalWeightK(boot) = sum(sum(k));

    ind = 0;
    for i = 1:128
        for j = 1:128
            ind = ind+1;
            coordsK(ind,1,boot) = k(i,j)*i;
            coordsK(ind,2,boot) = k(i,j)*j;
        end
    end
end

    kSum = sum(coordsK,1);

    for z = 1:kN
        kCOM(1,:,z) = kSum(1,:,z) ./ totalWeightK(z);
    end
    
    % Now let's find radial distance in DVA from center/fovea. Every 9.15
    % pixels is a degree of visual angle (128pixels/14 degree diameter)
    fovea = [64 64];
    for z = 1:kN
        kDist(1,:,z) = (kCOM(1,:,z) - fovea) ./ 9.15;
    end
    
    % Now let's convert to radial distance:
    for z = 1:kN
        coordK = kDist(1,:,z);
        kDistRadial(z) = sqrt(coordK(1)^2 + coordK(2)^2);
    end

% Rename for plotting:;
kDist_rots_avg = nanmean(kDistRadial); kDist_rots_ste = 1.96*(nanstd(kDistRadial)/sqrt(kN));
kDist_rots = kDistRadial;
clear kDistRadial;

% Now do adults
for boot = 1:aN
    aboot = [1:1:aN]; aboot(boot)=[];
    a = nanmean(aCov(:,:,aboot),3);
   
    totalWeightA(boot) = sum(sum(a));

    ind = 0;
    for i = 1:128
        for j = 1:128
            ind = ind+1;
            coordsA(ind,1,boot) = a(i,j)*i;
            coordsA(ind,2,boot) = a(i,j)*j;
        end
    end
end

    aSum = sum(coordsA,1);

    for z = 1:aN
        aCOM(1,:,z) = aSum(1,:,z) ./ totalWeightA(z);
    end
    
    % Now let's find radial distance in DVA from center/fovea. Every 9.15
    % pixels is a degree of visual angle (128pixels/14 degree diameter)
    fovea = [64 64];
    for z = 1:aN
        aDist(1,:,z) = (aCOM(1,:,z) - fovea) ./ 9.15;
    end
    
    % Now let's convert to radial distance:
    for z = 1:aN
        coordA = aDist(1,:,z);
        aDistRadial(z) = sqrt(coordA(1)^2 + coordA(2)^2);
    end

% Rename for plotting:;
aDist_rots_avg = nanmean(aDistRadial); aDist_rots_ste = 1.96*(nanstd(aDistRadial)/sqrt(aN));
aDist_rots = aDistRadial;
clear aDistRadial    

%% left OTS1-words
kN =  sum(~isnan(squeeze(lots1Cov(1,1,kidI))));
aN =  sum(~isnan(squeeze(lots1Cov(1,1,adI))));

% Let's get subjects that have coverage
covI = ~isnan(squeeze(lots1Cov(1,1,:)));
for q = 1:length(kidI), if kidI(q)==1 && covI(q) ==1, kInd(q)=1; else kInd(q)=0; end, end
for q = 1:length(adI), if adI(q)==1 && covI(q) ==1, aInd(q)=1; else aInd(q)=0; end, end    

kInd = logical(kInd); aInd = logical(aInd);

kCov = lots1Cov(:,:,kInd);
aCov = lots1Cov(:,:,aInd);

% now below we're going to jackknife each time drawing kN or aN from each
% distribution of kCov or aCov:

    coordsA = NaN(128*128,2,aN);
    coordsK = NaN(128*128,2,kN);

    
for boot = 1:kN
    kboot = [1:1:kN]; kboot(boot)=[];
    k = nanmean(kCov(:,:,kboot),3);
   
    totalWeightK(boot) = sum(sum(k));

    ind = 0;
    for i = 1:128
        for j = 1:128
            ind = ind+1;
            coordsK(ind,1,boot) = k(i,j)*i;
            coordsK(ind,2,boot) = k(i,j)*j;
        end
    end
end

    kSum = sum(coordsK,1);

    for z = 1:kN
        kCOM(1,:,z) = kSum(1,:,z) ./ totalWeightK(z);
    end
    
    % Now let's find radial distance in DVA from center/fovea. Every 9.15
    % pixels is a degree of visual angle (128pixels/14 degree diameter)
    fovea = [64 64];
    for z = 1:kN
        kDist(1,:,z) = (kCOM(1,:,z) - fovea) ./ 9.15;
    end
    
    % Now let's convert to radial distance:
    for z = 1:kN
        coordK = kDist(1,:,z);
        kDistRadial(z) = sqrt(coordK(1)^2 + coordK(2)^2);
    end

% Rename for plotting:;
kDist_lots_avg = nanmean(kDistRadial); kDist_lots_ste = 1.96*(nanstd(kDistRadial)/sqrt(kN));
kDist_lots = kDistRadial;
clear kDistRadial;

% Now do adults
for boot = 1:aN
    aboot = [1:1:aN]; aboot(boot)=[];
    a = nanmean(aCov(:,:,aboot),3);
   
    totalWeightA(boot) = sum(sum(a));

    ind = 0;
    for i = 1:128
        for j = 1:128
            ind = ind+1;
            coordsA(ind,1,boot) = a(i,j)*i;
            coordsA(ind,2,boot) = a(i,j)*j;
        end
    end
end

    aSum = sum(coordsA,1);

    for z = 1:aN
        aCOM(1,:,z) = aSum(1,:,z) ./ totalWeightA(z);
    end
    
    % Now let's find radial distance in DVA from center/fovea. Every 9.15
    % pixels is a degree of visual angle (128pixels/14 degree diameter)
    fovea = [64 64];
    for z = 1:aN
        aDist(1,:,z) = (aCOM(1,:,z) - fovea) ./ 9.15;
    end
    
    % Now let's convert to radial distance:
    for z = 1:aN
        coordA = aDist(1,:,z);
        aDistRadial(z) = sqrt(coordA(1)^2 + coordA(2)^2);
    end

% Rename for plotting:;
aDist_lots_avg = nanmean(aDistRadial); aDist_lots_ste = 1.96*(nanstd(aDistRadial)/sqrt(aN));
aDist_lots = aDistRadial;
clear aDistRadial;

%% Make a crossover plot
f = figure('Position',[10 10 600 700],'color','w');


plot([1 2.3],[kDist_lpfus_avg aDist_lpfus_avg],'r','linewidth',5); hold on;
plot([1 2.3],[kDist_lots_avg aDist_lots_avg] ,'b','linewidth',5); box off;
errorbar([1 2.3],[kDist_lpfus_avg aDist_lpfus_avg],[kDist_lpfus_ste aDist_lpfus_ste],'k.','linewidth',5); hold on;
errorbar([1 2.3],[kDist_lots_avg aDist_lots_avg] ,[kDist_lots_ste aDist_lots_ste],'k.','linewidth',5);


plot([3.7 5],[kDist_rpfus_avg aDist_rpfus_avg],'r','linewidth',5); hold on;
plot([3.7 5], [kDist_rots_avg aDist_rots_avg] ,'b','linewidth',5); box off;
errorbar([3.7 5],[kDist_rpfus_avg aDist_rpfus_avg],[kDist_rpfus_ste aDist_rpfus_ste],'k.','linewidth',5); hold on;
errorbar([3.7 5],[kDist_rots_avg aDist_rots_avg] ,[kDist_rots_ste aDist_rots_ste],'k.','linewidth',5);

set(gca,'fontsize',18,'ylim',[0 4],'xlim',[0.5 5.5],'xtick',[1 2.3 3.7 5],'xticklabel',{'C' 'A' 'C' 'A'}); 
ylabel('Center-of-mass distance from fovea (dva)')

text(1.3,3.75,'Left','fontsize',18)
text(4,3.75,'Right','fontsize',18)

savePathFig = fullfile(curdir,'output','pRF_figures'); if ~exist(savePathFig), mkdir(savePathFig); end
saveFile = fullfile(savePathFig,'foveality_CenterOfMass_words_faces_jackknife.fig');
saveas(gcf,saveFile)

