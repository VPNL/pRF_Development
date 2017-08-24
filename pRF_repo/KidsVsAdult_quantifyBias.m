% KidsVsAdults_quantifyBias.m
%
% This script will load the eyetracking data from the behavioral
% recognition experiment and quantify the direction of bias of children
% fixations compared to adults. 

clear
curdir = pwd; outputdir = fullfile(curdir,'output'); if ~exist(outputdir), mkdir(outputdir); end

%% First do faces
thresh = 0.3;
subs = {'kw11_3_eyedata_processed.mat'...
        'AOK08_3_eyedata_processed.mat'...
        'nw10_eyedata_processed.mat'...
        'os13_eyedata_processed.mat'...
        'LL11_3_eyedata_processed.mat'...
        'RBJ09_3_eyedata_processed.mat'...
        'RJM11_3_eyedata_processed.mat'...
        'CLC06_eyedata_processed.mat'...
        'SERA12_3_eyedata_processed.mat'...
        'GEJA_3_eyedata_processed.mat'...
        'RHSA08_3_eyedata_processed.mat'...
        'AW06_3_eyedata_processed.mat'...
        'jh22_eyedata_processed.mat'...
        'MW23_3_eyedata_processed.mat'...
        'ml23_eyedata_processed.mat'...
        'sl23_eyedata_processed.mat'...
        'CR24_3_eyedata_processed.mat'...
        'cb24_eyedata_processed.mat'...
        'LB23_3_eyedata_processed.mat'...
        'AD26_eyedata_processed.mat'...
        'MMC27_eyedata_processed.mat'...
        'MC26_eyedata_processed.mat'...
        'MH28_3_eyedata_processed.mat'};

ageFlag = [ones(12,1); zeros(11,1)];
    
% Define the data smoothing filter here for 2D convolution 
gauss = gauss2mf([1:1:50],[18.75 25 18.75 25]);
gaussX = repmat(gauss, 50, 1);
gy = gauss'; gaussY = repmat(gy, 1, 50);
filter = gaussX .* gaussY;

% Load data from the raw folder
dataDir = fullfile(curdir,'behavioral_data','processed');
timeDir = fullfile(curdir,'behavioral_data','stim_times');
load(fullfile(timeDir,'allSubNames.mat'));
load(fullfile(timeDir,'allFaceTimes.mat'));

fieldA = zeros(768, 1024, 16);
fieldK = zeros(768, 1024, 16);

numFixa = nan(23,16);
fixDura = nan(23,16);

for s = 1:length(subs)
load(fullfile(dataDir,subs{s}));

eyetime = data.eyetime;
xPoints = data.eyex;
yPoints = data.eyey; yPoints = 768 - yPoints; % We have to flip y

% Now we will loop through all 16 face stimuli for this subject and store
% their fixation points into a field matrix for each stim
for stim = 1:16
    % Indices corresponding to a stimulus' presentation
    time = eyetime>=allFaceTimes(s,1,stim) & eyetime<=allFaceTimes(s,2,stim);
    numFixa(s,stim) = allFaceTimes(s,3,stim);
    fixDura(s,stim) = allFaceTimes(s,4,stim);
    % Now put the points in the field on which the subject fixated into the
    % storage matrices depending on subject age
    x = xPoints(time);
    y = yPoints(time);
    x = round(x,0);
    y = round(y,0); 
    bad = x>1024 | x<=0; bady = y>768 | y<=0;
    bad(bady) = 1;
    x(bad) = NaN;
    y(bad) = NaN;
    
    if ageFlag(s)==0
        for fixation = 1:length(x)
            if isnan(x(fixation))
                continue
            end
            fieldA(y(fixation), x(fixation), stim) = fieldA(y(fixation), x(fixation), stim) + 1;
        end
    elseif ageFlag(s)==1
        for fixation = 1:length(x)
            if isnan(x(fixation))
                continue
            end
            fieldK(y(fixation), x(fixation), stim) = fieldK(y(fixation), x(fixation), stim) + 1;
        end 
    end
    
end
    
end


% Smooth the data a little (params defined above)
for f = 1:size(fieldA,3)
fieldSmoothA(:,:,f) = conv2(fieldA(:,:,f),filter);
fieldSmoothK(:,:,f) = conv2(fieldK(:,:,f),filter);
end

% Normalize by the max
for f = 1:size(fieldSmoothA,3)
fieldNormA(:,:,f) = fieldSmoothA(:,:,f) ./ max(max(fieldSmoothA(:,:,f)));
fieldNormK(:,:,f) = fieldSmoothK(:,:,f) ./ max(max(fieldSmoothK(:,:,f)));
end


% Threshold fixation density
fieldCutA = fieldNormA;  fieldCutA(fieldCutA<thresh)=0; fieldCutA(fieldCutA>=thresh)=1;
fieldCutK = fieldNormK;  fieldCutK(fieldCutK<thresh)=0; fieldCutK(fieldCutK>=thresh)=1;

% For each stimulus we calculate the center of mass of the thresholded
% fixations. 
xvals = [-536:1:536]; yvals = [-408:1:408]';
coordsX = repmat(xvals,817,1); coordsY = repmat(yvals,1,1073); coordsY = flipud(coordsY);
centersA = nan(16,2);
centersK = nan(16,2);
for zz = 1:16
    coordsXtmpK = coordsX; coordsXtmpK(fieldCutK(:,:,zz)==0)=NaN; 
    coordsYtmpK = coordsY; coordsYtmpK(fieldCutK(:,:,zz)==0)=NaN;
    avXK = coordsXtmpK .* fieldCutK(:,:,zz); avXK = nanmean(avXK,1); avXK = nanmean(avXK);
    avYK = coordsYtmpK .* fieldCutK(:,:,zz); avYK = nanmean(avYK,2); avYK = nanmean(avYK);
    centersK(zz,:) = [avXK avYK];
    
    coordsXtmpA = coordsX; coordsXtmpA(fieldCutA(:,:,zz)==0)=NaN; 
    coordsYtmpA = coordsY; coordsYtmpA(fieldCutA(:,:,zz)==0)=NaN;
    avXA = coordsXtmpA .* fieldCutA(:,:,zz); avXA = nanmean(avXA,1); avXA = nanmean(avXA);
    avYA = coordsYtmpA .* fieldCutA(:,:,zz); avYA = nanmean(avYA,2); avYA = nanmean(avYA);
    centersA(zz,:) = [avXA avYA];
end

% Now we will calculate vectors from adult centers to kid centers
% For faces our hypothesis is that kids fixate more to the right visual
% field so we will subtract from their X/Y value the adult values, so resulting
% values that are positive support our hypothesis:
centersDiff = centersK - centersA;
for i = 1:16
    distances(i) = sqrt((centersDiff(i,1))^2 + (centersDiff(i,2))^2);
end

% Now let's plot vectors on a grid that has an additional feature. First,
% we will include a gaussian gray central coloration that predicts where
% fixation centers should be if kids were "poor" fixaters and just fixated
% randomly outside of the adult centers. They should thus be distributed
% in a gaussian pattern around the center with a mean center of zero,zero.

f = figure('position',[10 10 500 500],'color','w');
maxval = 2 * max(distances); maxval = ceil(maxval); maxvalF=maxval;
gaus = fspecial('gaussian',[maxval maxval], maxval/6);
invgray = colormap(gray(100)); invgray = flipud(invgray);
imagesc(gaus); colormap(invgray); hold on; 
plot([0 maxval],[maxval/2 maxval/2],'-','linewidth',2,'color',[0.5 0.5 0.5]); 
plot([maxval/2 maxval/2],[0 maxval],'-','linewidth',2,'color',[0.5 0.5 0.5])

% Now plot difference vectors
for i = 1:16
    plot([maxval/2 centersDiff(i,1)+maxval/2],[maxval/2 centersDiff(i,2)+maxval/2],'r-','linewidth',2);
    scatter(centersDiff(:,1)+maxval/2,centersDiff(:,2)+maxval/2,'filled','MarkerFaceColor','r','MarkerEdgeColor','r');
    end
%set(gca,'xtick',[],'ytick',[])
% convert pixels into degrees visual angle
dva = (4.5*maxval)/200;
set(gca,'xtick',[1 maxval/2 maxval],'xticklabel',{num2str(-dva/2),'0', num2str(dva/2)})
set(gca,'ytick',[1 maxval/2 maxval],'yticklabel',{num2str(dva/2),'0', num2str(-dva/2)})
set(gca,'fontsize',16)
xlabel('Degrees of visual angle')
title('Child fixation vectors relative to adults')

centersFace = centersDiff;

%% Now do words
clearvars -except maxvalF f centersFace curdir outputdir
thresh = 0.3;
subs = {'kw11_3_eyedata_processed.mat'...
        'AOK08_3_eyedata_processed.mat'...
        'nw10_eyedata_processed.mat'...
        'os13_eyedata_processed.mat'...
        'LL11_3_eyedata_processed.mat'...
        'RBJ09_3_eyedata_processed.mat'...
        'RJM11_3_eyedata_processed.mat'...
        'CLC06_eyedata_processed.mat'...
        'SERA12_3_eyedata_processed.mat'...
        'GEJA_3_eyedata_processed.mat'...
        'RHSA08_3_eyedata_processed.mat'...
        'AW06_3_eyedata_processed.mat'...
        'jh22_eyedata_processed.mat'...
        'MW23_3_eyedata_processed.mat'...
        'ml23_eyedata_processed.mat'...
        'sl23_eyedata_processed.mat'...
        'CR24_3_eyedata_processed.mat'...
        'cb24_eyedata_processed.mat'...
        'LB23_3_eyedata_processed.mat'...
        'AD26_eyedata_processed.mat'...
        'MMC27_eyedata_processed.mat'...
        'MC26_eyedata_processed.mat'...
        'MH28_3_eyedata_processed.mat'};

ageFlag = [ones(12,1); zeros(11,1)];
    
% Define the data smoothing filter here for 2D convolution 
gauss = gauss2mf([1:1:50],[18.75 25 18.75 25]);
gaussX = repmat(gauss, 50, 1);
gy = gauss'; gaussY = repmat(gy, 1, 50);
filter = gaussX .* gaussY;


% Load data from the raw folder
dataDir = fullfile(curdir,'behavioral_data','processed');
timeDir = fullfile(curdir,'behavioral_data','stim_times');
load(fullfile(timeDir,'allSubNames.mat'));
load(fullfile(timeDir,'allWordTimes.mat'));

fieldA = zeros(768, 1024, 16);
fieldK = zeros(768, 1024, 16);

numFixa = nan(23,16);
fixDura = nan(23,16);

for s = 1:length(subs)
load(fullfile(dataDir,subs{s}));

eyetime = data.eyetime;
xPoints = data.eyex;
yPoints = data.eyey; yPoints = 768 - yPoints; % We have to flip y

% Now we will loop through all 16 face stimuli for this subject and store
% their fixation points into a field matrix for each stim
for stim = 1:16
    % Indices corresponding to a stimulus' presentation
    time = eyetime>=allWordTimes(s,1,stim) & eyetime<=allWordTimes(s,2,stim);
    numFixa(s,stim) = allWordTimes(s,3,stim);
    fixDura(s,stim) = allWordTimes(s,4,stim);
    % Now put the points in the field on which the subject fixated into the
    % storage matrices depending on subject age
    x = xPoints(time);
    y = yPoints(time);
    x = round(x,0);
    y = round(y,0); 
    bad = x>1024 | x<=0; bady = y>768 | y<=0;
    bad(bady) = 1;
    x(bad) = NaN;
    y(bad) = NaN;
    
    if ageFlag(s)==0
        for fixation = 1:length(x)
            if isnan(x(fixation))
                continue
            end
            fieldA(y(fixation), x(fixation), stim) = fieldA(y(fixation), x(fixation), stim) + 1;
        end
    elseif ageFlag(s)==1
        for fixation = 1:length(x)
            if isnan(x(fixation))
                continue
            end
            fieldK(y(fixation), x(fixation), stim) = fieldK(y(fixation), x(fixation), stim) + 1;
        end 
    end
    
end
    
end


% Smooth the data a little (params defined above)
for f = 1:size(fieldA,3)
fieldSmoothA(:,:,f) = conv2(fieldA(:,:,f),filter);
fieldSmoothK(:,:,f) = conv2(fieldK(:,:,f),filter);
end

% Normalize by the max
for f = 1:size(fieldSmoothA,3)
fieldNormA(:,:,f) = fieldSmoothA(:,:,f) ./ max(max(fieldSmoothA(:,:,f)));
fieldNormK(:,:,f) = fieldSmoothK(:,:,f) ./ max(max(fieldSmoothK(:,:,f)));
end


% Threshold fixation density
fieldCutA = fieldNormA;  fieldCutA(fieldCutA<thresh)=0; fieldCutA(fieldCutA>=thresh)=1;
fieldCutK = fieldNormK;  fieldCutK(fieldCutK<thresh)=0; fieldCutK(fieldCutK>=thresh)=1;

% For each stimulus we calculate the center of mass of the thresholded
% fixations. 
xvals = [-536:1:536]; yvals = [-408:1:408]';
coordsX = repmat(xvals,817,1); coordsY = repmat(yvals,1,1073); coordsY = flipud(coordsY);
centersA = nan(16,2);
centersK = nan(16,2);
for zz = 1:16
    coordsXtmpK = coordsX; coordsXtmpK(fieldCutK(:,:,zz)==0)=NaN; 
    coordsYtmpK = coordsY; coordsYtmpK(fieldCutK(:,:,zz)==0)=NaN;
    avXK = coordsXtmpK .* fieldCutK(:,:,zz); avXK = nanmean(avXK,1); avXK = nanmean(avXK);
    avYK = coordsYtmpK .* fieldCutK(:,:,zz); avYK = nanmean(avYK,2); avYK = nanmean(avYK);
    centersK(zz,:) = [avXK avYK];
    
    coordsXtmpA = coordsX; coordsXtmpA(fieldCutA(:,:,zz)==0)=NaN; 
    coordsYtmpA = coordsY; coordsYtmpA(fieldCutA(:,:,zz)==0)=NaN;
    avXA = coordsXtmpA .* fieldCutA(:,:,zz); avXA = nanmean(avXA,1); avXA = nanmean(avXA);
    avYA = coordsYtmpA .* fieldCutA(:,:,zz); avYA = nanmean(avYA,2); avYA = nanmean(avYA);
    centersA(zz,:) = [avXA avYA];
end

% Now we will calculate vectors from adult centers to kid centers
% For words our hypothesis is that kids fixate more to the left visual
% field so we will subtract from their X/Y value the adult values, so resulting
% x values that are negative support our hypothesis:
centersDiff = centersK - centersA;
for i = 1:16
    distances(i) = sqrt((centersDiff(i,1))^2 + (centersDiff(i,2))^2);
end

% Now let's plot vectors on a grid 
maxvalW = 2 * max(distances); maxvalW = ceil(maxvalW); maxes = [maxvalF maxvalW];
maxval = maxvalF;

% Now plot difference vectors
for i = 1:16
    plot([maxval/2 centersDiff(i,1)+maxval/2],[maxval/2 centersDiff(i,2)+maxval/2],'b-','linewidth',2);
    scatter(centersDiff(:,1)+maxval/2,centersDiff(:,2)+maxval/2,'filled','MarkerFaceColor','b','MarkerEdgeColor','b');
    %text(centersDiff(i,1)+maxval/2,centersDiff(i,2)+maxval/2,num2str(i),'fontsize',14,'color','b')
end
centersWord = centersDiff;
%% stats
% We will quantify this metric by asking if the vectors from each category
% significantly deviate away from quadrant in which the visual field
% coverage bias exists. For example, right pFus-faces in children has a VFC
% constrained almost entirely to the lower left quadrant. We thus predict
% that fixations should land outside of this quadrant. By chance 25% of
% fixations should be in the lower-left quadrant if kids fixate randomly 
% as suggested by the "poor fixaters" model. So the distance in theta of
% these points from the lower-left qudrant would be zero, and the mean of
% all others should be centered at 135 degrees, so the null distance would
% be 56.3 degrees. 

% First words
thetas = [];
centers = centersWord; centers(:,2) = -1 .* centers(:,2);
for n=1:length(centers)
    [th,rho] = cart2pol(centers(n,1),centers(n,2));
    thetas(n) = th;
end
% convert radians to degrees
thetas = 57.3 * thetas; 
%thetas = abs(thetas);
% word coverage is in lower-right quadrant in kids, so the null vector is
% 315 degrees. So we can rotate the word vectos by 45 degrees (by
% adding) and then any value over 180 degrees we will subtract from it 360,
% thereby giving us its distance to the new 0 degree line, and then it is
% this new vector of thetas that we can compare to our null
thetas =  thetas+45; 
thetas(find(thetas>180))=360-thetas(find(thetas>180));
testThetas = thetas-56.3;
[hW,pW,ciW,stW] = ttest(testThetas)

% Now faces
thetas = [];
centers = centersFace; centers(:,2) = -1 .* centers(:,2);
for n=1:length(centers)
    [th,rho] = cart2pol(centers(n,1),centers(n,2));
    thetas(n) = th;
end
thetas = 57.3 * thetas; 
% face coverage is in lower-left quadrant in kids, so the null vector is
% 225 degrees. We'll subtract it from each vector to get its angluar
% distance. 

% face coverage is in lower-left quadrant in kids, so the null vector is
% 225 degrees. So we can rotate the face vectors by 135 degrees (by
% adding) and then any value over 180 degrees we will subtract from it 360,
% thereby giving us its distance to the new 0 degree line, and then it is
% this new vector of thetas that we can compare to our null
thetas =  thetas+135;
thetas(find(thetas>180))=360-thetas(find(thetas>180));
testThetas = thetas-56.3;
[hF,pF,ciF,stF] = ttest(testThetas)

gcf;
text(2,41,['t(15) = ' num2str(stW.tstat,2) ', p=' num2str(pW,2)],'fontsize',16,'color','b')
text(2,44,['t(15) = ' num2str(stF.tstat,2) ', p=' num2str(pF,2)],'fontsize',16,'color','r')
text(38,3,'Words','color','b','fontsize',18)
text(38,5,'Faces','color','r','fontsize',18)
set(gcf,'PaperPositionMode','auto')

