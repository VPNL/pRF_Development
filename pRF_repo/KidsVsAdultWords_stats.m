% KidsVsAdultWords_stats.m
%
% This script will take the processed eyetracking data from the recognition
% memory behavioral experiment and plot the fixation density for children
% and adults on each stimulus. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
curdir = pwd; outputdir = fullfile(curdir,'output'); if ~exist(outputdir), mkdir(outputdir); end

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
gauss = gauss2mf([1:1:50],[18.75 15 18.75 15]);
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

% Now we will loop through all 16 word stimuli for this subject and store
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

% Find where 70% of adults fixate on the face.
% We will use this as a mask to quantify how much of child fixations are
% within this adult range and how much outside of it
fieldCut = fieldNormA; 
fieldCut(fieldCut<0.3)=0; 
fieldCut(fieldCut>=0.3)=1;


% Now we will calculate the ratio of each child's fixations inside and
% outside the adulthood fixation zone. We will take the amount of fixation
% time outside the adult zone divided by the total fixation time, such that
% 1 means all fixations were outside the adult region, and zero is
% adult-like. In this way, we can test if kids are significantly different
% from zero with a t-test. 

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
        'AW06_3_eyedata_processed.mat'};


% Where we will store how much of the kid's fixation is in the adultlike
% zone
ratios = [];

for s = 1:length(subs)

load(fullfile(dataDir,subs{s}));

eyetime = data.eyetime;
xPoints = data.eyex;
yPoints = data.eyey; yPoints = 768 - yPoints; % We have to flip y

for stim = 1:16
    
    field = zeros(768, 1024);
    
     % Indices corresponding to a stimulus' presentation
    time = eyetime>=allWordTimes(s,1,stim) & eyetime<=allWordTimes(s,2,stim);
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
    
    for fixation = 1:length(x)
            if isnan(x(fixation))
               continue
            end
        field(y(fixation), x(fixation)) = field(y(fixation), x(fixation)) + 1;
    end 
    
% Smooth the data a little (params defined above)
fieldSmooth = conv2(field,filter);
    
total = sum(sum(fieldSmooth));
out = sum(fieldSmooth(fieldCut(:,:,stim)==0));
ratios(s, stim) = out/total;
end

end

% We can average across all stimuli
avgRatio = nanmean(ratios,2);

% The "zero" ratio here means we need to subtract 0.3 from the kid numbers,
% because by definition 30% of adult fixations are outside the 70% cutoff
% zone:
zeroRatio = avgRatio - 0.3;

% Now calculate if kids are significantly different from zero (adults) in
% terms of their fixation patterns
%[h,p,ci,st] = ttest(zeroRatio)

% Let's plot
f = figure('Position',[50 50 1400 300]);
for x = 1:16
    subplot(2,8,x)
    imagesc(fieldNormK(:,:,x))
    set(gca,'ydir','Normal', 'xlim', [300 800], 'ylim', [200 600])
    set(f,'PaperPositionMode', 'auto')
end

g = figure('Position',[50 50 1400 300]);
for x = 1:16
    subplot(2,8,x)
    imagesc(fieldNormA(:,:,x))
    set(gca,'ydir','Normal', 'xlim', [300 800], 'ylim', [200 600])
    set(g,'PaperPositionMode', 'auto')
end



