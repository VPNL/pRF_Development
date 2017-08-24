% plotBarsForFixationData.m
%
% Will load and analyze the fixation data from the behavioral recognition
% memory experiment and analyze the number of fixations in children outside
% the adult-like zone, the duration, and the number of fixations. 

%% Faces
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

% Now we will loop through all 16 word stimuli for this subject and store
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
    time = eyetime>=allFaceTimes(s,1,stim) & eyetime<=allFaceTimes(s,2,stim);
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
avgRatioFace = nanmean(ratios,2); semRatioFace = (nanstd(avgRatioFace)/sqrt(12));
avFixNumFace = nanmean(numFixa,2); semFixNumFaceK = (nanstd(avFixNumFace(1:12))/sqrt(12)); semFixNumFaceA = (nanstd(avFixNumFace(13:end))/sqrt(11));
avFixDurFace = nanmean(fixDura,2); semFixDurFaceK = (nanstd(avFixDurFace(1:12))/sqrt(12)); semFixDurFaceA = (nanstd(avFixDurFace(13:end))/sqrt(11));
clearvars -except avgRatioFace semRatioFace semFixNumFaceK semFixDurFaceK avFixDurFace avFixNumFace semFixNumFaceA semFixDurFaceA
%% Words
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
avgRatioWord = nanmean(ratios,2); semRatioWord = (nanstd(avgRatioWord)/sqrt(12));
avFixNumWord = nanmean(numFixa,2); semFixNumWordK = (nanstd(avFixNumWord(1:12))/sqrt(12)); semFixNumWordA = (nanstd(avFixNumWord(13:end))/sqrt(11));
avFixDurWord = nanmean(fixDura,2); semFixDurWordK = (nanstd(avFixDurWord(1:12))/sqrt(12)); semFixDurWordA = (nanstd(avFixDurWord(13:end))/sqrt(11));


%% plot bar graphs
f = figure('Position',[50 50 800 300]);
subplot(1,3,1)
avRatios = [nanmean(avgRatioFace) nanmean(avgRatioWord)];
semRatios = [semRatioFace semRatioWord];
bar(avRatios, 0.4); hold on; errorbar(avRatios, semRatios,'k.','LineWidth',3); set(gca,'ylim',[0 1], 'xlim', [0.5 2.5], 'fontsize', 18);
hold on; plot([1 2], [avgRatioFace avgRatioWord], 'o','MarkerEdgeColor',[0.6 0.6 0.6])
plot([0 2.5],[0.3 0.3],'k--','LineWidth',3)
set(f,'PaperPositionMode', 'auto'); box off
ylabel('Ratio of child fixations outside adult zone')
set(gca,'xticklabel',{'Faces' 'Words'});

subplot(1,3,2)
avFixNums = [nanmean(avFixNumFace(1:12)) nanmean(avFixNumFace(13:end)) nanmean(avFixNumWord(1:12)) nanmean(avFixNumWord(13:end))];
semFixNums = [semFixNumFaceK semFixNumFaceA semFixNumWordK semFixNumWordA];
bar([1 1.3 2 2.3], avFixNums, 0.6); hold on; errorbar([1 1.3 2 2.3], avFixNums, semFixNums, 'k.','LineWidth',3); set(gca,'ylim',[0 16], 'xlim', [0.5 2.7], 'fontsize', 18);
plot([1 1.3 2 2.3],[avFixNumFace(1:12) [avFixNumFace(13:end); NaN] avFixNumWord(1:12) [avFixNumWord(13:end); NaN]],'o','MarkerEdgeColor',[0.6 0.6 0.6])
set(f,'PaperPositionMode', 'auto'); box off
set(gca,'xtick',[1 1.3 2 2.3],'xticklabel',{'C' 'A' 'C' 'A'});
ylabel('# of fixations')

subplot(1,3,3)
avFixDurs = [nanmean(avFixDurFace(1:12)) nanmean(avFixDurFace(13:end)) nanmean(avFixDurWord(1:12)) nanmean(avFixDurWord(13:end))];
semFixDurs = [semFixDurFaceK semFixDurFaceA semFixDurWordK semFixDurWordA];
bar([1 1.3 2 2.3], avFixDurs, 0.6); hold on; errorbar([1 1.3 2 2.3], avFixDurs, semFixDurs,'k.','LineWidth',3); set(gca,'ylim',[200 600], 'xlim', [0.5 2.7], 'fontsize', 18);
plot([1 1.3 2 2.3],[avFixDurFace(1:12) [avFixDurFace(13:end); NaN] avFixDurWord(1:12) [avFixDurWord(13:end); NaN]],'o','MarkerEdgeColor',[0.6 0.6 0.6]);
set(f,'PaperPositionMode', 'auto'); box off
set(gca,'xtick',[1 1.3 2 2.3],'xticklabel',{'C' 'A' 'C' 'A'});
ylabel('Fixation duration (ms)')

%% plot boxplots if you so desire (uncomment if needed)

% f = figure('Position',[50 50 800 300]);
% subplot(1,3,1)
% boxplot([avgRatioFace avgRatioWord], 'Widths',0.5)
% set(gca,'FontSize',18,'ylim',[0 1])
% hold on; plot([0 2.5],[0.3 0.3],'k--','LineWidth',3)
% 
% subplot(1,3,2)
% boxplot([avFixNumFace(1:12) [avFixNumFace(13:end); NaN] avFixNumWord(1:12) [avFixNumWord(13:end); NaN]],'Widths',0.5)
% set(gca,'FontSize',18)
% 
% subplot(1,3,3)
% boxplot([avFixDurFace(1:12) [avFixDurFace(13:end); NaN] avFixDurWord(1:12) [avFixDurWord(13:end); NaN]],'Widths',0.5)
% set(gca,'FontSize',18)

