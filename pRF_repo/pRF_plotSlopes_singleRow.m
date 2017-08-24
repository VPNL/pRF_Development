% pRF_compareSlopes.m
%
% This will load data structures produced from pRF_loopSigmaVsEcc.m and
% save out average linear fits between Eccentricity and Sigma comparing
% kids and adults. Produces figure 2C
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
curdir = pwd; outputdir = fullfile(curdir,'output'); if ~exist(outputdir), mkdir(outputdir); end

% match variance? don't make this false
matchFlag = true;

% Which file would you like to visualize? 
fileName = 'bi_EccVsSigma_lineData_anyHemi_vThresh_05.mat';

dataDir = fullfile(curdir,'voxel_data');
saveDir = outputdir; if ~exist(saveDir), mkdir(saveDir); end

load(fullfile(dataDir,fileName));

% Set up data structures
slope = NaN(numel(lineData),12);
intercept = NaN(numel(lineData),12);
age = [];
varexp = NaN(numel(lineData),12);

for i = 1:numel(lineData)
    age(i,1) = lineData{1,i}(1,1).age;
    
    for m = 1:numel(lineData{1,i})
        
        
        % If the no line could be fit in the roi, then the line will have a
        % NaN, but we'll change it to [NaN NaN] because the code expects a
        % two unit vector
        if isnan(lineData{1,i}(1,m).line)
            lineData{1,i}(1,m).line = [NaN NaN];
        end
        
        
        if strmatch('V1',lineData{1,i}(1,m).roi(end-1:end))
            variance = lineData{1,i}(1,m).variance;
            variance = variance(variance>=0.1);
            line = lineData{1,i}(1,m).line;
            if lineData{1,i}(1,1).age < 18
                varexp(i,1) = nanmean(variance);
                slope(i,1)  = line(1);
                intercept(i,1) = line(2);
            else
                varexp(i,2) = nanmean(variance);
                slope(i,2)  = line(1);
                intercept(i,2) = line(2);
            end
        
        elseif strmatch('V2',lineData{1,i}(1,m).roi(end-1:end))
            variance = lineData{1,i}(1,m).variance;
            variance = variance(variance>=0.1);
            line = lineData{1,i}(1,m).line;
            if lineData{1,i}(1,1).age < 18
                varexp(i,3) = nanmean(variance);
                slope(i,3)  = line(1);
                intercept(i,3) = line(2);
            else
                varexp(i,4) = nanmean(variance);
                slope(i,4)  = line(1);
                intercept(i,4) = line(2);
            end
          
        elseif strmatch('V3',lineData{1,i}(1,m).roi(end-1:end))
            variance = lineData{1,i}(1,m).variance;
            variance = variance(variance>=0.1);
            line = lineData{1,i}(1,m).line;
            if lineData{1,i}(1,1).age < 18
                varexp(i,5) = nanmean(variance);
                slope(i,5)  = line(1);
                intercept(i,5) = line(2);
            else
                varexp(i,6) = nanmean(variance);
                slope(i,6)  = line(1);
                intercept(i,6) = line(2);
            end
            
        elseif strmatch('V4',lineData{1,i}(1,m).roi(end-1:end))
            variance = lineData{1,i}(1,m).variance;
            variance = variance(variance>=0.1);
            line = lineData{1,i}(1,m).line;
            if lineData{1,i}(1,1).age < 18
                varexp(i,7) = nanmean(variance);
                slope(i,7)  = line(1);
                intercept(i,7) = line(2);
            else
                varexp(i,8) = nanmean(variance);
                slope(i,8)  = line(1);
                intercept(i,8) = line(2);
            end
            
        elseif strmatch('VO1',lineData{1,i}(1,m).roi(end-2:end))
            variance = lineData{1,i}(1,m).variance;
            variance = variance(variance>=0.1);
            line = lineData{1,i}(1,m).line;
            if lineData{1,i}(1,1).age < 18
                varexp(i,9) = nanmean(variance);
                slope(i,9)  = line(1);
                intercept(i,9) = line(2);
            else
                varexp(i,10) = nanmean(variance);
                slope(i,10)  = line(1);
                intercept(i,10) = line(2);
            end
  end
    end
end

% We will do a sanity check and remove anyone with negative slope fits:
varexp(slope<0)=NaN;
intercept(slope<0)=NaN;
slope(slope<0)=NaN;

% And remove negative intercept fits:
varexp(intercept<0)=NaN;
slope(intercept<0)=NaN;
intercept(intercept<0)=NaN;

% Now we will plot lines and patches of STE

%%%%%%%% plot ecc vs. sigma for V1 
if matchFlag
    % These subjects are matched for variance explained by the pRF model in
    % V1. These subjects are the those present in the main figures.
    load(fullfile(curdir,'voxel_data','varMatched_indices.mat'));
else
    kidI = zeros(1,53); kidI(1:26)=1;
    adI  = zeros(1,53); adI(28:end)=1;
end

xPoints = [0:0.5:7];

%% First plot adults

    mV1  = slope(adI,2);
    intV1  = intercept(adI,2); 
    colorV1 = [153 0 0]; colorV1=colorV1./255;

    mV2  = slope(adI,4);
    intV2  = intercept(adI,4);
    colorV2 = [153 76 0]; colorV2=colorV2./255;
    
    mV3  = slope(adI,6);
    intV3  = intercept(adI,6);
    colorV3 = [153 153 0]; colorV3=colorV3./255;

    mV4  = slope(adI,8);
    intV4  = intercept(adI,8);
    colorV4 = [76 153 0]; colorV4=colorV4./255;

    mVO1  = slope(adI,10);
    intVO1  = intercept(adI,10);
    colorVO1 = [0 153 153]; colorVO1=colorVO1./255;

% get slope (m) and intercept (i)
mV1(isnan(mV1))=[]; intV1(isnan(intV1))=[];
mMeanV1 = mean(mV1); intMeanV1 = mean(intV1);
mSTEV1 = std(mV1)/sqrt(length(mV1));
intSTEV1 = std(intV1)/sqrt(length(intV1));
yPointsV1 = (mMeanV1 .* xPoints) + intMeanV1;

mV2(isnan(mV2))=[]; intV2(isnan(intV2))=[];
mMeanV2 = mean(mV2); intMeanV2 = mean(intV2);
mSTEV2 = std(mV2)/sqrt(length(mV2));
intSTEV2 = std(intV2)/sqrt(length(intV2));
yPointsV2 = (mMeanV2 .* xPoints) + intMeanV2;

mV3(isnan(mV3))=[]; intV3(isnan(intV3))=[];
mMeanV3 = mean(mV3); intMeanV3 = mean(intV3);
mSTEV3 = std(mV3)/sqrt(length(mV3));
intSTEV3 = std(intV3)/sqrt(length(intV3));
yPointsV3 = (mMeanV3 .* xPoints) + intMeanV3;

mV4(isnan(mV4))=[]; intV4(isnan(intV4))=[];
mMeanV4 = mean(mV4); intMeanV4 = mean(intV4);
mSTEV4 = std(mV4)/sqrt(length(mV4));
intSTEV4 = std(intV4)/sqrt(length(intV4));
yPointsV4 = (mMeanV4 .* xPoints) + intMeanV4;

mVO1(isnan(mVO1))=[]; intVO1(isnan(intVO1))=[];
mMeanVO1 = mean(mVO1); intMeanVO1 = mean(intVO1);
mSTEVO1 = std(mVO1)/sqrt(length(mVO1));
intSTEVO1 = std(intVO1)/sqrt(length(intVO1));
yPointsVO1 = (mMeanVO1 .* xPoints) + intMeanVO1;

% set up patch error information
% upper and lower limits
UpperSlopeV1   = mMeanV1 + mSTEV1;
LowerSlopeV1   = mMeanV1 - mSTEV1;
yPointsUpperV1 = (UpperSlopeV1 .* xPoints) + (intMeanV1 + intSTEV1);
yPointsLowerV1 = (LowerSlopeV1 .* xPoints) + (intMeanV1 - intSTEV1);

UpperSlopeV2   = mMeanV2 + mSTEV2;
LowerSlopeV2   = mMeanV2 - mSTEV2;
yPointsUpperV2 = (UpperSlopeV2 .* xPoints) + (intMeanV2 + intSTEV2);
yPointsLowerV2 = (LowerSlopeV2 .* xPoints) + (intMeanV2 - intSTEV2);

UpperSlopeV3   = mMeanV3 + mSTEV3;
LowerSlopeV3   = mMeanV3 - mSTEV3;
yPointsUpperV3 = (UpperSlopeV3 .* xPoints) + (intMeanV3 + intSTEV3);
yPointsLowerV3 = (LowerSlopeV3 .* xPoints) + (intMeanV3 - intSTEV3);

UpperSlopeV4   = mMeanV4 + mSTEV4;
LowerSlopeV4   = mMeanV4 - mSTEV4;
yPointsUpperV4 = (UpperSlopeV4 .* xPoints) + (intMeanV4 + intSTEV4);
yPointsLowerV4 = (LowerSlopeV4 .* xPoints) + (intMeanV4 - intSTEV4);

UpperSlopeVO1   = mMeanVO1 + mSTEVO1;
LowerSlopeVO1   = mMeanVO1 - mSTEVO1;
yPointsUpperVO1 = (UpperSlopeVO1 .* xPoints) + (intMeanVO1 + intSTEVO1);
yPointsLowerVO1 = (LowerSlopeVO1 .* xPoints) + (intMeanVO1 - intSTEVO1);

% Now let's plot
f = figure('Position',[100 100 1550 475],'Color','w');
subplot_tight(1,5,1,[0.08, 0.02]); 
set(gca,'xlim',[0 7],'ylim',[0 7],'box','off','ycolor','k','tickdir','out','fontsize',18);

%%%%%%%%%%% Plot V1 information %%%%%%%%%%%
vV1 = [0 yPointsLowerV1(1); xPoints(end) yPointsLowerV1(end); xPoints(end) yPointsUpperV1(end); 0 yPointsUpperV1(1)];
f = [1 2 3 4];
patch('Vertices',vV1,'Faces',f,'FaceColor',colorV1,'FaceAlpha',0.4,'EdgeAlpha',0)
hold on; axis square;
% Now plot kids
mV1 = slope(kidI,1); 
intV1 = intercept(kidI,1); 
colorV1 = [255 51 51]; colorV1=colorV1./255;
mV1(isnan(mV1))=[]; intV1(isnan(intV1))=[];
mMeanV1 = mean(mV1); intMeanV1 = mean(intV1);
mSTEV1 = std(mV1)/sqrt(length(mV1));
intSTEV1 = std(intV1)/sqrt(length(intV1));
yPointsV1 = (mMeanV1 .* xPoints) + intMeanV1;
UpperSlopeV1   = mMeanV1 + mSTEV1;
LowerSlopeV1   = mMeanV1 - mSTEV1;
yPointsUpperV1 = (UpperSlopeV1 .* xPoints) + (intMeanV1 + intSTEV1);
yPointsLowerV1 = (LowerSlopeV1 .* xPoints) + (intMeanV1 - intSTEV1);
vV1 = [0 yPointsLowerV1(1); xPoints(end) yPointsLowerV1(end); xPoints(end) yPointsUpperV1(end); 0 yPointsUpperV1(1)];
f = [1 2 3 4];
patch('Vertices',vV1,'Faces',f,'FaceColor',colorV1,'FaceAlpha',0.4,'EdgeAlpha',0)
plot(xPoints,yPointsV1,'Color',colorV1,'LineWidth',2)
xlabel('Eccentricity (dva)','fontsize',18);
ylabel('pRF size (dva)','fontsize',18)
myTit = title('V1','fontsize',36); set(myTit,'interpreter','none');
set(myTit,'Position',[5 5.5])
set(gca,'fontsize',18);

%%%%%%%%%%% Plot V2 %%%%%%%%%%%%
subplot_tight(1,5,2,[0.08, 0.02]); 
vV2 = [0 yPointsLowerV2(1); xPoints(end) yPointsLowerV2(end); xPoints(end) yPointsUpperV2(end); 0 yPointsUpperV2(1)];
patch('Vertices',vV2,'Faces',f,'FaceColor',colorV2,'FaceAlpha',0.4,'EdgeAlpha',0)
hold on; axis square;
set(gca,'xlim',[0 7],'ylim',[0 7],'box','off','ycolor','w','tickdir','out','fontsize',14);
plot(xPoints,yPointsV2,'Color',colorV2,'LineWidth',2)
% Now plot kids
mV2 = slope(kidI,3);
intV2 = intercept(kidI,3);
colorV2 = [255 153 51]; colorV2=colorV2./255;
mV2(isnan(mV2))=[]; intV2(isnan(intV2))=[];
mMeanV2 = mean(mV2); intMeanV2 = mean(intV2);
mSTEV2 = std(mV2)/sqrt(length(mV2));
intSTEV2 = std(intV2)/sqrt(length(intV2));
yPointsV2 = (mMeanV2 .* xPoints) + intMeanV2;
UpperSlopeV2   = mMeanV2 + mSTEV2;
LowerSlopeV2   = mMeanV2 - mSTEV2;
yPointsUpperV2 = (UpperSlopeV2 .* xPoints) + (intMeanV2 + intSTEV2);
yPointsLowerV2 = (LowerSlopeV2 .* xPoints) + (intMeanV2 - intSTEV2);
vV2 = [0 yPointsLowerV2(1); xPoints(end) yPointsLowerV2(end); xPoints(end) yPointsUpperV2(end); 0 yPointsUpperV2(1)];
patch('Vertices',vV2,'Faces',f,'FaceColor',colorV2,'FaceAlpha',0.4,'EdgeAlpha',0)
plot(xPoints,yPointsV2,'Color',colorV2,'LineWidth',2)
myTit = title('V2','fontsize',36); set(myTit,'interpreter','none');
set(myTit,'Position',[5 5.5])
set(gca,'fontsize',18);

%%%%%%%%%%% Plot V3 %%%%%%%%%%%
subplot_tight(1,5,3,[0.08, 0.02]);
vV3 = [0 yPointsLowerV3(1); xPoints(end) yPointsLowerV3(end); xPoints(end) yPointsUpperV3(end); 0 yPointsUpperV3(1)];
patch('Vertices',vV3,'Faces',f,'FaceColor',colorV3,'FaceAlpha',0.4,'EdgeAlpha',0)
hold on; axis square; 
set(gca,'xlim',[0 7],'ylim',[0 7],'box','off','ycolor','w','tickdir','out','fontsize',14);
plot(xPoints,yPointsV3,'Color',colorV3,'LineWidth',2)
% now plot kids
mV3 = slope(kidI,5); 
colorV3 = [255 255 51]; colorV3=colorV3./255;
mV3(isnan(mV3))=[]; intV3(isnan(intV3))=[];
mMeanV3 = mean(mV3); intMeanV3 = mean(intV3);
mSTEV3 = std(mV3)/sqrt(length(mV3));
intSTEV3 = std(intV3)/sqrt(length(intV3));
yPointsV3 = (mMeanV3 .* xPoints) + intMeanV3;
UpperSlopeV3   = mMeanV3 + mSTEV3;
LowerSlopeV3   = mMeanV3 - mSTEV3;
yPointsUpperV3 = (UpperSlopeV3 .* xPoints) + (intMeanV3 + intSTEV3);
yPointsLowerV3 = (LowerSlopeV3 .* xPoints) + (intMeanV3 - intSTEV3);
vV3 = [0 yPointsLowerV3(1); xPoints(end) yPointsLowerV3(end); xPoints(end) yPointsUpperV3(end); 0 yPointsUpperV3(1)];
patch('Vertices',vV3,'Faces',f,'FaceColor',colorV3,'FaceAlpha',0.4,'EdgeAlpha',0)
plot(xPoints,yPointsV3,'Color',colorV3,'LineWidth',2)
myTit = title('V3','fontsize',36); set(myTit,'interpreter','none');
set(myTit,'Position',[5 5.5])
set(gca,'fontsize',18);

%%%%%%%%% Plot V4 %%%%%%%%%%
subplot_tight(1,5,4,[0.08, 0.02]);
vV4 = [0 yPointsLowerV4(1); xPoints(end) yPointsLowerV4(end); xPoints(end) yPointsUpperV4(end); 0 yPointsUpperV4(1)];
patch('Vertices',vV4,'Faces',f,'FaceColor',colorV4,'FaceAlpha',0.4,'EdgeAlpha',0)
hold on; axis square; 
set(gca,'xlim',[0 7],'ylim',[0 7],'box','off','ycolor','w','tickdir','out','fontsize',14);
plot(xPoints,yPointsV4,'Color',colorV4,'LineWidth',2)
% now plot kids
mV4 = slope(kidI,7); 
intV4 = intercept(kidI,7); 
colorV4 = [153 255 51]; colorV4=colorV4./255;
mV4(isnan(mV4))=[]; intV4(isnan(intV4))=[];
mMeanV4 = mean(mV4); intMeanV4 = mean(intV4);
mSTEV4 = std(mV4)/sqrt(length(mV4));
intSTEV4 = std(intV4)/sqrt(length(intV4));
yPointsV4 = (mMeanV4 .* xPoints) + intMeanV4;
UpperSlopeV4   = mMeanV4 + mSTEV4;
LowerSlopeV4   = mMeanV4 - mSTEV4;
yPointsUpperV4 = (UpperSlopeV4 .* xPoints) + (intMeanV4 + intSTEV4);
yPointsLowerV4 = (LowerSlopeV4 .* xPoints) + (intMeanV4 - intSTEV4);
vV4 = [0 yPointsLowerV4(1); xPoints(end) yPointsLowerV4(end); xPoints(end) yPointsUpperV4(end); 0 yPointsUpperV4(1)];
patch('Vertices',vV4,'Faces',f,'FaceColor',colorV4,'FaceAlpha',0.4,'EdgeAlpha',0)
plot(xPoints,yPointsV4,'Color',colorV4,'LineWidth',2)
myTit = title('hV4','fontsize',36); set(myTit,'interpreter','none');
set(myTit,'Position',[5 5.5])
set(gca,'fontsize',18);

%%%%%%%%%%% Plot VO1 %%%%%%%%%%
subplot_tight(1,5,5,[0.08, 0.02]); 
vVO1 = [0 yPointsLowerVO1(1); xPoints(end) yPointsLowerVO1(end); xPoints(end) yPointsUpperVO1(end); 0 yPointsUpperVO1(1)];
patch('Vertices',vVO1,'Faces',f,'FaceColor',colorVO1,'FaceAlpha',0.4,'EdgeAlpha',0)
hold on; axis square;
set(gca,'xlim',[0 7],'ylim',[0 7],'box','off','ycolor','w','tickdir','out','fontsize',14);
plot(xPoints,yPointsVO1,'Color',colorVO1,'LineWidth',2)
% now plot kids
mVO1 = slope(kidI,9);
    intVO1 = intercept(kidI,9);
    colorVO1 = [102 255 255]; colorVO1=colorVO1./255;
mVO1(isnan(mVO1))=[]; intVO1(isnan(intVO1))=[];
mMeanVO1 = mean(mVO1); intMeanVO1 = mean(intVO1);
mSTEVO1 = std(mVO1)/sqrt(length(mVO1));
intSTEVO1 = std(intVO1)/sqrt(length(intVO1));
yPointsVO1 = (mMeanVO1 .* xPoints) + intMeanVO1;
UpperSlopeVO1   = mMeanVO1 + mSTEVO1;
LowerSlopeVO1   = mMeanVO1 - mSTEVO1;
yPointsUpperVO1 = (UpperSlopeVO1 .* xPoints) + (intMeanVO1 + intSTEVO1);
yPointsLowerVO1 = (LowerSlopeVO1 .* xPoints) + (intMeanVO1 - intSTEVO1);
vVO1 = [0 yPointsLowerVO1(1); xPoints(end) yPointsLowerVO1(end); xPoints(end) yPointsUpperVO1(end); 0 yPointsUpperVO1(1)];
patch('Vertices',vVO1,'Faces',f,'FaceColor',colorVO1,'FaceAlpha',0.4,'EdgeAlpha',0)
plot(xPoints,yPointsVO1,'Color',colorVO1,'LineWidth',2)
myTit = title('VO1','fontsize',36); set(myTit,'interpreter','none');
set(myTit,'Position',[5 5.5])
set(gca,'fontsize',18);

%% Now save
saveFile = fullfile(saveDir,'slopes_allSubsMatched_V1-VO1_oneRow.fig');
saveas(gcf,saveFile)



