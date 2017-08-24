function figHandle = s_createCoveragePlot_averaged(RFcov, name)

% plotting subroutine for rmPlotCoverage. Broken off by ras 10/2009.
% And edited for scripting by JG 05/2016
% All you have to do is feed it a 128x128 coverage

% Load dummy filler data for your study. MAKE SURE THESE VARIABLE MATCH
% YOUR SPECIFIC STUDY
curdir=pwd;
fileDir = fullfile(curdir,'voxel_data');
load(fullfile(fileDir,'dummyVariables.mat'));


figHandle = gcf;

rfMax = max(RFcov(:)); 


% uncomment the line below if you want to normalize your averaged plot. As
% of now it will just plot the input field coverage matrix without
% normalizing. 
%img = RFcov ./ rfMax;
img = RFcov;
mask = makecircle(length(img));
img = img .* mask;
imagesc(data.X(1,:), data.Y(:,1), img);
set(gca, 'YDir', 'normal');
grid on

colormap(vfc.cmap);
colorbar;

% start plotting
hold on;


% add polar grid on top
p.ringTicks = (1:3)/3*vfc.fieldRange;
p.color = 'w';
polarPlot([], p);


% scale z-axis
% if vfc.normalizeRange
% 	if isequal( lower(vfc.method), 'maximum profile' )
% 		caxis([.3 1]);
% 	else
	    caxis([0 1]);
% 	end
% else
%     if min(RFcov(:))>=0
%         caxis([0 ceil(max(RFcov(:)))]);
%     else
%         caxis([-1 1] * ceil(max(abs(RFcov(:)))));
%     end
% end
axis image;   % axis square;
xlim([-vfc.fieldRange vfc.fieldRange])
ylim([-vfc.fieldRange vfc.fieldRange])

title(name, 'FontSize', 16, 'Interpreter', 'none');

return;