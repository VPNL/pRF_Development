function y = nanste(varargin)
%Replacement for Matlab NANSTD Standard deviation, ignoring NaNs.
%

vals = varargin{1};
goodVals = ~isnan(vals);
y = sqrt(nanvar(varargin{:})) / (sqrt(length(goodVals)));

return;
